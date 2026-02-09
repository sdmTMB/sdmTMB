#pragma once

namespace sdmTMB {

// Bundles all distributed-lag inputs passed from R/TMB DATA and PARAMETER
// blocks so that helper functions take a single context reference instead of
// 20+ individual arguments.
template <class Type>
struct DistributedLagContext {
  // Dimensions
  int n_terms;
  int n_covariates;
  int n_i;
  int n_t;

  // Per-term metadata (length n_terms, 0-based)
  const vector<int>& term_component;
  const vector<int>& term_covariate;

  // Covariate field on mesh vertices x time x covariates
  array<Type>& covariate_vertex_time;

  // Mesh projection + observation mapping
  const Eigen::SparseMatrix<Type>& A_st;
  const vector<int>& A_spatial_index;
  const vector<int>& year_i;

  // SPDE mass/stiffness matrices (from main mesh)
  const Eigen::SparseMatrix<Type>& M0;
  const Eigen::SparseMatrix<Type>& M1;

  // Transformed lag scale parameters (already on natural scale)
  Type kappaS;
  Type kappaT;
  Type kappaST;

  // Fixed-effect coefficients (full vector; lag slots are at the tail)
  const vector<Type>& b_j;

  // Column of eta_fixed_i/proj_fe to accumulate into (0-based model index)
  int model_col;
};

enum DistributedLagComponent {
  dl_space = 0,
  dl_time = 1,
  dl_spacetime = 2
};

inline bool dl_is_valid_component(int component) {
  return component == dl_space ||
    component == dl_time ||
    component == dl_spacetime;
}

// Extract column t of covariate cov_i from 3D array (vertices x time x covariates)
template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, 1> dl_get_covariate_col(
    array<Type>& covariate_vertex_time, int cov_i, int t, int n_vertices) {
  Eigen::Matrix<Type, Eigen::Dynamic, 1> col(n_vertices);
  for (int v = 0; v < n_vertices; v++) col(v) = covariate_vertex_time(v, t, cov_i);
  return col;
}

template <class Type>
// Solve the per-term distributed-lag transport on vertices x time.
// Uses sparse solves for spatial/spatiotemporal components and a
// recursive update for temporal components. Returns false on invalid
// component codes or sparse solver failures.
bool dl_solve_transformed_vertex_time(
    int component,
    array<Type>& covariate_vertex_time,
    int cov_i,
    int n_vertices,
    int n_t,
    const Eigen::SparseMatrix<Type>& M0_dl,
    const Eigen::SparseMatrix<Type>& M1_dl,
    Type kappaS_scale,
    Type kappaT_dl,
    Type kappaST_scale,
    bool has_spatial_solver,
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> >& lu_spatial,
    bool has_m0_solver,
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> >& lu_m0,
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& transformed_vertex_time) {
  transformed_vertex_time.setZero();
  if (!dl_is_valid_component(component)) return false;

  if (component == dl_space) { // spatial
    if (!has_spatial_solver || lu_spatial.info() != Eigen::Success) return false;
    for (int t = 0; t < n_t; t++) {
      Eigen::Matrix<Type, Eigen::Dynamic, 1> rhs =
        M0_dl * dl_get_covariate_col(covariate_vertex_time, cov_i, t, n_vertices);
      Eigen::Matrix<Type, Eigen::Dynamic, 1> solved = lu_spatial.solve(rhs);
      if (lu_spatial.info() != Eigen::Success) return false;
      transformed_vertex_time.col(t) = solved;
    }
    return true;
  }

  if (component == dl_time) { // temporal
    for (int v = 0; v < n_vertices; v++) {
      transformed_vertex_time(v, 0) = covariate_vertex_time(v, 0, cov_i);
    }
    for (int t = 1; t < n_t; t++) {
      for (int v = 0; v < n_vertices; v++) {
        transformed_vertex_time(v, t) =
          covariate_vertex_time(v, t, cov_i) + kappaT_dl * transformed_vertex_time(v, t - 1);
      }
    }
    return true;
  }

  if (component == dl_spacetime) { // spatiotemporal
    if (!has_m0_solver || lu_m0.info() != Eigen::Success) return false;
    for (int v = 0; v < n_vertices; v++) {
      transformed_vertex_time(v, 0) = covariate_vertex_time(v, 0, cov_i);
    }
    for (int t = 1; t < n_t; t++) {
      Eigen::Matrix<Type, Eigen::Dynamic, 1> rhs =
        M0_dl * dl_get_covariate_col(covariate_vertex_time, cov_i, t, n_vertices) -
        kappaST_scale * (M1_dl * transformed_vertex_time.col(t - 1));
      Eigen::Matrix<Type, Eigen::Dynamic, 1> solved = lu_m0.solve(rhs);
      if (lu_m0.info() != Eigen::Success) return false;
      transformed_vertex_time.col(t) = solved;
    }
    return true;
  }

  return false;
}

template <class Type>
// Project transformed lag fields from vertices to observations.
// First projects each time slice through A_st to unique spatial rows,
// then indexes those projections to observation rows via year_i and
// A_spatial_index.
Eigen::Matrix<Type, Eigen::Dynamic, 1> dl_project_vertex_time_to_observations(
    const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& transformed_vertex_time,
    const Eigen::SparseMatrix<Type>& A_st,
    const vector<int>& A_spatial_index,
    const vector<int>& year_i,
    int n_i,
    int n_t) {
  std::vector<Eigen::Matrix<Type, Eigen::Dynamic, 1> > projected_by_t;
  projected_by_t.reserve(n_t);
  for (int t = 0; t < n_t; t++) {
    projected_by_t.push_back(A_st * transformed_vertex_time.col(t));
  }
  Eigen::Matrix<Type, Eigen::Dynamic, 1> term_i(n_i);
  term_i.setZero();
  for (int i = 0; i < n_i; i++) {
    term_i(i) = projected_by_t[year_i(i)](A_spatial_index(i));
  }
  return term_i;
}

// Add distributed-lag term contributions to eta_fixed_i.
// All dimension/index consistency is enforced by the R layer;
// C++ only checks conditions that can arise at solve time
// (sparse factorization failures, parameter presence).
template <class Type>
void add_distributed_lags_to_eta_fixed(
    array<Type>& eta_fixed_i,
    DistributedLagContext<Type>& ctx) {
  if (ctx.n_terms <= 0) return;

  int n_vertices_dl = ctx.covariate_vertex_time.dim[0];
  int n_t_dl = ctx.covariate_vertex_time.dim[1];

  // Determine which scale parameters are needed by scanning components
  bool needs_kappaS = false;
  bool needs_kappaST = false;
  bool need_spatial_solver = false;
  bool need_m0_solver = false;
  for (int term = 0; term < ctx.n_terms; term++) {
    int component = ctx.term_component(term);
    if (!dl_is_valid_component(component)) {
      error("Distributed lag metadata error: invalid component code (expected spatial=0, temporal=1, spatiotemporal=2).");
    }
    if (component == dl_space) { needs_kappaS = true; need_spatial_solver = true; }
    if (component == dl_spacetime) {
      needs_kappaS = true;
      needs_kappaST = true;
      need_m0_solver = true;
    }
  }

  // Compute derived scales
  Type kappaS_scale = Type(0.0);
  Type kappaST_scale = Type(0.0);
  if (needs_kappaS) kappaS_scale = Type(1.0) / (ctx.kappaS * ctx.kappaS);
  if (needs_kappaST) kappaST_scale = ctx.kappaST * kappaS_scale;

  // Factorize spatial system once for all spatial terms
  Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu_spatial;
  if (need_spatial_solver) {
    Eigen::SparseMatrix<Type> spatial_system = ctx.M0 + kappaS_scale * ctx.M1;
    lu_spatial.compute(spatial_system);
    if (lu_spatial.info() != Eigen::Success) {
      error("Distributed lag sparse solve failed while factorizing spatial system (M0 + kappa^{-2} M1).");
    }
  }

  // Factorize M0 once for spatiotemporal terms
  Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu_m0;
  if (need_m0_solver) {
    lu_m0.compute(ctx.M0);
    if (lu_m0.info() != Eigen::Success) {
      error("Distributed lag sparse solve failed while factorizing M0.");
    }
  }

  // Solve and project each term, accumulate into eta_fixed_i
  int dl_coef_start = ctx.b_j.size() - ctx.n_terms;
  for (int term = 0; term < ctx.n_terms; term++) {
    int component = ctx.term_component(term);
    int cov_i = ctx.term_covariate(term);

    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> transformed_vertex_time(n_vertices_dl, n_t_dl);
    bool solved = dl_solve_transformed_vertex_time(
      component,
      ctx.covariate_vertex_time,
      cov_i,
      n_vertices_dl,
      n_t_dl,
      ctx.M0,
      ctx.M1,
      kappaS_scale,
      ctx.kappaT,
      kappaST_scale,
      need_spatial_solver,
      lu_spatial,
      need_m0_solver,
      lu_m0,
      transformed_vertex_time
    );
    if (!solved) {
      error("Distributed lag sparse solve failed while transforming a lagged covariate.");
    }

    Eigen::Matrix<Type, Eigen::Dynamic, 1> term_i = dl_project_vertex_time_to_observations(
      transformed_vertex_time,
      ctx.A_st,
      ctx.A_spatial_index,
      ctx.year_i,
      ctx.n_i,
      ctx.n_t
    );

    Type beta_dl = ctx.b_j(dl_coef_start + term);
    for (int i = 0; i < ctx.n_i; i++) eta_fixed_i(i, ctx.model_col) += beta_dl * term_i(i);
  }
}

}  // namespace sdmTMB
