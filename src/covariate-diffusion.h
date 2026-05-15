#pragma once

namespace sdmTMB {

template <class Type>
struct covariate_diffusion_data_t {
  int n_terms;
  int n_covariates;
  array<Type> covariate_vertex_time;
  array<Type> proj_covariate_vertex_time;
  vector<int> term_component;
  vector<int> term_covariate;
  matrix<int> has; // [n_covariates x 3], cols = {space, time, spacetime}

  covariate_diffusion_data_t(SEXP x) {
    n_terms      = CppAD::Integer(asVector<Type>(getListElement(x, "n_terms"))[0]);
    n_covariates = CppAD::Integer(asVector<Type>(getListElement(x, "n_covariates"))[0]);
    covariate_vertex_time =
      tmbutils::asArray<Type>(getListElement(x, "covariate_vertex_time"));
    proj_covariate_vertex_time =
      tmbutils::asArray<Type>(getListElement(x, "proj_covariate_vertex_time"));
    term_component = asVector<int>(getListElement(x, "term_component"));
    term_covariate = asVector<int>(getListElement(x, "term_covariate"));
    has.resize(n_covariates, 3);
    has.setZero();
    for (int t = 0; t < n_terms; t++) {
      int cov_i = term_covariate(t);
      int comp = term_component(t);
      if (cov_i >= 0 && cov_i < n_covariates && comp >= 0 && comp < 3) {
        has(cov_i, comp) = 1;
      }
    }
  }
};

// Bundle of covariate-diffusion inputs; lag coefficients live at the tail of b_j.
template <class Type>
struct CovariateDiffusionContext {
  int n_terms;
  int n_covariates;
  int n_i;
  int n_t;
  const vector<int>& term_component;
  const vector<int>& term_covariate;
  array<Type>& covariate_vertex_time;
  const Eigen::SparseMatrix<Type>& A_st;
  const vector<int>& A_spatial_index;
  const vector<int>& year_i;
  const Eigen::SparseMatrix<Type>& M0;
  const Eigen::SparseMatrix<Type>& M1;
  const vector<Type>& kappaS_by_covariate;
  const vector<Type>& kappaT_by_covariate;
  const vector<Type>& kappaST_by_covariate;
  const vector<Type>& b_j;
  int model_col;
};

enum CovariateDiffusionComponent {
  dl_space = 0,
  dl_time = 1,
  dl_spacetime = 2
};

inline bool dl_is_valid_component(int component) {
  return component == dl_space ||
    component == dl_time ||
    component == dl_spacetime;
}

template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, 1> dl_get_covariate_col(
    array<Type>& covariate_vertex_time, int cov_i, int t, int n_vertices) {
  Eigen::Matrix<Type, Eigen::Dynamic, 1> col(n_vertices);
  for (int v = 0; v < n_vertices; v++) col(v) = covariate_vertex_time(v, t, cov_i);
  return col;
}

template <class Type>
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
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& transformed_vertex_time) {
  transformed_vertex_time.setZero();
  if (!dl_is_valid_component(component)) return false;

  if (component == dl_space) {
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

  if (component == dl_time) {
    Type denom = Type(1.0) + kappaT_dl;
    for (int v = 0; v < n_vertices; v++) {
      transformed_vertex_time(v, 0) = covariate_vertex_time(v, 0, cov_i) / denom;
    }
    for (int t = 1; t < n_t; t++) {
      for (int v = 0; v < n_vertices; v++) {
        transformed_vertex_time(v, t) =
          (covariate_vertex_time(v, t, cov_i) + kappaT_dl * transformed_vertex_time(v, t - 1)) / denom;
      }
    }
    return true;
  }

  if (component == dl_spacetime) {
    Eigen::SparseMatrix<Type> spacetime_system =
      M0_dl + (kappaS_scale - kappaST_scale) * M1_dl;
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu_spacetime;
    lu_spacetime.compute(spacetime_system);
    if (lu_spacetime.info() != Eigen::Success) return false;
    for (int t = 0; t < n_t; t++) {
      Eigen::Matrix<Type, Eigen::Dynamic, 1> rhs =
        M0_dl * dl_get_covariate_col(covariate_vertex_time, cov_i, t, n_vertices);
      if (t > 0) rhs -= kappaST_scale * (M1_dl * transformed_vertex_time.col(t - 1));
      Eigen::Matrix<Type, Eigen::Dynamic, 1> solved = lu_spacetime.solve(rhs);
      if (lu_spacetime.info() != Eigen::Success) return false;
      transformed_vertex_time.col(t) = solved;
    }
    return true;
  }

  return false;
}

template <class Type>
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

template <class Type>
void add_covariate_diffusion_to_eta_fixed(
    array<Type>& eta_fixed_i,
    CovariateDiffusionContext<Type>& ctx) {
  if (ctx.n_terms <= 0) return;
  if (ctx.n_covariates <= 0) {
    error("Distributed lag metadata error: n_covariates must be > 0 when n_terms > 0.");
  }
  if (ctx.kappaS_by_covariate.size() != ctx.n_covariates ||
      ctx.kappaT_by_covariate.size() != ctx.n_covariates ||
      ctx.kappaST_by_covariate.size() != ctx.n_covariates) {
    error("Distributed lag parameter length mismatch with n_covariates.");
  }

  int n_vertices_dl = ctx.covariate_vertex_time.dim[0];
  int n_t_dl = ctx.covariate_vertex_time.dim[1];

  // Determine required solvers/scales by scanning terms
  std::vector<int> cov_needs_spatial_scale(ctx.n_covariates, 0);
  std::vector<int> cov_needs_spatial_solver(ctx.n_covariates, 0);
  for (int term = 0; term < ctx.n_terms; term++) {
    int component = ctx.term_component(term);
    int cov_i = ctx.term_covariate(term);
    if (!dl_is_valid_component(component)) {
      error("Distributed lag metadata error: invalid component code (expected spatial=0, temporal=1, spatiotemporal=2).");
    }
    if (cov_i < 0 || cov_i >= ctx.n_covariates) {
      error("Distributed lag metadata error: term covariate index out of bounds.");
    }
    if (component == dl_space) {
      cov_needs_spatial_scale[cov_i] = 1;
      cov_needs_spatial_solver[cov_i] = 1;
    }
    if (component == dl_spacetime) {
      cov_needs_spatial_scale[cov_i] = 1;
    }
  }

  // Compute per-covariate derived scales
  vector<Type> kappaS_scale(ctx.n_covariates);
  vector<Type> kappaST_scale(ctx.n_covariates);
  kappaS_scale.setZero();
  kappaST_scale.setZero();
  for (int cov_i = 0; cov_i < ctx.n_covariates; cov_i++) {
    if (cov_needs_spatial_scale[cov_i]) {
      kappaS_scale(cov_i) = Type(1.0) / (ctx.kappaS_by_covariate(cov_i) * ctx.kappaS_by_covariate(cov_i));
      kappaST_scale(cov_i) = ctx.kappaST_by_covariate(cov_i) * kappaS_scale(cov_i);
    }
  }

  // Factorize per-covariate spatial systems once
  std::vector< Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > >
    lu_spatial_by_covariate(ctx.n_covariates);
  for (int cov_i = 0; cov_i < ctx.n_covariates; cov_i++) {
    if (!cov_needs_spatial_solver[cov_i]) continue;
    Eigen::SparseMatrix<Type> spatial_system = ctx.M0 + kappaS_scale(cov_i) * ctx.M1;
    lu_spatial_by_covariate[cov_i].compute(spatial_system);
    if (lu_spatial_by_covariate[cov_i].info() != Eigen::Success) {
      error("Distributed lag sparse solve failed while factorizing spatial system (M0 + kappa^{-2} M1).");
    }
  }

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
      kappaS_scale(cov_i),
      ctx.kappaT_by_covariate(cov_i),
      kappaST_scale(cov_i),
      cov_needs_spatial_solver[cov_i] == 1,
      lu_spatial_by_covariate[cov_i],
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
