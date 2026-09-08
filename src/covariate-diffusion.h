#pragma once

namespace sdmTMB {

template <class Type>
struct covariate_diffusion_data_t {
  int n_terms;
  int n_covariates;
  tmbutils::array<Type> covariate_vertex_time;
  tmbutils::array<Type> proj_covariate_vertex_time;
  vector<int> term_component;
  vector<int> term_covariate;
  matrix<int> has; // [n_covariates x 2], cols = {space, time}

  covariate_diffusion_data_t(SEXP x) {
    n_terms      = CppAD::Integer(asVector<Type>(getListElement(x, "n_terms"))[0]);
    n_covariates = CppAD::Integer(asVector<Type>(getListElement(x, "n_covariates"))[0]);
    covariate_vertex_time =
      tmbutils::asArray<Type>(getListElement(x, "covariate_vertex_time"));
    proj_covariate_vertex_time =
      tmbutils::asArray<Type>(getListElement(x, "proj_covariate_vertex_time"));
    term_component = asVector<int>(getListElement(x, "term_component"));
    term_covariate = asVector<int>(getListElement(x, "term_covariate"));
    has.resize(n_covariates, 2);
    has.setZero();
    for (int t = 0; t < n_terms; t++) {
      int cov_i = term_covariate(t);
      int comp = term_component(t);
      if (cov_i >= 0 && cov_i < n_covariates) {
        if (comp == 0 || comp == 2) has(cov_i, 0) = 1;
        if (comp == 1 || comp == 2) has(cov_i, 1) = 1;
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
  tmbutils::array<Type>& covariate_vertex_time;
  const Eigen::SparseMatrix<Type>& A_st;
  const vector<int>& A_spatial_index;
  const vector<int>& year_i;
  const Eigen::SparseMatrix<Type>& M0;
  const Eigen::SparseMatrix<Type>& M1;
  const vector<Type>& kappaS_by_covariate;
  const vector<Type>& kappaT_by_covariate;
  const vector<Type>& b_j;
  int model_col;
};

enum CovariateDiffusionComponent {
  nl_space = 0,
  nl_time = 1,
  nl_joint = 2
};

inline bool nl_is_valid_component(int component) {
  return component == nl_space ||
    component == nl_time ||
    component == nl_joint;
}

template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, 1> nl_get_covariate_col(
    tmbutils::array<Type>& covariate_vertex_time, int cov_i, int t, int n_vertices) {
  Eigen::Matrix<Type, Eigen::Dynamic, 1> col(n_vertices);
  for (int v = 0; v < n_vertices; v++) col(v) = covariate_vertex_time(v, t, cov_i);
  return col;
}

template <class Type>
bool nl_solve_transformed_vertex_time(
    int component,
    tmbutils::array<Type>& covariate_vertex_time,
    int cov_i,
    int n_vertices,
    int n_t,
    const Eigen::SparseMatrix<Type>& M0_nl,
    const Eigen::SparseMatrix<Type>& M1_nl,
    Type kappaS_scale,
    Type kappaT_nl,
    bool has_system_solver,
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> >& lu_system,
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& transformed_vertex_time) {
  transformed_vertex_time.setZero();
  if (!nl_is_valid_component(component)) return false;

  if (component == nl_space) {
    if (!has_system_solver || lu_system.info() != Eigen::Success) return false;
    for (int t = 0; t < n_t; t++) {
      Eigen::Matrix<Type, Eigen::Dynamic, 1> rhs =
        M0_nl * nl_get_covariate_col(covariate_vertex_time, cov_i, t, n_vertices);
      Eigen::Matrix<Type, Eigen::Dynamic, 1> solved = lu_system.solve(rhs);
      if (lu_system.info() != Eigen::Success) return false;
      transformed_vertex_time.col(t) = solved;
    }
    return true;
  }

  if (component == nl_time) {
    Type denom = Type(1.0) + kappaT_nl;
    for (int v = 0; v < n_vertices; v++) {
      transformed_vertex_time(v, 0) = covariate_vertex_time(v, 0, cov_i) / denom;
    }
    for (int t = 1; t < n_t; t++) {
      for (int v = 0; v < n_vertices; v++) {
        transformed_vertex_time(v, t) =
          (covariate_vertex_time(v, t, cov_i) + kappaT_nl * transformed_vertex_time(v, t - 1)) / denom;
      }
    }
    return true;
  }

  if (component == nl_joint) {
    if (!has_system_solver || lu_system.info() != Eigen::Success) return false;
    for (int t = 0; t < n_t; t++) {
      Eigen::Matrix<Type, Eigen::Dynamic, 1> rhs =
        M0_nl * nl_get_covariate_col(covariate_vertex_time, cov_i, t, n_vertices);
      if (t > 0) rhs += kappaT_nl * M0_nl * transformed_vertex_time.col(t - 1);
      Eigen::Matrix<Type, Eigen::Dynamic, 1> solved = lu_system.solve(rhs);
      if (lu_system.info() != Eigen::Success) return false;
      transformed_vertex_time.col(t) = solved;
    }
    return true;
  }

  return false;
}

template <class Type>
Eigen::Matrix<Type, Eigen::Dynamic, 1> nl_project_vertex_time_to_observations(
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
    tmbutils::array<Type>& eta_fixed_i,
    CovariateDiffusionContext<Type>& ctx,
    matrix<Type>* term_values_out = nullptr) {
  if (ctx.n_terms <= 0) return;
  if (ctx.n_covariates <= 0) {
    error("Nonlocal metadata error: n_covariates must be > 0 when n_terms > 0.");
  }
  if (ctx.kappaS_by_covariate.size() != ctx.n_covariates ||
      ctx.kappaT_by_covariate.size() != ctx.n_covariates) {
    error("Nonlocal parameter length mismatch with n_covariates.");
  }

  int n_vertices_nl = ctx.covariate_vertex_time.dim[0];
  int n_t_nl = ctx.covariate_vertex_time.dim[1];

  // Determine required solvers/scales by scanning terms
  std::vector<int> cov_needs_spatial_scale(ctx.n_covariates, 0);
  std::vector<int> cov_needs_system_solver(ctx.n_covariates, 0);
  std::vector<int> cov_uses_joint_system(ctx.n_covariates, 0);
  for (int term = 0; term < ctx.n_terms; term++) {
    int component = ctx.term_component(term);
    int cov_i = ctx.term_covariate(term);
    if (!nl_is_valid_component(component)) {
      error("Nonlocal metadata error: invalid component code (expected spatial=0, temporal=1, or joint=2).");
    }
    if (cov_i < 0 || cov_i >= ctx.n_covariates) {
      error("Nonlocal metadata error: term covariate index out of bounds.");
    }
    if (component == nl_space || component == nl_joint) {
      cov_needs_spatial_scale[cov_i] = 1;
      cov_needs_system_solver[cov_i] = 1;
      if (component == nl_joint) cov_uses_joint_system[cov_i] = 1;
    }
  }

  // Compute per-covariate derived scales
  vector<Type> kappaS_scale(ctx.n_covariates);
  kappaS_scale.setZero();
  for (int cov_i = 0; cov_i < ctx.n_covariates; cov_i++) {
    if (cov_needs_spatial_scale[cov_i]) {
      kappaS_scale(cov_i) = Type(1.0) / (ctx.kappaS_by_covariate(cov_i) * ctx.kappaS_by_covariate(cov_i));
    }
  }

  // Factorize each spatial or joint system once per covariate
  std::vector< Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > >
    lu_system_by_covariate(ctx.n_covariates);
  for (int cov_i = 0; cov_i < ctx.n_covariates; cov_i++) {
    if (!cov_needs_system_solver[cov_i]) continue;
    Type temporal_scale = cov_uses_joint_system[cov_i] == 1 ?
      Type(1.0) + ctx.kappaT_by_covariate(cov_i) : Type(1.0);
    Eigen::SparseMatrix<Type> system = temporal_scale * ctx.M0 + kappaS_scale(cov_i) * ctx.M1;
    lu_system_by_covariate[cov_i].compute(system);
    if (lu_system_by_covariate[cov_i].info() != Eigen::Success) {
      error("Nonlocal sparse solve failed while factorizing a spatial or joint system.");
    }
  }

  int nl_coef_start = ctx.b_j.size() - ctx.n_terms;
  for (int term = 0; term < ctx.n_terms; term++) {
    int component = ctx.term_component(term);
    int cov_i = ctx.term_covariate(term);

    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> transformed_vertex_time(n_vertices_nl, n_t_nl);
    bool solved = nl_solve_transformed_vertex_time(
      component,
      ctx.covariate_vertex_time,
      cov_i,
      n_vertices_nl,
      n_t_nl,
      ctx.M0,
      ctx.M1,
      kappaS_scale(cov_i),
      ctx.kappaT_by_covariate(cov_i),
      cov_needs_system_solver[cov_i] == 1,
      lu_system_by_covariate[cov_i],
      transformed_vertex_time
    );
    if (!solved) {
      error("Nonlocal sparse solve failed while transforming a lagged covariate.");
    }

    Eigen::Matrix<Type, Eigen::Dynamic, 1> term_i = nl_project_vertex_time_to_observations(
      transformed_vertex_time,
      ctx.A_st,
      ctx.A_spatial_index,
      ctx.year_i,
      ctx.n_i,
      ctx.n_t
    );

    if (term_values_out != nullptr) {
      term_values_out->col(term) = term_i;
    }

    Type beta_nl = ctx.b_j(nl_coef_start + term);
    for (int i = 0; i < ctx.n_i; i++) eta_fixed_i(i, ctx.model_col) += beta_nl * term_i(i);
  }
}

}  // namespace sdmTMB
