#pragma once

// Preferential-sampling grid-cell presence/absence sub-model.
//
// The catch/observation model elsewhere in sdmTMB.cpp is completely
// unaffected by this module: it only *reads* the already-computed fixed
// effects (`b_j`, model column 0 only -- not delta-model-aware) and
// spatial/spatiotemporal fields (`omega_s`, `epsilon_st`, model column 0 only),
// and re-projects them onto a user-supplied grid via `A_pref` to drive a
// Bernoulli "was this grid cell/row sampled" sub-model.
//
// The grid-level linear predictor could be extended to include lots of additional
// features - smooths, random effects, spatiotemporal fields, etc

namespace sdmTMB {

template <class Type>
struct preferential_data_t {
  int n_pref;                        // number of preferential_grid rows; 0 when the feature is off
  vector<Type> R_i;                  // 0/1 "was this row sampled" indicator; length n_pref
  matrix<Type> X_pref_ij;            // n_pref x ncol(X_ij[[1]]) fixed-effect design matrix
  Eigen::SparseMatrix<Type> A_pref;  // n_pref x n_s mesh-projection matrix
  vector<int> year_i_pref;           // 0-indexed year per row; length n_pref
  int b_pref_type;                   // 0 = constant, 1 = iid, 2 = rw

  preferential_data_t(SEXP x) {
    n_pref      = CppAD::Integer(asVector<Type>(getListElement(x, "n_pref"))[0]);
    R_i         = asVector<Type>(getListElement(x, "R_i"));
    X_pref_ij   = asMatrix<Type>(getListElement(x, "X_pref_ij"));
    A_pref      = tmbutils::asSparseMatrix<Type>(getListElement(x, "A_pref"));
    year_i_pref = asVector<int>(getListElement(x, "year_i_pref"));
    b_pref_type = CppAD::Integer(asVector<Type>(getListElement(x, "b_pref_type"))[0]);
  }
};

}  // namespace sdmTMB
