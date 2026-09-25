# Rows with missing values in variables the model uses are omitted before
# fitting (see `sdmTMB()`); prediction data must be complete instead.

.formula_list <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "formula")) return(list(x))
  if (is.list(x) && all(vapply(x, inherits, logical(1), what = "formula"))) {
    return(x)
  }
  cli_abort("Expected a formula or a list of formulas.")
}

.covariate_model_frame <- function(formula, data) {
  formula_no_smooths <- remove_s_and_t2(formula)
  formula_no_response <- stats::formula(
    stats::delete.response(
      stats::terms(reformulas::subbars(formula_no_smooths))
    )
  )
  mf <- stats::model.frame(
    formula_no_response,
    data = data,
    na.action = stats::na.pass
  )

  # Smooth terms are removed before constructing the model frame because their
  # values are not ordinary model-frame columns. Check their data variables
  # directly so smooth covariates follow the same strict rule.
  formula_no_response_full <- stats::formula(
    stats::delete.response(stats::terms(formula))
  )
  smooth_vars <- setdiff(
    all.vars(formula_no_response_full),
    all.vars(stats::formula(stats::delete.response(stats::terms(formula_no_smooths))))
  )
  smooth_vars <- intersect(smooth_vars, names(data))
  complete <- stats::complete.cases(mf)
  for (variable in smooth_vars) {
    complete <- complete & stats::complete.cases(data[[variable]])
  }

  complete
}

.check_no_missing_covariates <- function(
    data,
    formulas,
    shared_formulas = list(),
    required_columns = character(0),
    stage = "data") {
  all_formulas <- c(.formula_list(formulas), .formula_list(shared_formulas))
  for (formula in all_formulas) {
    complete <- .covariate_model_frame(formula, data)
    if (any(!complete)) {
      rows <- which(!complete)
      cli_abort(c(
        "Missing covariate values are not supported in {stage}.",
        "x" = "Rows with missing covariates: {paste(rows, collapse = ', ')}.",
        "i" = "Remove those rows or impute the predictors before fitting or predicting."
      ))
    }
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    cli_abort(
      "Required column(s) missing from {stage}: {paste(missing_columns, collapse = ', ')}."
    )
  }
  for (column in unique(required_columns)) {
    if (anyNA(data[[column]])) {
      cli_abort(c(
        "Missing covariate values are not supported in {stage}.",
        "x" = "Column `{column}` has missing values.",
        "i" = "Remove those rows or impute the predictors before fitting or predicting."
      ))
    }
  }
  invisible(NULL)
}

.complete_model_rows <- function(data, formulas, row_vectors = list()) {
  vars <- intersect(unique(unlist(lapply(formulas, all.vars))), names(data))
  row_vectors <- Filter(function(x) NROW(x) == nrow(data), row_vectors)
  do.call(stats::complete.cases, c(list(data[vars]), row_vectors))
}

.subset_rows <- function(x, rows, n) {
  if (NROW(x) != n) return(x)
  if (is.matrix(x) || is.data.frame(x)) x[rows, , drop = FALSE] else x[rows]
}

# Mesh rows are the data rows (`make_mesh()`); areal domains are built from
# `data` in `sdmTMB()` and need no subsetting.
.subset_mesh_rows <- function(mesh, rows, n) {
  if (is_areal_domain(mesh) || NROW(mesh$loc_xy) != n) return(mesh)
  mesh$loc_xy <- mesh$loc_xy[rows, , drop = FALSE]
  mesh$A_st <- mesh$A_st[rows, , drop = FALSE]
  mesh$sdm_spatial_id <- seq_along(rows)
  mesh
}
