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

.check_no_missing_covariates <- function(
    data,
    formulas,
    shared_formulas = list(),
    required_columns = character(0),
    stage = "data") {
  all_formulas <- c(.formula_list(formulas), .formula_list(shared_formulas))
  for (formula in all_formulas) {
    # raw variables, as in fitting; evaluating terms such as poly() or
    # scale() on `newdata` would recompute them from the prediction rows;
    # the response is not a covariate, so drop it
    if (length(formula) == 3L) formula <- formula[-2]
    complete <- .complete_model_rows(data, list(formula))
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
