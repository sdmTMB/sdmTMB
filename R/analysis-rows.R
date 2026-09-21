# Keep the observation row contract explicit without allowing component-specific
# missing predictors to create different implicit model-frame row sets.

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

.establish_analysis_rows <- function(data, formulas, family_spec) {
  formulas <- .formula_list(formulas)
  if (length(formulas) < family_spec$n_m) {
    cli_abort("Internal row-identity error: not enough component formulas were supplied.")
  }
  row_family_id <- family_spec$family_id_i
  active <- .family_spec_component_active(family_spec, row_family_id)
  response_complete <- rep(TRUE, nrow(data))

  for (component in seq_len(family_spec$n_m)) {
    formula_no_smooths <- remove_s_and_t2(formulas[[component]])
    mf <- stats::model.frame(
      formula_no_smooths,
      data = data,
      na.action = stats::na.pass
    )
    component_response <- stats::model.response(mf, type = "any")
    component_complete <- stats::complete.cases(component_response)
    response_complete <- response_complete &
      (!active[, component] | component_complete)
  }

  used <- which(response_complete)
  if (!length(used)) {
    cli_abort("No rows with a non-missing response remain for fitting.")
  }
  omitted <- setdiff(seq_len(nrow(data)), used)
  original_to_analysis <- rep.int(NA_integer_, nrow(data))
  original_to_analysis[used] <- seq_along(used)
  list(
    original_n = as.integer(nrow(data)),
    used = as.integer(used),
    omitted = as.integer(omitted),
    original_to_analysis = as.integer(original_to_analysis)
  )
}

.subset_analysis_rows <- function(x, analysis_rows, name = "value") {
  if (is.null(x)) return(x)
  n <- if (is.matrix(x) || is.data.frame(x)) nrow(x) else length(x)
  if (n == analysis_rows$original_n) {
    if (is.matrix(x) || is.data.frame(x)) {
      return(x[analysis_rows$used, , drop = FALSE])
    }
    return(x[analysis_rows$used])
  }
  if (n == length(analysis_rows$used)) return(x)
  cli_abort(
    "`{name}` has {n} rows/elements, but expected {analysis_rows$original_n}."
  )
}
