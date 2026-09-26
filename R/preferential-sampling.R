# Preferential sampling: a Bernoulli model for which cells of an eligible
# sampling frame were visited in each time step, linked to a standardized
# expected-catch surface from the main model. The surface is the main
# model's fitted design evaluated on the sampling frame, like prediction
# `newdata`, so the frame's covariate values define the standardization.
#
# Only the RTMB backend will support it. The specification and frame/design
# preparation live here; the joint likelihood is not implemented yet, so
# `rtmb_validate()` still rejects preferential data.

#' Preferential-sampling specification
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Specify a model for which cells of an eligible sampling frame were sampled
#' in each time step, linked to the main model's standardized expected catch.
#' Pass the result to the `preferential` argument of [sdmTMB()].
#'
#' **The joint likelihood is not implemented yet.** This constructor and the
#' preparation of its inputs are in place, but [sdmTMB()] currently stops
#' before fitting a model with a `preferential` specification.
#'
#' @param formula A two-sided formula for the sampling model, e.g.
#'   `sampled ~ 0 + factor(year) + distance_to_port`. The response column
#'   must be 0 (eligible but not sampled), 1 (sampled), or `NA` (unknown).
#'   These coefficients belong to the sampling model alone. Only ordinary
#'   fixed effects are supported: no smoothers, random effects, or `offset()`
#'   terms.
#' @param data A data frame with one row per eligible cell and time step. It
#'   must contain the mesh coordinates, the `time` column of the main model,
#'   the sampling response and covariates, and every covariate in the main
#'   model's fixed-effect formula. The main model is evaluated on these rows
#'   as with `newdata` in [predict.sdmTMB()], so set catchability covariates
#'   (e.g., gear) to reference values here.
#' @param re_form_iid `NA` (default) excludes the main model's IID random
#'   effects from the expected-catch surface. `NULL` would include them, but
#'   this is not supported yet for models with IID random effects.
#' @param offset Offset for the main model's expected-catch surface on the
#'   link scale (not an offset for the sampling model). A single value or one
#'   value per row of `data`. The default `0` means unit exposure.
#'   Observation offsets are not copied.
#' @param coefficient How the preference coefficient varies. Only
#'   `"constant"` (one coefficient) is currently supported.
#' @param spatial Whether to add a time-invariant spatial random field to the
#'   sampling model (`"off"` or `"on"`). This does not affect the main
#'   model's fields.
#'
#' @return A list of class `sdmTMB_preferential`.
#' @export
#' @examples
#' grid <- data.frame(
#'   X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1), year = 2020,
#'   sampled = c(1, 0, NA, 1), depth = 100
#' )
#' preferential_sampling(sampled ~ 1, data = grid)
preferential_sampling <- function(formula, data, re_form_iid = NA, offset = 0,
                                  coefficient = "constant",
                                  spatial = c("off", "on")) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    cli_abort("`formula` must be a two-sided formula such as `sampled ~ 1`.")
  }
  if (!is.name(formula[[2L]])) {
    cli_abort("The left side of `formula` must be the name of the sampling response column.")
  }
  rhs <- stats::delete.response(stats::terms(formula))
  if (length(reformulas::findbars(formula))) {
    cli_abort("Random effects are not supported in the sampling `formula` yet.")
  }
  if (length(get_smooth_terms(all_terms(rhs)))) {
    cli_abort("Smoothers are not supported in the sampling `formula` yet.")
  }
  if (!is.null(attr(rhs, "offset"))) {
    cli_abort(c(
      "`offset()` terms are not supported in the sampling `formula`.",
      "i" = "The `offset` argument is an offset for the main model's expected-catch surface, not for the sampling model."
    ))
  }
  if (!inherits(data, "data.frame")) {
    cli_abort("`data` must be a data frame.")
  }
  data <- as.data.frame(data)
  if (!nrow(data)) cli_abort("`data` has no rows.")
  response <- as.character(formula[[2L]])
  if (!response %in% names(data)) {
    cli_abort("`data` is missing the sampling response column {.field {response}}.")
  }
  r <- data[[response]]
  if (!(is.logical(r) || is.numeric(r)) || !all(r %in% c(0, 1, NA))) {
    cli_abort(c(
      "The sampling response {.field {response}} must be 0, 1, or `NA` (or logical).",
      "i" = "1 = sampled, 0 = eligible but not sampled, `NA` = unknown."
    ))
  }
  if (!(is.null(re_form_iid) || identical(re_form_iid, NA))) {
    cli_abort("`re_form_iid` must be `NA` (exclude IID random effects) or `NULL` (include them).")
  }
  if (!is.numeric(offset) || !all(is.finite(offset)) ||
      !length(offset) %in% c(1L, nrow(data))) {
    cli_abort("`offset` must be finite, with length 1 or `nrow(data)`.")
  }
  if (!identical(coefficient, "constant")) {
    cli_abort("Only `coefficient = \"constant\"` is currently supported.")
  }
  spatial <- match.arg(spatial)
  structure(
    list(
      formula = formula, data = data, response = response,
      include_iid = is.null(re_form_iid),
      offset = rep_len(as.numeric(offset), nrow(data)),
      coefficient = coefficient, spatial = spatial
    ),
    class = "sdmTMB_preferential"
  )
}

#' Reject model features the preferential-sampling likelihood doesn't support
#'
#' Called before any mesh projection or objective construction so that an
#' unsupported combination errors instead of silently fitting a model whose
#' shared catch surface omits part of the main model. Relax each guard only
#' once the corresponding support is implemented and tested.
#' @noRd
.validate_preferential_scope <- function(spec, formula, delta, multi_family,
                                         areal, mesh, mesh_missing,
                                         anisotropy, time_varying,
                                         spatial_varying, nonlocal_formula,
                                         normalize, backend, no_spatial) {
  if (!inherits(spec, "sdmTMB_preferential")) {
    cli_abort("`preferential` must be created with `preferential_sampling()`.")
  }
  if (!identical(backend, "rtmb")) {
    cli_abort("Preferential sampling requires backend = \"rtmb\"; set control = sdmTMBcontrol(backend = \"rtmb\").")
  }
  if (mesh_missing) {
    cli_abort("`mesh` must be supplied when using `preferential`.")
  }
  formulas <- .formula_list(formula)
  term_labels <- unlist(lapply(formulas, all_terms))
  has_iid <- any(vapply(formulas, function(f) {
    length(reformulas::findbars(f)) > 0L
  }, logical(1L)))
  unsupported <- c(
    "delta models" = delta,
    "multi-family models" = multi_family,
    "areal (SAR/CAR) models" = areal,
    "barrier meshes" = "spde_barrier" %in% names(mesh),
    "anisotropy" = isTRUE(anisotropy),
    "`time_varying`" = !is.null(time_varying),
    "`spatial_varying`" = !is.null(spatial_varying),
    "`nonlocal_formula`" = !is.null(nonlocal_formula),
    "threshold (`breakpt()`/`logistic()`) terms" =
      any(grepl("^(breakpt|logistic)\\(", term_labels)),
    "smoothers in `formula`" = length(get_smooth_terms(term_labels)) > 0L,
    "including IID random effects (`re_form_iid = NULL`)" =
      has_iid && spec$include_iid,
    "`normalize = TRUE`" = isTRUE(normalize)
  )
  if (any(unsupported)) {
    cli_abort(c(
      "Preferential sampling does not yet support some requested model features.",
      "x" = "Unsupported: {names(unsupported)[unsupported]}."
    ))
  }
  if (no_spatial) {
    cli_abort(c(
      "Preferential sampling requires a spatial or spatiotemporal field in the main model.",
      "i" = "The shared field is what identifies the preference coefficient separately from the sampling covariates."
    ))
  }
  invisible(NULL)
}

# Functions that summarize the data they are given. In a formula they are
# recomputed on new data unless the fitted values are kept in `predvars`, as
# `poly()` and spline bases do.
.data_dependent_calls <- c("scale", "mean", "sd", "var", "min", "max",
  "median", "range", "quantile")

.check_prediction_safe_terms <- function(terms) {
  predvars <- attr(terms, "predvars")
  if (is.null(predvars)) predvars <- attr(terms, "variables")
  calls <- all.names(predvars)
  bad <- intersect(calls, .data_dependent_calls)
  if (length(bad)) {
    cli_abort(c(
      "The main model formula uses {.fn {bad}}, which would be recomputed on the sampling `data`.",
      "i" = "Precompute the transformed covariate as a column in both `data` and the sampling `data`."
    ))
  }
}

.check_frame_columns <- function(data, vars, what) {
  missing <- setdiff(vars, names(data))
  if (length(missing)) {
    cli_abort(c(
      "The sampling `data` is missing column(s) required by {what}.",
      "x" = "Missing: {.field {missing}}"
    ))
  }
  has_na <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1L))]
  if (length(has_na)) {
    cli_abort(c(
      "Columns required by {what} can't contain `NA` values in the sampling `data`.",
      "x" = "Column(s) with `NA`: {.field {has_na}}"
    ))
  }
  not_finite <- vars[vapply(vars, function(v) {
    is.numeric(data[[v]]) && any(!is.finite(data[[v]]))
  }, logical(1L))]
  if (length(not_finite)) {
    cli_abort(c(
      "Numeric columns required by {what} must be finite in the sampling `data`.",
      "x" = "Column(s) with Inf/-Inf: {.field {not_finite}}"
    ))
  }
}

.first_rows <- function(i) cli::cli_vec(i, list("vec-trunc" = 5L))

# Sampling fixed-effect design, with its terms, factor levels, and contrasts
# kept for later sampling predictions. The rank is checked on rows with an
# observed indicator: other rows add no information.
.preferential_sampling_design <- function(spec) {
  data <- spec$data
  terms <- stats::delete.response(stats::terms(spec$formula))
  .check_frame_columns(data, all.vars(terms), "the sampling `formula`")
  mf <- stats::model.frame(terms, data, na.action = stats::na.pass)
  Z <- stats::model.matrix(terms, mf)
  if (!ncol(Z)) {
    cli_abort("The sampling `formula` must have at least one term, e.g., `sampled ~ 1`.")
  }
  r <- as.numeric(data[[spec$response]])
  observed <- !is.na(r)
  Z_obs <- Z[observed, , drop = FALSE]
  q <- qr(Z_obs)
  if (q$rank < ncol(Z)) {
    aliased <- colnames(Z)[q$pivot[-seq_len(q$rank)]]
    cli_abort(c(
      "The sampling design is not full rank on rows with an observed sampling indicator.",
      "x" = "Unidentified column(s): {.code {aliased}}",
      "i" = "A time step whose indicators are all `NA` gives its intercept no information. Drop that level from the sampling `data`, or use a pooled baseline (e.g., `sampled ~ 1`)."
    ))
  }
  # A 0/1 column whose observed rows have a constant indicator where it is 1
  # has an infinite maximum likelihood estimate.
  for (k in seq_len(ncol(Z))) {
    z <- Z_obs[, k]
    if (!all(z %in% c(0, 1)) || !any(z == 1)) next
    value <- unique(r[observed][z == 1])
    if (length(value) == 1L) {
      cli_abort(c(
        "Observed sampling indicators are all {value} wherever column {.code {colnames(Z)[k]}} is 1.",
        "i" = "Its coefficient has no finite estimate. A time step with all or no eligible cells sampled can't have a free intercept."
      ))
    }
  }
  list(
    Z = Z, terms = attr(mf, "terms"),
    xlevels = stats::.getXlevels(attr(mf, "terms"), mf),
    contrasts = attr(Z, "contrasts")
  )
}

#' Prepare the preferential-sampling frame and designs
#'
#' Validates the sampling `data` against the fitted main model and builds the
#' nested `preferential` block of the model data plus R-only metadata. The
#' shared catch design is the main model's fitted fixed-effect design
#' evaluated on the sampling rows (as prediction does), using the fitted
#' `terms` (with `predvars`), `xlevels`, and `contrasts`. Indices are
#' zero-based, as elsewhere in the model data; `rtmb_prepare()` converts them.
#' Row order is kept throughout.
#' @noRd
.prepare_preferential <- function(spec, terms, xlevels, contrasts, X_ij,
                                  mesh, time, time_df) {
  data <- spec$data
  n <- nrow(data)
  r <- as.numeric(data[[spec$response]])
  observed <- !is.na(r)
  if (!any(r[observed] == 0) || !any(r[observed] == 1)) {
    cli_abort("The observed sampling indicators must include both 0s and 1s.")
  }

  xy_cols <- mesh$xy_cols
  if (length(xy_cols) != 2L) {
    cli_abort("Preferential sampling requires a mesh built with known `xy_cols` (e.g., from `make_mesh()`).")
  }
  .check_frame_columns(data, xy_cols, "the mesh coordinates")
  if (!all(vapply(xy_cols, function(v) is.numeric(data[[v]]), logical(1L)))) {
    cli_abort("The sampling `data` coordinates must be numeric.")
  }
  if (identical(time, "_sdmTMB_time")) {
    time_values <- rep(0L, n)
  } else {
    .check_frame_columns(data, time, "the main model's `time`")
    time_values <- data[[time]]
  }
  year_i <- time_df$year_i[match(time_values, time_df$time_from_data)]
  if (anyNA(year_i)) {
    rows <- .first_rows(which(is.na(year_i)))
    cli_abort(c(
      "The sampling `data` has time value(s) outside the fitted time steps (including `extra_time`).",
      "x" = "Row(s): {rows}"
    ))
  }
  key <- data.frame(data[, xy_cols, drop = FALSE], time = time_values)
  if (anyDuplicated(key)) {
    rows <- .first_rows(which(duplicated(key)))
    cli_abort(c(
      "The sampling `data` must have one row per cell and time step.",
      "x" = "Duplicated coordinate/time row(s): {rows}"
    ))
  }

  # Project unique locations once; rows index them by station and time.
  loc_key <- paste(data[[xy_cols[1L]]], data[[xy_cols[2L]]], sep = "\r")
  first <- !duplicated(loc_key)
  station_i <- match(loc_key, loc_key[first]) - 1L
  A_station <- fmesher::fm_basis(mesh$mesh,
    loc = as.matrix(data[first, xy_cols, drop = FALSE]))
  outside <- Matrix::rowSums(A_station) == 0
  if (any(outside)) {
    rows <- .first_rows(which(outside[station_i + 1L]))
    cli_abort(c(
      "Some sampling `data` locations are outside the mesh.",
      "x" = "Row(s): {rows}",
      "i" = "Remove them from the sampling frame or extend the mesh."
    ))
  }

  sampling <- .preferential_sampling_design(spec)

  X_pref <- lapply(seq_along(X_ij), function(m) {
    .check_prediction_safe_terms(terms[[m]])
    fixed_terms <- stats::delete.response(terms[[m]])
    .check_frame_columns(data, all.vars(fixed_terms),
      "the main model's fixed effects")
    X <- tryCatch(
      .fixed_effect_design(fixed_terms, data, xlevels[[m]], contrasts[[m]],
        na.action = stats::na.pass),
      error = function(e) {
        cli_abort(c(
          "Failed to evaluate the main model's fixed effects on the sampling `data`.",
          "x" = conditionMessage(e),
          "i" = "Factor levels must be levels present in the fitted `data`."
        ))
      }
    )
    if (!identical(colnames(X), colnames(X_ij[[m]]))) {
      cli_abort("Internal error: the shared catch design doesn't match the fitted design.")
    }
    X
  })

  list(
    data = list(
      n_pref = n,
      R_i = r,
      Z_ij = sampling$Z,
      X_ij = X_pref,
      offset_i = spec$offset,
      A_station = A_station,
      station_i = station_i,
      year_i = as.integer(year_i),
      include_iid = as.integer(spec$include_iid),
      spatial_xi = as.integer(spec$spatial == "on")
    ),
    info = list(
      spec = spec,
      sampling_terms = sampling$terms,
      sampling_xlevels = sampling$xlevels,
      sampling_contrasts = sampling$contrasts,
      n_observed = sum(observed),
      n_unknown = sum(!observed),
      n_sampled = sum(r == 1, na.rm = TRUE)
    )
  )
}

#' Placeholder preferential-sampling TMB data when the feature is off
#'
#' The C++ template still reads the retired experimental preferential block,
#' so ordinary models pass it with zero-length fields. Remove this with the
#' C++ preferential code.
#' @noRd
.default_preferential_tmb <- function(n_b_j) {
  list(
    n_pref = 0L,
    R_i = numeric(0),
    X_pref_ij = matrix(0, nrow = 0, ncol = n_b_j),
    A_pref = Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, 0L)),
    year_i_pref = integer(0),
    b_pref_type = 0L
  )
}

#' Zero-length parameters of the retired C++ preferential block
#'
#' The C++ template still declares these parameters. Remove this with the
#' C++ preferential code.
#' @noRd
.default_preferential_params <- function() {
  list(
    gamma_0 = numeric(0),
    b_pref = numeric(0),
    log_sigma_b_pref = numeric(0),
    ln_tau_xi = numeric(0),
    ln_kappa_xi = numeric(0),
    xi_s = numeric(0)
  )
}
