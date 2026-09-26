# Preferential sampling: a Bernoulli model for which cells of an eligible
# sampling frame were visited in each time step, linked to a standardized
# expected-catch surface from the main model. The surface is the main
# model's fitted design evaluated on the sampling frame, like prediction
# `newdata`, so the frame's covariate values define the standardization.
#
# Only the RTMB backend supports it. The specification and frame/design
# preparation live here; the likelihood is in `R/rtmb-preferential.R`.

#' Preferential-sampling specification
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Specify a model for which cells of an eligible sampling frame were sampled
#' in each time step, linked to the main model's standardized expected catch.
#' Pass the result to the `preferential` argument of [sdmTMB()].
#'
#' The sampling indicators are modeled jointly with the catch data:
#' \deqn{\mathrm{logit}(p) = Z\gamma + \alpha_t + b h + (b_t - b)(h -
#' \bar{h}_t) + \xi,}
#' where \eqn{Z\gamma} is the sampling `formula`, \eqn{h} is the main model's
#' log expected catch evaluated on the sampling `data` (including its spatial
#' and spatiotemporal fields), \eqn{b} is the preference coefficient
#' `b_pref` and \eqn{b_t} its value in time step \eqn{t} (\eqn{b_t = b} by
#' default; see `coefficient`, which also defines \eqn{\bar{h}_t}),
#' \eqn{\alpha_t} is an optional temporal baseline process (see
#' `baseline`), and \eqn{\xi} is an optional sampling-only spatial field.
#' For a delta model, \eqn{h} is the log of the encounter probability times the positive mean,
#' \eqn{\log(\mathrm{logit}^{-1}(\eta_1)) + \eta_2}. For a Poisson-link delta
#' model (`type = "poisson-link"`), it is \eqn{o + \eta_1 + \eta_2}, where
#' \eqn{o} is `offset`: as in the observation model (but unlike
#' [predict.sdmTMB()], which ignores offsets for these models), the offset
#' scales the expected catch.
#'
#' This feature is under development: it requires the RTMB backend
#' (`control = sdmTMBcontrol(backend = "rtmb")`), a main model with a spatial
#' or spatiotemporal field, and either a single log-link family (Poisson,
#' NB2, Gamma, Tweedie, or lognormal) or a [delta_gamma()] or
#' [delta_lognormal()] family with default links (standard or Poisson-link). Smoothers in the main model must be univariate
#' `s()` terms without `by` variables.
#'
#' @section Fitted models:
#' * [print()] adds a sampling-model section. `tidy(fit, model = "sampling")`
#'   gives the sampling coefficients and `b_pref`,
#'   `tidy(fit, "ran_pars", model = "sampling")` the sampling field's SD and
#'   range and the SDs of any temporal processes, and
#'   `tidy(fit, "ran_vals", model = "sampling")` the preference coefficient
#'   and baseline deviation by time step. [predict_sampling()] gives fitted
#'   sampling probabilities for the sampling frame.
#' * [predict.sdmTMB()], [get_index()], and related functions predict catch
#'   as usual. Their prediction rows are not sampling observations, and their
#'   uncertainty comes from the joint model, including the sampling
#'   likelihood. Use `bias_correct = TRUE` (the default) in [get_index()]:
#'   in simulations, the plug-in index was biased low for delta models.
#' * [logLik()] and [AIC()] use the joint likelihood of the catch data and
#'   the sampling indicators. Don't compare them with a catch-only model or a
#'   model with a different sampling frame. [nobs()] still counts catch
#'   observations only.
#' * [simulate.sdmTMB()] and [residuals.sdmTMB()] are conditional on the
#'   fitted fields and describe the catch data only: they don't check the
#'   sampling model or draw new sampling locations. `simulate()` with
#'   `re_form = NA`, [sdmTMB_cv()], and [project()] aren't supported yet.
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
#'   model's fixed effects and smoothers. The main model is evaluated on these
#'   rows as with `newdata` in [predict.sdmTMB()], so set catchability
#'   covariates (e.g., gear) to reference values here. Holding a smoothed
#'   covariate at a reference value holds its smoother at that value.
#' @param re_form_iid `NA` (default) excludes the main model's IID random
#'   effects from the expected-catch surface, so their grouping columns aren't
#'   needed. `NULL` includes the fitted random intercepts at the group levels
#'   in `data`; these must be levels in the fitted data, and random slopes
#'   aren't supported.
#' @param offset Offset for the main model's expected-catch surface on the
#'   link scale (not an offset for the sampling model). A single value or one
#'   value per row of `data`. The default `0` means unit exposure.
#'   Observation offsets are not copied.
#' @param coefficient How the preference coefficient \eqn{b} varies over
#'   the main model's time steps:
#'   * `"constant"` (default): one coefficient, `b_pref`.
#'   * `"iid"`: \eqn{b_t = \bar{b} + u_t} with independent
#'     \eqn{u_t \sim \mathrm{Normal}(0, \sigma_b^2)}. `b_pref` is
#'     \eqn{\bar{b}}.
#'   * `"rw"`: a random walk, \eqn{b_1} = `b_pref` and
#'     \eqn{b_t = b_{t-1} + d_t} with independent
#'     \eqn{d_t \sim \mathrm{Normal}(0, \sigma_b^2)}. Time steps must be
#'     equally spaced; fill gaps with `extra_time` in [sdmTMB()].
#'
#'   The deviations \eqn{b_t - \bar{b}} (with \eqn{\bar{b}} = `b_pref`)
#'   multiply \eqn{h - \bar{h}_t}, where \eqn{\bar{h}_t} is the mean over
#'   the rows of `data` in time step \eqn{t} of \eqn{h} without its spatial
#'   and spatiotemporal fields (i.e., of the fixed effects, smoothers, IID
#'   effects, and `offset`). The deviations therefore change how strongly
#'   sampling concentrates on high expected catch within a time step, not
#'   (to the extent the fields average to about 0 over the frame) the time
#'   step's sampling rate. With free time-step intercepts in `formula`, this
#'   changes only what the intercepts mean, but it keeps the deviations' SD
#'   estimable: an uncentered deviation also shifts the rate that the
#'   intercept already pins down, and its SD then tends to be estimated as 0.
#'   The fields are left out of \eqn{\bar{h}_t} because their frame mean
#'   would link every sampling row to every mesh vertex of its time step and
#'   make fitting slow. So if the main model's time-step means come from a
#'   random-walk or AR(1) spatiotemporal field rather than fixed effects,
#'   consider adding time-step fixed effects (e.g., `0 + factor(year)`) to the
#'   main model.
#'
#'   Every time step of the main model (including `extra_time`) gets a
#'   coefficient; time steps without observed sampling indicators get theirs
#'   from the IID or random-walk distribution alone. At least two time steps
#'   are required, but a handful of time steps carries little information
#'   about \eqn{\sigma_b}: check its estimate and interval, and prefer
#'   `"constant"` unless the data support more.
#' @param baseline An optional temporal process for the sampling intercept,
#'   added to the sampling `formula`: `"off"` (default), `"iid"`, or `"rw"`,
#'   defined as for `coefficient`. The `formula`'s intercept is the mean
#'   (`"iid"`) or first-time-step value (`"rw"`) of the baseline. This
#'   shrinks the time-step baselines towards each other instead of
#'   estimating them freely with `0 + factor(year)`, and so requires a
#'   `formula` with an intercept and without free time-step effects.
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
                                  coefficient = c("constant", "iid", "rw"),
                                  baseline = c("off", "iid", "rw"),
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
  coefficient <- match.arg(coefficient)
  baseline <- match.arg(baseline)
  spatial <- match.arg(spatial)
  structure(
    list(
      formula = formula, data = data, response = response,
      include_iid = is.null(re_form_iid),
      offset = rep_len(as.numeric(offset), nrow(data)),
      coefficient = coefficient, baseline = baseline, spatial = spatial
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
                                         family, areal, mesh, mesh_missing,
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
  bars <- unlist(lapply(formulas, reformulas::findbars))
  random_slopes <- !all(vapply(bars, function(b) identical(b[[2L]], 1),
    logical(1L)))
  smooths <- term_labels[get_smooth_terms(term_labels)]
  # Ordinary univariate smooths only: no `by` variables, tensor products,
  # or Markov random fields.
  univariate <- vapply(smooths, function(x) {
    s <- eval(str2expression(x))
    length(s$term) == 1L && identical(s$by, "NA") &&
      !inherits(s, c("t2.smooth.spec", "tensor.smooth.spec", "mrf.smooth.spec"))
  }, logical(1L))
  log_link_family <- !delta && !multi_family &&
    family$family[[1L]] %in% .preferential_families &&
    identical(family$link[[1L]], "log")
  # Standard deltas use logit/log links; Poisson-link deltas log/log.
  supported_delta <- delta && !multi_family &&
    identical(family$family[[1L]], "binomial") &&
    family$family[[2L]] %in% c("Gamma", "lognormal") && (
      (identical(family$type, "standard") &&
        identical(unname(family$link), c("logit", "log"))) ||
      (identical(family$type, "poisson_link_delta") &&
        identical(unname(family$link), c("log", "log"))))
  unsupported <- c(
    "delta models other than `delta_gamma()` or `delta_lognormal()` with default links" =
      delta && !multi_family && !supported_delta,
    "families other than log-link Poisson, NB2, Gamma, Tweedie, or lognormal" =
      !delta && !multi_family && !log_link_family,
    "multi-family models" = multi_family,
    "areal (SAR/CAR) models" = areal,
    "barrier meshes" = "spde_barrier" %in% names(mesh),
    "anisotropy" = isTRUE(anisotropy),
    "`time_varying`" = !is.null(time_varying),
    "`spatial_varying`" = !is.null(spatial_varying),
    "`nonlocal_formula`" = !is.null(nonlocal_formula),
    "threshold (`breakpt()`/`logistic()`) terms" =
      any(grepl("^(breakpt|logistic)\\(", term_labels)),
    "smoothers other than univariate `s()` terms without `by`" =
      !all(univariate),
    "including random slopes (`re_form_iid = NULL`)" =
      random_slopes && spec$include_iid,
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

# Families whose log-link predictor is the log of the response mean.
.preferential_families <- c("poisson", "nbinom2", "Gamma", "tweedie",
  "lognormal")

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

# Check the time grid and sampling design for the temporal preference and
# baseline processes, and return the number of random deviations of each: one
# per time step for IID, one per step after the first for a random walk.
.check_preferential_temporal <- function(spec, Z_obs, year_obs, time_df) {
  n_t <- nrow(time_df)
  type <- c(coefficient = spec$coefficient, baseline = spec$baseline)
  type[type == "off"] <- "constant"
  temporal <- type != "constant"
  if (any(temporal) && n_t < 2L) {
    cli_abort(c(
      "A temporal preference coefficient or baseline needs at least two time steps.",
      "i" = "Use `coefficient = \"constant\"` and `baseline = \"off\"`."
    ))
  }
  if (any(type == "rw")) {
    t <- time_df$time_from_data
    if (!is.numeric(t) || length(unique(diff(t))) > 1L) {
      missed <- if (is.numeric(t)) find_missing_time(t)
      cli_abort(c(
        "A random-walk preference coefficient or baseline needs numeric, equally spaced time steps.",
        "i" = "Fill gaps in time with `extra_time` in `sdmTMB()`.",
        if (length(missed)) {
          "i" = paste0("`extra_time = c(", paste(missed, collapse = ", "), ")`")
        }
      ))
    }
  }
  if (temporal[["baseline"]]) {
    rank <- function(x) qr(x)$rank
    k <- rank(Z_obs)
    if (rank(cbind(Z_obs, 1)) > k) {
      cli_abort(c(
        "`baseline` needs a sampling `formula` with an intercept.",
        "i" = "The intercept is the mean (IID) or first value (random walk) of the baseline."
      ))
    }
    steps <- stats::model.matrix(~ 0 + factor(year_obs))
    if (ncol(steps) > 1L && rank(cbind(Z_obs, steps)) == k) {
      cli_abort(c(
        "The sampling `formula` already has a free baseline for each time step.",
        "i" = "Use either `baseline` or time-step effects such as `0 + factor(year)`, not both."
      ))
    }
  }
  ifelse(type == "iid", n_t, ifelse(type == "rw", n_t - 1L, 0L))
}

#' Prepare the preferential-sampling frame and designs
#'
#' Validates the sampling `data` against the fitted main model and builds the
#' nested `preferential` block of the model data, its parameters, and R-only
#' metadata. The shared catch design is the main model evaluated on the
#' sampling rows by `.newdata_design()`, as prediction evaluates `newdata`.
#' `model` holds the fitted pieces that it needs (a fitted object, or the same
#' pieces while fitting); `X_ij` is the fitted fixed-effect design. Indices are
#' zero-based, as elsewhere in the model data; `rtmb_prepare()` converts them.
#' Row order is kept throughout.
#' @noRd
.prepare_preferential <- function(spec, model, X_ij, mesh, time, time_df) {
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
  locations <- .project_unique_locations(mesh$mesh, data, xy_cols)
  A_station <- locations$A
  station_i <- locations$index
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

  # Every covariate of the main model must be complete on every row: the
  # shared surface is evaluated on all of them.
  lapply(model$terms, .check_prediction_safe_terms)
  bars <- reformulas::findbars(model$smoothers$formula_no_sm)
  include_iid <- spec$include_iid && length(bars) > 0L
  .check_frame_columns(data, unique(c(
    unlist(lapply(model$terms, function(x) all.vars(stats::delete.response(x)))),
    all.vars(stats::delete.response(
      stats::terms(model$smoothers$formula_no_bars))),
    if (include_iid) barnames(bars)
  )), "the main model")
  design <- tryCatch(
    .newdata_design(model, data, include_iid = include_iid,
      new_levels = "error", na.action = stats::na.pass),
    error = function(e) {
      cli_abort(c(
        "Failed to evaluate the main model on the sampling `data`.",
        "x" = conditionMessage(e),
        "i" = "Factor levels must be levels present in the fitted `data`."
      ))
    }
  )
  for (m in seq_along(X_ij)) {
    if (!identical(colnames(design$X_ij[[m]]), colnames(X_ij[[m]]))) {
      cli_abort("Internal error: the shared catch design doesn't match the fitted design.")
    }
  }

  # Start the sampling coefficients at the sampling-only logistic regression
  # estimates; separation was ruled out above.
  start <- tryCatch(
    stats::glm.fit(sampling$Z[observed, , drop = FALSE], r[observed],
      family = stats::binomial())$coefficients,
    error = function(e) rep(0, ncol(sampling$Z)),
    warning = function(w) rep(0, ncol(sampling$Z))
  )
  xi <- spec$spatial == "on"
  n_dev <- .check_preferential_temporal(spec, sampling$Z[observed, , drop = FALSE],
    year_i[observed], time_df)

  list(
    data = list(
      n_pref = n,
      R_i = r,
      Z_ij = sampling$Z,
      X_ij = design$X_ij,
      Zs = design$Zs,
      Xs = design$Xs,
      Zt_list = design$Zt_list,
      offset_i = spec$offset,
      A_station = A_station,
      station_i = station_i,
      year_i = as.integer(year_i),
      include_iid = as.integer(spec$include_iid),
      spatial_xi = as.integer(xi),
      # 0 = none, 1 = IID, 2 = random walk
      coefficient_type = match(spec$coefficient, c("constant", "iid", "rw")) - 1L,
      baseline_type = match(spec$baseline, c("off", "iid", "rw")) - 1L
    ),
    parameters = c(
      list(gamma_pref = unname(start), b_pref = 0),
      if (xi) {
        list(ln_tau_xi = 0, ln_kappa_xi = 0, xi_s = rep(0, ncol(A_station)))
      },
      if (spec$coefficient != "constant") {
        list(ln_sigma_b_pref = 0, b_pref_dev = rep(0, n_dev[["coefficient"]]))
      },
      if (spec$baseline != "off") {
        list(ln_sigma_alpha_pref = 0,
          alpha_pref_dev = rep(0, n_dev[["baseline"]]))
      }
    ),
    random = c(if (xi) "xi_s",
      if (spec$coefficient != "constant") "b_pref_dev",
      if (spec$baseline != "off") "alpha_pref_dev"),
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

# Post-fit methods ------------------------------------------------------------

.check_preferential_fit <- function(object, what) {
  if (is.null(object$preferential)) {
    cli_abort("{what} requires a model fit with `preferential` (see `preferential_sampling()`).")
  }
}

# Sampling-model rows for tidy(x, model = "sampling"). The coefficients and
# `b_pref` are fixed effects. The sampling field's SD and range and the SDs of
# the temporal processes are random-effect parameters, with intervals on the
# log scale. "ran_vals" gives the preference coefficient and baseline
# deviation by time step.
.tidy_sampling <- function(x, effects, conf.int, crit, trans) {
  .check_preferential_fit(x, "`model = \"sampling\"`")
  if (!effects %in% c("fixed", "ran_pars", "ran_vals")) {
    cli_abort("With `model = \"sampling\"`, `effects` must be \"fixed\", \"ran_pars\", or \"ran_vals\".")
  }
  spec <- x$preferential$spec
  if (effects == "fixed") {
    est <- as.list(x$sd_report, "Estimate")
    se <- as.list(x$sd_report, "Std. Error")
    out <- data.frame(
      term = c(colnames(x$tmb_data$preferential$Z_ij), "b_pref"),
      estimate = c(est$gamma_pref, est$b_pref),
      std.error = c(se$gamma_pref, se$b_pref),
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      out$conf.low <- as.numeric(trans(out$estimate - crit * out$std.error))
      out$conf.high <- as.numeric(trans(out$estimate + crit * out$std.error))
    }
    out$estimate <- as.numeric(trans(out$estimate))
    if (!identical(trans, I)) out$std.error <- NULL
  } else {
    est <- as.list(x$sd_report, "Estimate", report = TRUE)
    se <- as.list(x$sd_report, "Std. Error", report = TRUE)
    if (effects == "ran_pars") {
      terms <- c(
        if (spec$spatial == "on") c("range_xi", "sigma_xi"),
        if (spec$coefficient != "constant") "sigma_b_pref",
        if (spec$baseline != "off") "sigma_alpha_pref"
      )
      log_est <- unlist(est[paste0("log_", terms)])
      log_se <- unlist(se[paste0("log_", terms)])
      lower <- exp(log_est - crit * log_se)
      upper <- exp(log_est + crit * log_se)
    } else {
      terms <- c(
        if (spec$coefficient != "constant") "b_pref_t",
        if (spec$baseline != "off") "alpha_pref_t"
      )
      lower <- unlist(est[terms]) - crit * unlist(se[terms])
      upper <- unlist(est[terms]) + crit * unlist(se[terms])
    }
    out <- data.frame(
      term = rep(terms, lengths(est[terms])),
      estimate = as.numeric(unlist(est[terms])),
      std.error = as.numeric(unlist(se[terms])),
      conf.low = as.numeric(lower), conf.high = as.numeric(upper),
      stringsAsFactors = FALSE
    )
    if (effects == "ran_vals") {
      time <- x$time_lu$time_from_data
      out <- data.frame(out[1L], time = rep(time, length.out = nrow(out)),
        out[-1L])
    }
    if (!conf.int) out$conf.low <- out$conf.high <- NULL
  }
  row.names(out) <- NULL
  tibble::as_tibble(out)
}

# Sampling-model section of print.sdmTMB().
print_sampling <- function(x) {
  info <- x$preferential
  spec <- info$spec
  offset <- unique(range(spec$offset))
  cat("\nSampling model (preferential sampling): ----------------------\n")
  cat("Formula: ", deparse1(spec$formula), "\n", sep = "")
  cat("Sampling frame: ", nrow(spec$data), " rows; ", info$n_observed,
    " observed (", info$n_sampled, " sampled), ", info$n_unknown,
    " unknown\n", sep = "")
  cat("Shared target: log standardized expected catch (IID effects ",
    if (spec$include_iid) "included" else "excluded", "; offset ",
    paste(format(offset, digits = 3L), collapse = " to "), ")\n", sep = "")
  process <- c(iid = "IID by time step", rw = "random walk over time steps")
  cat("Preference coefficient: ", if (spec$coefficient == "constant") {
    "constant"
  } else {
    process[[spec$coefficient]]
  }, "\n", sep = "")
  if (spec$baseline != "off") {
    cat("Baseline: ", process[[spec$baseline]], "\n", sep = "")
  }
  cat("\n")
  b <- tidy(x, model = "sampling", silent = TRUE)
  mm <- cbind(round(b$estimate, 2L), round(b$std.error, 2L))
  dimnames(mm) <- list(b$term, c("coef.est", "coef.se"))
  print(mm)
  cat("\n")
  r <- tidy(x, "ran_pars", model = "sampling", silent = TRUE)
  labels <- c(range_xi = "Sampling field range",
    sigma_xi = "Sampling field SD",
    sigma_b_pref = "Preference coefficient SD over time",
    sigma_alpha_pref = "Baseline SD over time")
  for (i in seq_len(nrow(r))) {
    cat(labels[[r$term[i]]], ": ", mround(r$estimate[i], 2L), "\n", sep = "")
  }
  cat("The criterion below includes the sampling likelihood.\n")
}

#' Predict sampling probabilities from a preferential-sampling model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fitted sampling probabilities, and their predictor pieces, for the
#' sampling frame of a model fit with `preferential` (see
#' [preferential_sampling()]). Predictions for a new sampling frame or new
#' time steps are not supported yet.
#'
#' @param object A model fit by [sdmTMB()] with `preferential`.
#' @param type `"response"` for sampling probabilities or `"link"` for the
#'   logit scale. Applies to `est` and to the draws with `nsim > 0`; the
#'   predictor pieces are always on the logit scale.
#' @param nsim If `> 0`, return a matrix of `nsim` draws (columns) of `est`
#'   for each frame row (rows), from the joint precision matrix of all fixed
#'   and random effects. Use these draws for uncertainty intervals.
#'
#' @return
#' With `nsim = 0`, the sampling `data` in its original row order with these
#' columns added:
#' * `est`: the sampling probability (or its logit, for `type = "link"`).
#' * `est_target`: the shared target \eqn{h}, the main model's log
#'   standardized expected catch.
#' * `est_fixed`: the sampling formula's contribution \eqn{Z\gamma}.
#' * `est_baseline`: the baseline deviation \eqn{\alpha_t}, if `baseline`
#'   was used.
#' * `est_preference`: the preference contribution, \eqn{b h} plus
#'   \eqn{(b_t - b)(h - \bar{h}_t)} for a temporal `coefficient`.
#' * `est_xi`: the sampling field \eqn{\xi}, if estimated.
#'
#' On the link scale, `est` is the sum of `est_fixed`, `est_baseline`,
#' `est_preference`, and `est_xi`. With `nsim > 0`, a matrix with one row per frame row.
#' @export
predict_sampling <- function(object, type = c("response", "link"), nsim = 0) {
  assert_that(inherits(object, "sdmTMB"))
  .check_preferential_fit(object, "`predict_sampling()`")
  type <- match.arg(type)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 0 ||
      nsim != floor(nsim)) {
    cli_abort("`nsim` must be one non-negative whole number.")
  }
  reinitialize(object)
  inv <- if (type == "response") stats::plogis else identity
  lp <- object$tmb_obj$env$last.par.best
  if (nsim > 0) {
    draws <- .joint_par_draws(object, lp, nsim)
    out <- apply(draws, 2L, function(par) {
      object$tmb_obj$report(par)$sampling_eta_i
    })
    return(inv(matrix(out, ncol = nsim)))
  }
  r <- object$tmb_obj$report(lp)
  nd <- object$preferential$spec$data
  nd$est <- inv(r$sampling_eta_i)
  nd$est_target <- r$sampling_target_i
  nd$est_fixed <- r$sampling_fixed_i
  if (object$preferential$spec$baseline != "off") {
    nd$est_baseline <- r$sampling_baseline_i
  }
  nd$est_preference <- r$sampling_preference_i
  if (object$preferential$spec$spatial == "on") nd$est_xi <- r$sampling_field_i
  nd
}
