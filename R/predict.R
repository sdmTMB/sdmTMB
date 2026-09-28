# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Make predictions from an \pkg{sdmTMB} model. Predictions can be made on the
#' original data or on new data.
#'
#' @param object A model fitted with [sdmTMB()].
#' @param newdata A data frame to make predictions on. It should contain the
#'   same predictor columns as the fitted data and, for spatiotemporal models,
#'   a time column with the same name as in the fitted data.
#' @param type Should predictions be returned in link space (default) or
#'   response space? Standard errors (`se_fit = TRUE`) are only available in
#'   link space.
#' @param se_fit Should standard errors on predictions be calculated? Warning:
#'   can be slow for large datasets or high-resolution projections when random
#'   fields are included. For faster uncertainty estimation, either use
#'   `re_form = NA` to exclude random fields or use the `nsim` argument to
#'   simulate from the joint precision matrix. Requires `type = "link"`; for
#'   response-scale uncertainty, use `nsim`.
#' @param return_tmb_object `r lifecycle::badge("deprecated")` Logical. If
#'   `TRUE`, include the TMB object in a list-format output. Instead, pass the
#'   fitted model and `newdata` directly to [get_index()] or [get_cog()].
#' @param re_form `NULL` to include all spatial/spatiotemporal random fields in
#'   predictions. `~0` or `NA` for population-level predictions (predictions
#'   excluding spatial/spatiotemporal random fields). Often used with
#'   `se_fit = TRUE` to visualize marginal effects. Does not affect
#'   [get_index()] calculations.
#' @param re_form_iid `NULL` to include all IID random intercepts/slopes in the
#'   predictions. `~0` or `NA` for population-level predictions. No other
#'   options (e.g., some but not all random intercepts) are not yet implemented.
#'   Only affects predictions with `newdata`. This *does* affect [get_index()].
#' @param allow_new_levels Logical or `NULL`. Similar to \pkg{glmmTMB}'s
#'   `allow.new.levels`.
#'   Allows predictions for previously unobserved levels in random effect
#'   grouping variables. If `NULL` (default), new levels are allowed when
#'   `re_form_iid = NA` or `re_form_iid = ~0` and a warning is issued
#'   otherwise. If `TRUE`, new levels are explicitly allowed. If `FALSE`, a
#'   warning is issued if new levels are found. New levels are always treated
#'   as population-level predictions for the IID random effects
#'   (i.e., random effect value = 0).
#' @param nsim If `> 0`, simulate from the joint precision matrix with `nsim`
#'   draws. Returns a matrix with one row per prediction location and one column
#'   per draw. By default, each column represents one draw of the linear predictor
#'   in link space; use `type = "response"` for response-space draws. Simulating
#'   from the joint precision matrix accounts for uncertainty in both fixed and
#'   random effects. Use this to derive uncertainty on predictions (e.g.,
#'   `apply(x, 1, sd)`) or propagate uncertainty to derived quantities. This is
#'   the fastest way to characterize spatial uncertainty with sdmTMB.
#' @param sample_fe Logical. When `nsim > 0`, sample uncertainty in the fixed
#'   effects and other estimated parameters? If `FALSE`, these are held at
#'   their estimated values and only the random effects (random fields, IID
#'   random effects, time-varying coefficients, and smoother coefficients) are
#'   drawn from their distribution conditional on the estimated parameters
#'   (similar to `obj$MC()` in \pkg{TMB}). Fixed effects are held at their
#'   estimates even with REML. See also the same argument in [project()].
#' @param sims_var Experimental: Which TMB reported variable from the model
#'   should be extracted from the joint precision matrix simulation draws?
#'   Defaults to link-space predictions. Other options are `"omega_s"`,
#'   `"zeta_s"`, `"epsilon_st"`, `"est_rf"`, and `"est_non_rf"` (as described
#'   below); the model must include the term. Options other than `"est"` are
#'   returned in link space and cannot be combined with `type = "response"`.
#'   For other reported variables, use `return_tmb_report = TRUE`.
#' @param mcmc_samples See `extract_mcmc()` in the
#'   \href{https://github.com/sdmTMB/sdmTMBextra}{sdmTMBextra} package for
#'   more details and the
#'   \href{https://sdmTMB.github.io/sdmTMB/articles/bayesian.html}{Bayesian vignette}.
#'   If specified, the predict function will return a matrix of a similar form
#'   as if `nsim > 0` but representing Bayesian posterior samples from the Stan
#'   model.
#' @param nonlocal_newdata An optional data frame overriding the
#'   `nonlocal_formula` covariate field used for prediction (e.g., for a
#'   counterfactual/scenario surface), with the same requirements as
#'   `nonlocal_data` in [sdmTMB()]. `newdata`'s x/y and time columns
#'   always determine *where* predictions are projected to; this argument only
#'   controls where the underlying diffused covariate values come from.
#'   Defaults to `NULL`: if a grid was supplied at fit time, the fitted field
#'   is reused as-is (so `newdata` need not contain the diffusion covariate
#'   columns); otherwise the field is rebuilt from `newdata`'s own covariate
#'   columns, as before. Also applies when `newdata = NULL`, in which case
#'   predictions are projected onto the fitted data.
#' @param model Which component to predict from delta/hurdle models when `nsim >
#'   0` or `mcmc_samples` is supplied. `NA` (default) returns the combined
#'   prediction from both components; `1` returns the binomial component only; `2`
#'   returns the positive component only. Predictions are on the link or response
#'   scale depending on `type`. For regular predictions (without simulation),
#'   both components are returned. See the [delta-model
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html).
#' @param offset A numeric vector of optional offset values, one per row of
#'   `newdata`. If `NULL` (default), predictions on the fitted data (`newdata =
#'   NULL`) use the offset from the fitted model and predictions with `newdata`
#'   use an offset of 0.
#' @param return_tmb_report Logical: return the output from the TMB
#'   report? For regular prediction, this is all the reported variables
#'   at the MLE parameter values. For `nsim > 0` or when `mcmc_samples`
#'   is supplied, this is a list with one element per sample; each element
#'   contains the report output for that sample. With `newdata = NULL`, this
#'   is the fitted model's report unless another argument requires projecting
#'   onto the fitted data, in which case it is the projection report.
#' @param return_tmb_data Logical: return formatted data for TMB? Used
#'   internally. With `newdata = NULL`, the fitted data are prepared as
#'   prediction data (with the fitted offset unless `offset` is supplied).
#' @param ... Unused.
#'
#' @return
#' If `return_tmb_object = FALSE` (and `nsim = 0` and `mcmc_samples = NULL`):
#'
#' A data frame:
#' * `est`: Estimate in link or response space, depending on `type`
#' * `est_non_rf`: Estimate from everything except spatial/spatiotemporal random fields (fixed effects, random intercepts, time-varying effects, etc.)
#' * `est_rf`: Estimate from all random fields combined
#' * `omega_s`: Spatial random field (models consistent spatial patterns)
#' * `zeta_s`: Spatially varying coefficient field (models how effects vary across space)
#' * `epsilon_st`: Spatiotemporal random field (models spatial patterns that vary over time)
#' * `nl_*`: Nonlocal transformed covariate values (one column per
#'   nonlocal term; available when `nonlocal_formula` terms were fitted)
#'
#' Delta/hurdle models return component-specific columns with `1` and `2`
#' suffixes for the binomial and positive components, respectively (e.g.,
#' `est1`, `est2`, `omega_s1`, `omega_s2`). With `type = "response"`,
#' `est` is the combined response-scale prediction.
#'
#' If `return_tmb_object = TRUE` (and `nsim = 0` and `mcmc_samples = NULL`):
#'
#' A list:
#' * `data`: The data frame described above
#' * `report`: The TMB report on parameter values
#' * `obj`: The TMB object returned from the prediction run
#' * `fit_obj`: The original TMB model object
#'
#' In this case, you likely only need the `data` element as an end user.
#' The other elements are included for other functions.
#'
#' If `nsim > 0` or `mcmc_samples` is not `NULL`:
#'
#' A matrix:
#'
#' * Columns represent samples
#' * Rows represent predictions, with one row per row of `newdata`
#'
#' @export
#'
#' @examplesIf ggplot2_installed()
#'
#' d <- pcod_2011
#' mesh <- make_mesh(d, c("X", "Y"), cutoff = 30) # a coarse mesh for example speed
#' m <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # Predictions at original data locations -------------------------------
#'
#' predictions <- predict(m)
#' head(predictions)
#'
#' predictions$resids <- residuals(m) # randomized quantile residuals
#'
#' library(ggplot2)
#' ggplot(predictions, aes(X, Y, col = resids)) + scale_colour_gradient2() +
#'   geom_point() + facet_wrap(~year)
#' hist(predictions$resids)
#' qqnorm(predictions$resids); abline(a = 0, b = 1)
#'
#' # Predictions on new data ----------------------------------------------
#'
#' qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
#' predictions <- predict(m, newdata = qcs_grid_2011)
#'
#' \donttest{
#' # A short function for plotting predictions:
#' plot_map <- function(dat, column = est) {
#'   ggplot(dat, aes(X, Y, fill = {{ column }})) +
#'     geom_raster() +
#'     facet_wrap(~year) +
#'     coord_fixed()
#' }
#'
#' plot_map(predictions, exp(est)) +
#'   scale_fill_viridis_c(trans = "sqrt") +
#'   ggtitle("Prediction (fixed effects + all random effects)")
#'
#' plot_map(predictions, exp(est_non_rf)) +
#'   ggtitle("Prediction (fixed effects and any time-varying effects)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(predictions, est_rf) +
#'   ggtitle("All random field estimates") +
#'   scale_fill_gradient2()
#'
#' plot_map(predictions, omega_s) +
#'   ggtitle("Spatial random effects only") +
#'   scale_fill_gradient2()
#'
#' plot_map(predictions, epsilon_st) +
#'   ggtitle("Spatiotemporal random effects only") +
#'   scale_fill_gradient2()
#'
#' # Visualizing a marginal effect ----------------------------------------
#'
#' # See the visreg package or the ggeffects::ggeffect() or
#' # ggeffects::ggpredict() functions
#' # To do this manually:
#'
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$depth_scaled2 <- nd$depth_scaled^2
#'
#' # Because this is a spatiotemporal model, you'll need at least one time
#' # value. For these population-level predictions, if time isn't also a fixed
#' # effect, it doesn't matter what you pick:
#' nd$year <- 2011L # L: integer to match original data
#' p <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Plotting marginal effect of a spline ---------------------------------
#'
#' m_gam <- sdmTMB(
#'  data = d, formula = density ~ 0 + as.factor(year) + s(depth_scaled, k = 5),
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#' if (require("visreg", quietly = TRUE)) {
#'   visreg::visreg(m_gam, "depth_scaled")
#' }
#'
#' # or manually:
#' nd <- data.frame(depth_scaled =
#'   seq(min(d$depth_scaled), max(d$depth_scaled), length.out = 100))
#' nd$year <- 2011L
#' p <- predict(m_gam, newdata = nd, se_fit = TRUE, re_form = NA)
#' ggplot(p, aes(depth_scaled, exp(est),
#'   ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
#'   geom_line() + geom_ribbon(alpha = 0.4)
#'
#' # Forecasting ----------------------------------------------------------
#' mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#'
#' unique(d$year)
#' m <- sdmTMB(
#'   data = d, formula = density ~ 1,
#'   spatiotemporal = "AR1", # using AR(1) to have something to forecast with
#'   extra_time = 2019L, # `L` for integer to match our data
#'   spatial = "off",
#'   time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # Add a year to our grid:
#' grid2019 <- qcs_grid_2011[qcs_grid_2011$year == max(qcs_grid_2011$year), ]
#' grid2019$year <- 2019L # `L` because `year` is an integer in the data
#' qcsgrid_forecast <- rbind(qcs_grid_2011, grid2019)
#'
#' predictions <- predict(m, newdata = qcsgrid_forecast)
#' plot_map(predictions, exp(est)) +
#'   scale_fill_viridis_c(trans = "log10")
#' plot_map(predictions, epsilon_st) +
#'   scale_fill_gradient2()
#'
#' # Estimating local trends ----------------------------------------------
#'
#' d <- pcod
#' d$year_scaled <- as.numeric(scale(d$year))
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#' m <- sdmTMB(data = d, formula = density ~ depth_scaled + depth_scaled2,
#'   mesh = mesh, family = tweedie(link = "log"),
#'   spatial_varying = ~ 0 + year_scaled, time = "year", spatiotemporal = "off")
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#' nd$year_scaled <- (nd$year - mean(d$year)) / sd(d$year)
#' p <- predict(m, newdata = nd)
#'
#' plot_map(subset(p, year == 2003), zeta_s_year_scaled) + # pick any year
#'   ggtitle("Spatial slopes") +
#'   scale_fill_gradient2()
#'
#' plot_map(p, est_rf) +
#'   ggtitle("Random field estimates") +
#'   scale_fill_gradient2()
#'
#' plot_map(p, exp(est_non_rf)) +
#'   ggtitle("Prediction (fixed effects only)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(p, exp(est)) +
#'   ggtitle("Prediction (fixed effects + all random effects)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#' }

predict.sdmTMB <- function(object, newdata = NULL,
  type = c("link", "response"),
  se_fit = FALSE,
  re_form = NULL,
  re_form_iid = NULL,
  allow_new_levels = NULL,
  nsim = 0,
  sims_var = "est",
  sample_fe = TRUE,
  model = c(NA, 1, 2),
  offset = NULL,
  mcmc_samples = NULL,
  nonlocal_newdata = NULL,
  return_tmb_object = deprecated(),
  return_tmb_report = FALSE,
  return_tmb_data = FALSE,
  ...) {

  # lifecycle only warns for direct user calls, so check here:
  if (is_present(return_tmb_object)) {
    lifecycle::deprecate_soft(
      "1.2.0",
      "predict.sdmTMB(return_tmb_object = )",
      details = paste(
        "Pass the fitted model and `newdata` directly to `get_index()`,",
        "`get_cog()`, `get_eao()`, or `get_weighted_average()`."
      )
    )
  } else {
    return_tmb_object <- FALSE
  }
  if (is_visreg_call()) {
    return(predict_visreg(object, newdata, se_fit = isTRUE(list(...)$se.fit)))
  }
  if (missing(model) && .has_delta_attr(object)) {
    model <- attr(object, "delta_model_predict") # for ggpredict
  }
  predict_sdmTMB(object,
    newdata = newdata, type = type, se_fit = se_fit, re_form = re_form,
    re_form_iid = re_form_iid, allow_new_levels = allow_new_levels,
    nsim = nsim, sims_var = sims_var, sample_fe = sample_fe, model = model,
    offset = offset, mcmc_samples = mcmc_samples,
    nonlocal_newdata = nonlocal_newdata,
    return_tmb_object = return_tmb_object,
    return_tmb_report = return_tmb_report, return_tmb_data = return_tmb_data
  )
}

# The body of predict.sdmTMB(), without the call-stack inspection for visreg.
predict_sdmTMB <- function(object, newdata = NULL, type = "link",
                           se_fit = FALSE, re_form = NULL, re_form_iid = NULL,
                           allow_new_levels = NULL, nsim = 0, sims_var = "est",
                           sample_fe = TRUE, model = NA, offset = NULL,
                           mcmc_samples = NULL,
                           nonlocal_newdata = NULL,
                           return_tmb_object = FALSE, return_tmb_report = FALSE,
                           return_tmb_data = FALSE) {
  req <- predict_request(
    object = object, newdata = newdata, type = type, se_fit = se_fit,
    re_form = re_form, re_form_iid = re_form_iid,
    allow_new_levels = allow_new_levels, nsim = nsim, model = model,
    offset = offset, mcmc_samples = mcmc_samples, sims_var = sims_var,
    sample_fe = sample_fe, return_tmb_data = return_tmb_data,
    nonlocal_newdata = nonlocal_newdata
  )

  tmb_data <- object$tmb_data
  if (is.null(tmb_data$link_pred) && !is.null(tmb_data$link)) {
    tmb_data$link_pred <- tmb_data$link
  }
  tmb_data$do_predict <- 1L
  has_nonlocal <- !is.null(object$nonlocal_parsed)

  if (req$project) {
    prep <- predict_prepare(object, req, tmb_data, nonlocal_newdata)
    tmb_data <- prep$tmb_data
    if (return_tmb_data) {
      return(tmb_data)
    }
    if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE

    objective <- predict_objective(object, tmb_data)
    new_tmb_obj <- objective$obj
    lp <- objective$lp

    if (req$nsim > 0 || !is.null(req$mcmc_samples)) {
      r <- predict_draw_reports(object, new_tmb_obj, lp, req)
      if (return_tmb_report) return(r)
      return(predict_draws(r, req, object, tmb_data, prep$nd))
    }

    r <- new_tmb_obj$report(lp)
    if (return_tmb_report) return(r)
    obj <- new_tmb_obj

    if (req$se_fit) {
      sr <- sdreport_sdmTMB(new_tmb_obj, bias.correct = FALSE)
      sr_est <- as.list(sr, "Estimate", report = TRUE)
      sr_se <- as.list(sr, "Std. Error", report = TRUE)
    }

    cols <- list(est = predict_est(if (req$se_fit) sr_est else r, req, tmb_data,
      req$type, req$model))
    components <- if (req$has_two_components) {
      predict_components(r, req, tmb_data, req$type)[c("est1", "est2")]
    }
    se <- if (req$se_fit) {
      list(est_se = predict_select(sr_se, req, tmb_data, req$model))
    }
    terms <- if (!req$pop_pred) {
      predict_terms(
        lapply(predict_term_reports, function(x) r[[x]]),
        object, req$family_spec, tmb_data$proj_family_id + 1L
      )
    }
    cols <- c(cols, components, terms, se)
    if (has_nonlocal) {
      cols <- c(predict_nonlocal_cols(object, r$proj_covariate_diffusion_values), cols)
    }
  } else { # Link prediction on the fitted data from the fitted report:
    reinitialize(object)
    lp <- object$tmb_obj$env$last.par.best
    r <- object$tmb_obj$report(lp)
    if (return_tmb_report) return(r)

    cols <- c(list(est = r$eta_i[, 1]), predict_terms(
      predict_fitted_term_reports(r, object),
      object, req$family_spec, req$family_spec$family_id_i
    ))
    if (has_nonlocal) {
      cols <- c(predict_nonlocal_cols(object, r$covariate_diffusion_values), cols)
    }
    obj <- object
  }

  predict_return(cols, req, object, r, obj, tmb_data, return_tmb_object)
}

# Build the prediction objective from `tmb_data` and return it with the fitted
# parameter vector `lp`. Fits saved with `parlist` and `last.par.best` need no
# live fitted objective; older fits recover parameters from `object$tmb_obj`.
predict_objective <- function(object, tmb_data) {
  has_saved_fit <- !is.null(object$parlist) && !is.null(object$last.par.best)
  if (!has_saved_fit) reinitialize(object)
  obj <- make_sdmTMB_adfun(
    data = tmb_data,
    profile = object$control$profile,
    parameters = if (has_saved_fit) object$parlist else get_pars(object),
    map = object$tmb_map,
    random = object$tmb_random,
    backend = backend_sdmTMB(object),
    silent = TRUE
  )
  if (has_saved_fit) {
    lp <- object$last.par.best
  } else {
    obj$fn(object$model$par)
    lp <- obj$env$last.par.best
  }
  list(obj = obj, lp = lp)
}

# Nonlocal covariate values as named prediction columns.
predict_nonlocal_cols <- function(object, values) {
  colnames(values) <- .nonlocal_predict_colnames(object$nonlocal_parsed$term_coef_name)
  as.list(as.data.frame(values))
}

# Bind the prediction columns `cols` onto the user's rows (`newdata`, or the
# fitted data) and return the data frame or the `return_tmb_object` list.
predict_return <- function(cols, req, object, r, obj, tmb_data, return_tmb_object) {
  nd <- if (is.null(req$newdata)) object$data else req$newdata
  for (col in names(cols)) nd[[col]] <- cols[[col]]
  nd[["_sdmTMB_time"]] <- NULL
  row.names(nd) <- NULL

  if (return_tmb_object) {
    return(list(data = nd, report = r, obj = obj, fit_obj = object, pred_tmb_data = tmb_data))
  }
  nd
}

# Link-scale population predictions for visreg, returned like predict.lm():
# a vector or `list(fit, se.fit)`. visreg_delta() sets `visreg_model` to pick
# the component; otherwise component 1.
predict_visreg <- function(object, newdata, se_fit) {
  if (!is.null(newdata) && !object$time %in% names(newdata)) {
    newdata[[object$time]] <- max(object$data[[object$time]], na.rm = TRUE)
  }
  model <- if ("visreg_model" %in% names(object)) object$visreg_model else 1L
  out <- predict_sdmTMB(object, newdata = newdata, se_fit = se_fit,
    re_form = NA, model = model)
  y_i <- object$tmb_data$y_i
  if (model == 2L && nrow(out) == nrow(y_i)) {
    out <- out[!is.na(y_i[, 2]), , drop = FALSE] # drop NAs from delta positive component
  }
  if (se_fit) list(fit = out$est, se.fit = out$est_se) else out$est
}

# Resolve all predict.sdmTMB() options once. Nothing downstream should modify
# the returned list.
predict_request <- function(object, newdata, type, se_fit, re_form,
                            re_form_iid, allow_new_levels, nsim, model,
                            offset, mcmc_samples, sims_var, sample_fe,
                            return_tmb_data, nonlocal_newdata) {
  if ("version" %in% names(object)) {
    check_sdmTMB_version(object$version)
  } else {
    cli_abort(c("This looks like a very old version of a model fit.",
        "Re-fit the model before predicting with it."))
  }
  family_spec <- .object_family_spec(object, caller = "`predict()`")
  multi_family <- family_spec$n_f > 1L
  has_two_components <- family_spec$n_m == 2L
  is_areal <- object$tmb_data$spatial_model %in% c(1L, 2L)
  if (!is_areal && is_areal_domain(object$spde)) {
    is_areal <- TRUE
  }
  xy_cols <- NULL
  if (!is_areal) {
    if (!"xy_cols" %in% names(object$spde)) {
      cli_warn(c("It looks like this model was fit with make_spde().",
      "Using `xy_cols`, but future versions of sdmTMB may not be compatible with this.",
      "Please replace make_spde() with make_mesh()."))
    } else {
      xy_cols <- object$spde$xy_cols
    }
  }

  assert_that(model[[1]] %in% c(NA, 1, 2),
    msg = "`model` argument not valid; should be one of NA, 1, 2")
  model <- model[[1]]
  type <- match.arg(type, c("link", "response"))
  if (!sims_var %in% c("est", names(predict_term_reports))) {
    cli_abort(c("`sims_var` must be one of {.val {c('est', names(predict_term_reports))}}.",
      "i" = "For other reported variables, use `return_tmb_report = TRUE` with `nsim`."))
  }
  if (sims_var != "est" && type == "response") {
    cli_abort("`type = 'response'` is only supported with `sims_var = 'est'`.")
  }
  if (!is.logical(sample_fe) || length(sample_fe) != 1L || is.na(sample_fe)) {
    cli_abort("`sample_fe` must be one non-missing logical value.")
  }
  if (!sample_fe && nsim > 0) {
    if (!is.null(mcmc_samples)) {
      cli_abort("`sample_fe = FALSE` cannot be combined with `mcmc_samples`.")
    }
    if (has_no_random_effects(object)) {
      cli_abort("`sample_fe = FALSE` requires a model with random effects.")
    }
  }
  if (isTRUE(se_fit) && type == "response") {
    cli_abort(c("Standard errors are only available on the link scale.",
      "i" = "Use `type = 'link'` with `se_fit = TRUE`, or use `nsim` for response-scale uncertainty."))
  }

  if (is.null(re_form) && isTRUE(se_fit)) {
    msg <- paste0("Prediction can be slow when `se_fit = TRUE` and random fields ",
      "are included (i.e., `re_form = NULL`). Consider using the `nsim` argument ",
      "to take draws from the joint precision matrix and summarizing the standard ",
      "deviation of those draws.")
    cli_inform(msg)
  }

  # Where the rows come from (`newdata_supplied`, which sets the default
  # offset) is separate from how they are predicted (`project`). Omitted
  # `newdata` uses the fitted report unless an option needs the projection
  # path, which then runs on the fitted data. Mixture families need the mean
  # adjustment only the projection path applies.
  newdata_supplied <- !is.null(newdata)
  project <- newdata_supplied || return_tmb_data ||
    !is.null(nonlocal_newdata) || multi_family || has_two_components ||
    any(object$family$family %in% rtmb_mixture_families) ||
    nsim > 0 || type == "response" || !is.null(mcmc_samples) || se_fit ||
    !is.null(re_form) || !is.null(re_form_iid) || !is.null(offset)
  if (project && !newdata_supplied) newdata <- object$data

  # from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  pop_pred_iid <- (!is.null(re_form_iid) && ((re_form_iid == ~0) || identical(re_form_iid, NA)))
  if (is.null(allow_new_levels)) {
    allow_new_levels <- pop_pred_iid
  }
  exclude_RE <- if (pop_pred_iid) 1L else object$tmb_data$exclude_RE

  named_list(
    family_spec, has_two_components, is_areal, xy_cols, model, type, se_fit,
    pop_pred, pop_pred_iid, allow_new_levels, exclude_RE, newdata,
    newdata_supplied, project, offset, nsim, sims_var, sample_fe, mcmc_samples
  )
}

# Offset rule: an explicit `offset` is used as is; otherwise omitted `newdata`
# carries over the fitted offset, and supplied `newdata` gets 0.
predict_offset <- function(req, object, n) {
  if (req$newdata_supplied && is.null(req$offset) && !all(object$offset == 0)) { # #372
    cli_inform(c(
      "Fitted object contains an offset but the offset is `NULL` in `predict.sdmTMB()` and `newdata` were supplied.",
      "Prediction will proceed assuming the offset vector is 0 in the prediction.",
      "Specify an offset vector in `predict.sdmTMB()` to override this."))
  }
  if (!is.null(req$offset)) {
    if (n != length(req$offset))
      cli_abort("Prediction offset vector does not equal number of rows in prediction dataset.")
    return(req$offset)
  }
  if (req$newdata_supplied) rep(0, n) else object$tmb_data$offset_i
}

# Each model component on `scale` (`est1`, `est2`), plus `est` for `model`,
# from the full (`proj_eta`) or population (`proj_fe`) linear predictor.
predict_components <- function(r, req, tmb_data, scale, model = NA) {
  .family_spec_component_prediction_output(
    x = r[[if (req$pop_pred) "proj_fe" else "proj_eta"]],
    family_spec = req$family_spec,
    row_family_id = tmb_data$proj_family_id + 1L,
    type = scale,
    model = model,
    offset = tmb_data$proj_offset_i
  )
}

# Select the link-scale `est` from one report `r` (a point estimate, a draw,
# or sdreport() estimates or standard errors) without transforming it.
# Two-component models with `model = NA` use the combined report; otherwise
# component `model`, NA on rows where that component is inactive.
predict_select <- function(r, req, tmb_data, model) {
  if (req$has_two_components && is.na(model)) {
    return(as.numeric(r[[if (req$pop_pred) "proj_fe_combined" else "proj_eta_combined"]]))
  }
  x <- as.matrix(r[[if (req$pop_pred) "proj_fe" else "proj_eta"]])
  m <- if (is.na(model)) 1L else as.integer(model)
  if (m > req$family_spec$n_m) return(rep(NA_real_, nrow(x)))
  out <- x[, m]
  if (m == 2L) {
    active <- .family_spec_component_active(req$family_spec, tmb_data$proj_family_id + 1L)
    out[!active[, 2L]] <- NA_real_
  }
  out
}

# The prediction `est` from one report `r` on `scale`. Response-scale
# components need both link predictors (e.g., Poisson-link delta models).
predict_est <- function(r, req, tmb_data, scale, model) {
  if (scale == "link") return(predict_select(r, req, tmb_data, model))
  if (req$has_two_components && is.na(model)) {
    return(as.numeric(r$proj_response_combined)) # C++ applies `pop_pred`
  }
  predict_components(r, req, tmb_data, "response", model)$est
}

# The fitted report's terms in the shapes predict_terms() expects, for the
# fast path (single-component, non-mixture models only). As in the projection
# path, `est_rf` holds the spatial, spatiotemporal, and SVC terms, and
# `est_non_rf` everything else (including the offset).
predict_fitted_term_reports <- function(r, object) {
  est <- r$eta_i[, 1]
  est_rf <- r$omega_s_A[, 1] + r$epsilon_st_A_vec[, 1]
  z_i <- unname(object$tmb_data$z_i)
  for (z in seq_len(ncol(z_i))) est_rf <- est_rf + r$zeta_s_A[, z, 1] * z_i[, z]
  list(
    est_non_rf = as.matrix(est - est_rf),
    est_rf = as.matrix(est_rf),
    omega_s = r$omega_s_A,
    zeta_s = r$zeta_s_A,
    epsilon_st = r$epsilon_st_A_vec
  )
}

# Linear predictor term columns and the projected report holding each.
predict_term_reports <- c(
  est_non_rf = "proj_fe",
  est_rf = "proj_rf",
  omega_s = "proj_omega_s_A",
  zeta_s = "proj_zeta_s_A",
  epsilon_st = "proj_epsilon_st_A_vec"
)

# Term columns from `x`, a list named like `predict_term_reports`
# with one column per component (`zeta_s`: rows x SVCs x components). Names get
# a component suffix only in two-component models, where component 2 is NA on
# rows whose family has no second component. Effects the model doesn't have
# are left out.
predict_terms <- function(x, object, family_spec, row_family_id) {
  d <- object$tmb_data
  n_m <- family_spec$n_m
  active2 <- if (n_m == 2L) {
    .family_spec_component_active(family_spec, row_family_id)[, 2L]
  }
  out <- list()
  for (col in names(x)) {
    svc <- col == "zeta_s"
    names_j <- if (svc) {
      paste0("zeta_s_", object$spatial_varying, recycle0 = TRUE)
    } else {
      col
    }
    for (j in seq_along(names_j)) {
      for (m in seq_len(n_m)) {
        include <- switch(col,
          omega_s = as.logical(d$include_spatial[m]),
          epsilon_st = !as.logical(d$spatial_only[m]),
          zeta_s = TRUE,
          n_m == 2L || !as.logical(d$no_spatial) # est_non_rf, est_rf
        )
        if (!include) next
        values <- if (svc) x[[col]][, j, m] else x[[col]][, m]
        if (m == 2L) values[!active2] <- NA_real_
        out[[if (n_m == 2L) paste0(names_j[j], m) else names_j[j]]] <- values
      }
    }
  }
  out
}

# Reports for each parameter draw: MCMC samples or draws from the joint
# precision matrix (or the random effects only, conditional on the fixed
# effects, if `sample_fe = FALSE`).
predict_draw_reports <- function(object, obj, lp, req) {
  nsim <- req$nsim
  if (!is.null(req$mcmc_samples)) {
    t_draws <- req$mcmc_samples
    if (nsim > 0) {
      if (nsim > ncol(t_draws)) {
        cli_abort("`nsim` must be <= number of MCMC samples.")
      }
      t_draws <- t_draws[, seq(ncol(t_draws) - nsim + 1, ncol(t_draws)), drop = FALSE]
    }
  } else {
    if (!"jointPrecision" %in% names(object$sd_report) && !has_no_random_effects(object)) {
      message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
      reinitialize(object)
      sd_report <- sdreport_sdmTMB(object$tmb_obj, getJointPrecision = TRUE)
    } else {
      sd_report <- object$sd_report
    }
    t_draws <- project_historical_draws(lp, sd_report, nsim,
      sample_fe = req$sample_fe, sample_historical_re = TRUE)
  }
  apply(t_draws, 2L, obj$report)
}

# Matrix of draws (rows x draws) of `sims_var` from the per-draw reports `r`.
predict_draws <- function(r, req, object, tmb_data, nd) {
  sims_var <- req$sims_var
  pred_row_family_id <- tmb_data$proj_family_id + 1L
  if (sims_var == "est") {
    out <- lapply(r, predict_est, req = req, tmb_data = tmb_data,
      scale = req$type, model = req$model)
    out <- do.call("cbind", out)
    rownames(out) <- nd[[object$time]] # for use in index calcs
    attr(out, "time") <- object$time
    attr(out, "link") <- if (req$type == "response") {
      "response"
    } else {
      .family_spec_prediction_link_name(
        family_spec = req$family_spec,
        row_family_id = pred_row_family_id,
        model = req$model
      )
    }
    return(out)
  }

  if (req$has_two_components && is.na(req$model)) {
    cli_warn("`model` argument was left as NA; defaulting to 1st model component.")
  }
  m <- if (req$has_two_components && !is.na(req$model)) as.integer(req$model) else 1L
  cols <- if (sims_var == "zeta_s") {
    if (length(object$spatial_varying)) paste0("zeta_s_", object$spatial_varying)
  } else {
    sims_var
  }
  if (req$has_two_components) cols <- paste0(cols, m)
  draws <- lapply(r, function(x) {
    predict_terms(
      stats::setNames(list(x[[predict_term_reports[[sims_var]]]]), sims_var),
      object, req$family_spec, pred_row_family_id
    )
  })
  if (!length(cols) || !all(cols %in% names(draws[[1]]))) {
    cli_abort("This model has no {.val {sims_var}} term{if (req$has_two_components) paste0(' in component ', m)} to draw from.")
  }
  out <- lapply(cols, function(col) do.call("cbind", lapply(draws, `[[`, col)))
  if (sims_var == "zeta_s") names(out) <- object$spatial_varying
  if (length(out) == 1L) out[[1]] else out
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

check_sdmTMB_version <- function(version) {
  fitted <- as.character(version)
  installed <- as.character(utils::packageVersion("sdmTMB"))
  if (fitted != installed) {
    cli_warn(c(
      "This model was fit with sdmTMB {fitted}, but sdmTMB {installed} is installed.",
      "i" = paste(
        "The TMB model may differ between versions, so prediction may fail",
        "or give different results. We recommend refitting the model with the installed version."
      )
    ), .frequency = "once", .frequency_id = paste0("sdmTMB_version_", fitted))
  }
}

check_time_class <- function(object, newdata) {
  cls1 <- class(object$data[[object$time]])
  cls2 <- class(newdata[[object$time]])
  if (!identical(cls1, cls2)) {
    if (!identical(sort(c(cls1, cls2)), c("integer", "numeric"))) {
      msg <- paste0(
        "Class of fitted time column (", cls1, ") does not match class of ",
        "`newdata` time column (", cls2 ,")."
        )
      cli_abort(msg)
    }
  }
}

# Is predict() being called by visreg (other than via residuals())?
# Match function names only; deparsing whole calls is slow when arguments
# are large objects (e.g., from `do.call()`).
is_visreg_call <- function() {
  fns <- vapply(sys.calls(), call_fn_name, character(1L))
  visreg_fns <- c("setupV", "visregPred", "build_visreg", "build_visreg2d", "visreg_pred")
  any(fns %in% visreg_fns) && !any(fns %in% "residuals")
}

# Name of the function in a call (`f` for `f()` or `pkg::f()`), else NA.
call_fn_name <- function(x) {
  f <- x[[1L]]
  if (is.call(f) && as.character(f[[1L]]) %in% c("::", ":::")) f <- f[[3L]]
  if (is.name(f)) as.character(f) else NA_character_
}
