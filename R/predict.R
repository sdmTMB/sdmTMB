# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Make predictions from a model fitted with [sdmTMB()], either for the fitted
#' data or for new data (e.g., a prediction grid). At new locations, the random
#' fields are interpolated from their estimated values at the mesh vertices.
#' Besides the overall prediction (`est`), the output separates the linear
#' predictor into the contribution of the spatial and spatiotemporal random
#' fields and everything else (see Value).
#'
#' @param object A model fitted with [sdmTMB()].
#' @param newdata A data frame to predict on. If `NULL` (default), predictions
#'   are for the fitted data. Must contain the columns used in the model
#'   formulas (including `spatial_varying` and `time_varying`), the coordinate
#'   columns used to build `mesh`, and, for spatiotemporal models, the `time`
#'   column. Time values must be among those in the fitted data or
#'   `extra_time`. Factor levels must have been seen in fitting, except for
#'   random effect grouping factors (see `allow_new_levels`).
#' @param type The scale of `est`: `"link"` (default; the linear predictor) or
#'   `"response"` (the expected value of the response). For delta models,
#'   `"response"` gives the expected value combining both components. Only
#'   `est`, `est1`, and `est2` depend on `type`; the other columns are always on
#'   the link scale. Standard errors (`se_fit`) require `"link"`.
#' @param se_fit Logical: calculate standard errors of `est` on the link scale
#'   (returned as `est_se`)? These account for uncertainty in both fixed and
#'   random effects but can be slow for many rows when random fields are
#'   included. For faster uncertainty, exclude the random fields with
#'   `re_form = NA` or summarize draws from `nsim`. Use `nsim` for
#'   response-scale uncertainty.
#' @param re_form Include the spatial and spatiotemporal random fields
#'   (including spatially varying coefficient fields)? `NULL` (default)
#'   includes them. `NA` or `~ 0` excludes them for population-level
#'   predictions, so that `est` equals `est_non_rf`; the random field columns
#'   are then omitted. Often used with `se_fit = TRUE` to plot covariate
#'   effects. IID random effects are set separately with `re_form_iid`.
#'   Both also apply to derived quantities such as [get_index()] when passed
#'   through its `predict_args`.
#' @param re_form_iid Include the IID random intercepts and slopes (e.g.,
#'   `(1 | g)` in `formula`)? `NULL` (default) includes them. `NA` or `~ 0`
#'   sets them all to zero. Excluding only some of them is not yet supported.
#' @param allow_new_levels Allow levels of random effect grouping factors in
#'   `newdata` that were not seen in fitting? A new level gets a random effect
#'   value of zero, while other grouping factors keep their random effects for
#'   those rows. `TRUE` allows them silently, `FALSE` gives an error (as with
#'   `allow.new.levels` in \pkg{lme4} and \pkg{glmmTMB}), and `NULL` (default)
#'   allows them with a warning. With `nsim`, `"sample"` instead draws each new
#'   level's random effects from their estimated (multivariate) normal
#'   distribution in every draw, as with `sample_new_levels = "gaussian"` in
#'   \pkg{brms}, so the draws include variation among levels. Rows sharing a new
#'   level share its draw, and fitted levels keep their usual draws. To include
#'   this variation for all rows, give the grouping factor a level not in the
#'   fitted data (see Examples). Not relevant if `re_form_iid` excludes the
#'   random effects. Grouping columns in `newdata` must be factors.
#' @param nsim Number of draws. If `> 0`, returns a matrix of draws instead of
#'   a data frame (see Value). Each draw takes the fixed and random effects
#'   from their approximate joint (multivariate normal) distribution and
#'   computes `est` (or `sims_var`). Summarize across draws (e.g.,
#'   `apply(x, 1, sd)` or quantiles) for uncertainty on any scale, or carry
#'   the draws through to derived quantities. Usually much faster than
#'   `se_fit = TRUE` for models with random fields.
#' @param sims_var Which quantity to return when `nsim > 0`: `"est"`
#'   (default; on the scale set by `type`) or one of the link-scale columns
#'   `"est_non_rf"`, `"est_rf"`, `"omega_s"`, `"zeta_s"`, or `"epsilon_st"`
#'   (see Value). The model must include the term. For delta models, these
#'   come from the component set by `model` (default: the first). With more
#'   than one spatially varying coefficient, `"zeta_s"` returns a list of
#'   matrices, one per coefficient. For other quantities, use
#'   `return_tmb_report = TRUE`.
#' @param sample_fe Logical. When `nsim > 0`, draw the fixed effects and other
#'   estimated parameters along with the random effects (`TRUE`, default)? If
#'   `FALSE`, these are held at their estimates and only the random effects
#'   (random fields, IID random effects, time-varying coefficients, and smoother
#'   coefficients) are drawn, conditional on those estimates. This ignores
#'   parameter uncertainty and so typically gives narrower intervals. Fixed
#'   effects are held at their estimates even with REML. See also the same
#'   argument in [project()].
#' @param model For delta models, what `est` (and `est_se` or the draws from
#'   `nsim` or `mcmc_samples`) refers to: `NA` (default) combines both
#'   components, `1` gives the first (binary) component, and `2` the second
#'   (positive) component. The data frame output always also includes each
#'   component as `est1` and `est2`. Ignored for other models. See the
#'   [delta-model
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html).
#' @param offset A numeric vector of offset values, one per row of `newdata`,
#'   or a single value for all rows (not a column name). If `NULL` (default),
#'   predictions for the fitted data
#'   use the fitted offset, and predictions with `newdata` use an offset of 0,
#'   i.e., predictions per unit of the offset (e.g., density rather than
#'   catch when the offset is log area swept).
#' @param mcmc_samples A matrix of posterior samples from a model passed to
#'   \pkg{tmbstan} (see `bayesian` in [sdmTMB()]), as returned by
#'   `sdmTMBextra::extract_mcmc()` in the
#'   \href{https://github.com/sdmTMB/sdmTMBextra}{sdmTMBextra} package. If
#'   supplied, returns a matrix of posterior draws in the same form as with
#'   `nsim`. If `nsim` is also supplied, the last `nsim` samples are used. See
#'   the \href{https://sdmTMB.github.io/sdmTMB/articles/bayesian.html}{Bayesian
#'   vignette}.
#' @param nonlocal_newdata An optional data frame of the `nonlocal_formula`
#'   covariates to predict with instead of those used in fitting (e.g., for a
#'   scenario with different conditions). Same requirements as `nonlocal_data`
#'   in [sdmTMB()]. The rows of `newdata` (or the fitted data) still set where
#'   and when predictions are made; this argument only supplies the covariate
#'   values that are diffused or lagged. If `NULL` (default), the covariates
#'   from `nonlocal_data` in [sdmTMB()] are reused if supplied (so `newdata`
#'   need not contain those columns); otherwise they come from `newdata`.
#' @param return_tmb_object `r lifecycle::badge("deprecated")` Logical. If
#'   `TRUE`, include the TMB object in a list-format output. Instead, pass the
#'   fitted model and `newdata` directly to [get_index()] or [get_cog()].
#' @param return_tmb_report Logical: return the TMB report (a list of all
#'   reported quantities at the estimated parameters) instead of a data frame?
#'   With `nsim > 0` or `mcmc_samples`, a list with one report per draw.
#'   Mainly for developers.
#' @param return_tmb_data Logical: return the data list passed to TMB instead
#'   of predicting? Used internally.
#' @param ... Unused.
#'
#' @return
#' By default, `newdata` (or the fitted data) with these columns added:
#'
#' * `est`: The prediction on the scale set by `type`.
#' * `est_se`: The standard error of `est` on the link scale, if
#'   `se_fit = TRUE`.
#' * `est_non_rf`: The linear predictor excluding the spatial and
#'   spatiotemporal random fields: fixed effects, smoothers, IID random
#'   effects, time-varying coefficients, and the offset.
#' * `est_rf`: The sum of all random field terms, including spatially varying
#'   coefficient fields multiplied by their covariates. On the link scale,
#'   `est_non_rf + est_rf` equals `est`.
#' * `omega_s`: The spatial random field.
#' * `zeta_s_<x>`: The spatially varying coefficient field for covariate `<x>`
#'   in `spatial_varying`: the local deviation from the average coefficient,
#'   not multiplied by the covariate.
#' * `epsilon_st`: The spatiotemporal random field.
#' * `nl_*`: The transformed covariate values for each `nonlocal_formula`
#'   term.
#'
#' Columns for terms not in the model are left out, as are the random field
#' columns with `re_form = NA`.
#'
#' Delta models instead return each column (other than `est` and `est_se`)
#' per component, with suffixes `1` and `2` for the first and second
#' components (e.g., `est1`, `est2`, `omega_s1`, `omega_s2`). The identity
#' above then holds within each component, and `est` combines the components
#' (or gives one of them; see `model`).
#'
#' If `nsim > 0` or `mcmc_samples` is supplied: a matrix with one row per row
#' of `newdata` (or the fitted data) and one column per draw. Row names are
#' the time values.
#'
#' If `return_tmb_object = TRUE` (deprecated): a list with elements `data`
#' (the data frame above), `report` (the TMB report), `obj` (the TMB object
#' from the prediction), and `fit_obj` (the fitted model).
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
#'   ggtitle("Prediction without random fields (fixed effects only here)") +
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
#'   ggtitle("Prediction without random fields (fixed effects only here)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' plot_map(p, exp(est)) +
#'   ggtitle("Prediction (fixed effects + all random effects)") +
#'   scale_fill_viridis_c(trans = "sqrt")
#'
#' # Including variation among random effect levels ----------------------
#'
#' d <- pcod_2011
#' d$fyear <- as.factor(d$year)
#' m <- sdmTMB(density ~ depth_scaled + (1 | fyear), data = d,
#'   family = tweedie(), spatial = "off")
#' # a level not in the data draws a new year effect in each draw:
#' nd <- data.frame(depth_scaled = 0, fyear = factor("new"))
#' p <- predict(m, newdata = nd, nsim = 200, type = "response",
#'   allow_new_levels = "sample")
#' quantile(p, c(0.1, 0.5, 0.9))
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
    sample_fe = sample_fe
  )

  tmb_data <- object$tmb_data
  if (is.null(tmb_data$link_pred) && !is.null(tmb_data$link)) {
    tmb_data$link_pred <- tmb_data$link
  }
  tmb_data$do_predict <- 1L
  has_nonlocal <- !is.null(object$nonlocal_parsed)

  prep <- predict_prepare(object, req, tmb_data, nonlocal_newdata)
  tmb_data <- prep$tmb_data
  if (return_tmb_data) {
    return(tmb_data)
  }

  # On RTMB, reports alone need no AD object, which for large grids dominates
  # prediction time and memory.
  lp <- object$last.par.best
  new_tmb_obj <- if (!predict_report_only(object, req, return_tmb_object)) {
    predict_objective(object, tmb_data,
      length(unlist(lapply(prep$new_re, `[[`, "rows"))))
  }
  if (req$nsim > 0 || !is.null(req$mcmc_samples)) {
    r <- predict_draw_reports(object, new_tmb_obj, tmb_data, lp, req, prep$new_re)
    if (return_tmb_report) return(r)
    return(predict_draws(r, req, object, tmb_data, prep$nd))
  }
  r <- predict_report(new_tmb_obj, tmb_data, object$parlist, lp)
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

  predict_return(cols, req, object, r, obj, tmb_data, return_tmb_object)
}

# Whether predictions need only the model's reports at saved parameter values
# (the fitted ones or draws), which the RTMB backend can evaluate without an
# AD object.
predict_report_only <- function(object, req, return_tmb_object) {
  backend_sdmTMB(object) == "rtmb" && !req$se_fit && !return_tmb_object
}

# Reports at the parameter list `pars` (equivalently, the parameter vector
# `x`): from the AD object `obj`, or with plain numbers if `obj` is NULL.
# Prediction needs no observation likelihood.
predict_report <- function(obj, tmb_data, pars, x) {
  if (is.null(obj)) {
    rtmb_report_values(tmb_data, pars, likelihood = FALSE)
  } else {
    obj$report(x)
  }
}

# Build the prediction objective from `tmb_data`, with `n_new_re` extra
# `re_b_pars` rows for sampled new levels, at the saved fitted parameters.
predict_objective <- function(object, tmb_data, n_new_re = 0L) {
  parameters <- object$parlist
  map <- object$tmb_map
  if (n_new_re > 0L) {
    # Rows for sampled new levels go after the fitted ones and are numbered
    # after them in the map, so they sit at the end of the `re_b_pars` block.
    b_map <- matrix(as.integer(map$re_b_pars), ncol = ncol(parameters$re_b_pars))
    new_map <- matrix(NA_integer_, n_new_re, ncol(b_map))
    has_re <- object$tmb_data$n_re_groups > 0
    new_map[, has_re] <- max(b_map, na.rm = TRUE) + seq_len(n_new_re * sum(has_re))
    map$re_b_pars <- factor(rbind(b_map, new_map))
    parameters$re_b_pars <- rbind(parameters$re_b_pars,
      matrix(0, n_new_re, ncol(b_map)))
  }
  obj <- make_sdmTMB_adfun(
    data = tmb_data,
    profile = object$control$profile,
    parameters = parameters,
    map = map,
    random = object$tmb_random,
    backend = backend_sdmTMB(object),
    silent = TRUE
  )
  obj
}

# Draws from N(0, Sigma) for the sampled new levels in `new_re`, given
# `re_cov_pars`: a matrix with a row per new-level row of `re_b_pars` and a
# column per component (zero for components without random intercepts).
predict_new_re_b <- function(object, re_cov_pars, new_re) {
  rows <- unlist(lapply(new_re, `[[`, "rows"))
  n_re_groups <- object$tmb_data$n_re_groups
  b <- matrix(0, length(rows), length(n_re_groups))
  for (m in which(n_re_groups > 0)) {
    groups <- rtmb_re_groups(object$tmb_data, m)
    for (term in new_re) {
      L <- t(chol(iid_re_covariance(re_cov_pars, groups[[term$term]], m)))
      # one column per level:
      z <- matrix(stats::rnorm(nrow(L) * term$n_levels), nrow(L))
      b[match(term$rows, rows), m] <- L %*% z
    }
  }
  b
}

# Nonlocal covariate values as named prediction columns.
predict_nonlocal_cols <- function(object, values) {
  colnames(values) <- .nonlocal_predict_colnames(object$nonlocal_parsed$term_coef_name)
  as.list(as.data.frame(values))
}

# Bind the prediction columns `cols` onto the user's rows (`newdata`, or the
# fitted data) and return the data frame or the `return_tmb_object` list.
predict_return <- function(cols, req, object, r, obj, tmb_data, return_tmb_object) {
  nd <- req$newdata
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
                            offset, mcmc_samples, sims_var, sample_fe) {
  check_fit_version(object)
  check_sdmTMB_version(object$version)
  family_spec <- .object_family_spec(object, caller = "`predict()`")
  has_two_components <- family_spec$n_m == 2L
  is_areal <- is_areal_fit(object)
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

  if (identical(allow_new_levels, "sample") && (se_fit ||
    (nsim == 0 && is.null(mcmc_samples)))) {
    cli_abort(c("`allow_new_levels = \"sample\"` requires draws.",
      "i" = "Use `nsim` (or `mcmc_samples`), and summarize the draws for uncertainty."))
  }

  if (is.null(re_form) && isTRUE(se_fit)) {
    msg <- paste0("Prediction can be slow when `se_fit = TRUE` and random fields ",
      "are included (i.e., `re_form = NULL`). Consider using the `nsim` argument ",
      "to take draws from the joint precision matrix and summarizing the standard ",
      "deviation of those draws.")
    cli_inform(msg)
  }

  # Omitted `newdata` predicts on the fitted data. `newdata_supplied` only
  # sets the default offset (see predict_offset()).
  newdata_supplied <- !is.null(newdata)
  if (!newdata_supplied) newdata <- object$data

  # from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  pop_pred_iid <- (!is.null(re_form_iid) && ((re_form_iid == ~0) || identical(re_form_iid, NA)))
  if (!is.null(allow_new_levels) && !isTRUE(allow_new_levels) &&
    !isFALSE(allow_new_levels) && !identical(allow_new_levels, "sample")) {
    cli_abort("`allow_new_levels` must be `NULL`, `TRUE`, `FALSE`, or `\"sample\"`.")
  }
  exclude_RE <- if (pop_pred_iid) 1L else object$tmb_data$exclude_RE

  named_list(
    family_spec, has_two_components, is_areal, xy_cols, model, type, se_fit,
    pop_pred, pop_pred_iid, allow_new_levels, exclude_RE, newdata,
    newdata_supplied, offset, nsim, sims_var, sample_fe, mcmc_samples,
    ln_phi = .object_par(object, "ln_phi"), psi = .object_par(object, "psi")
  )
}

# Offset rule: an explicit `offset` is used as is; otherwise omitted `newdata`
# carries over the fitted offset, and supplied `newdata` gets 0.
predict_offset <- function(req, object, n) {
  if (req$newdata_supplied && is.null(req$offset) && !all(object$offset == 0)) { # #372
    cli_inform(c(
      "i" = "This model has an offset, but no {.arg offset} was given for {.arg newdata}, so an offset of 0 is used.",
      " " = "With a log offset (e.g., log of effort, area, or time), predictions are then per unit of that quantity rather than totals.",
      " " = "That is usually what you want for maps and indices. Supply {.arg offset} for other values, or use {.code offset = 0} to hide this message."
    ), .frequency = "once", .frequency_id = "sdmTMB_predict_offset")
  }
  if (!is.null(req$offset)) {
    if (!is.numeric(req$offset) && !all(is.na(req$offset)))
      cli_abort("`offset` must be a numeric vector (not a column name).")
    if (length(req$offset) == 1L) return(rep(req$offset, n))
    if (n != length(req$offset))
      cli_abort("`offset` must have length one or one value per row of the prediction data.")
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
    offset = tmb_data$proj_offset_i,
    ln_phi = req$ln_phi,
    psi = req$psi
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
          est_non_rf = TRUE,
          n_m == 2L || !as.logical(d$no_spatial) # est_rf
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
# effects, if `sample_fe = FALSE`), with any sampled new levels. Draws are
# in the fitted parameter layout; `obj` (NULL on RTMB, see predict_report())
# has the new-level rows after the fitted ones in `re_b_pars`.
predict_draw_reports <- function(object, obj, tmb_data, lp, req, new_re = list()) {
  t_draws <- predict_parameter_draws(object, lp, req)
  last_b <- max(which(names(lp) == "re_b_pars"), 0L)
  has_re <- object$tmb_data$n_re_groups > 0
  pars <- lapply(seq_len(ncol(t_draws)), function(i) {
    object$tmb_obj$env$parList(par = t_draws[, i])
  })
  out <- lapply(seq_along(pars), function(i) {
    x <- t_draws[, i]
    p <- pars[[i]]
    if (length(new_re)) {
      new_b <- predict_new_re_b(object, p$re_cov_pars, new_re)
      p$re_b_pars <- rbind(p$re_b_pars, new_b)
      x <- append(x, new_b[, has_re], after = last_b)
    }
    predict_report(obj, tmb_data, p, x)
  })
  # Response-scale means also depend on `ln_phi` and `psi`, which vary by draw:
  attr(out, "draw_pars") <- lapply(pars, function(p) {
    p[intersect(c("ln_phi", "psi"), names(p))]
  })
  out
}

# Parameter draws (parameters x draws) for predict_draw_reports().
predict_parameter_draws <- function(object, lp, req) {
  nsim <- req$nsim
  if (!is.null(req$mcmc_samples)) {
    t_draws <- req$mcmc_samples
    if (nsim > 0) {
      if (nsim > ncol(t_draws)) {
        cli_abort("`nsim` must be <= number of MCMC samples.")
      }
      t_draws <- t_draws[, seq(ncol(t_draws) - nsim + 1, ncol(t_draws)), drop = FALSE]
    }
    return(t_draws)
  }
  if (!"jointPrecision" %in% names(object$sd_report) && !has_no_random_effects(object)) {
    message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
    reinitialize(object)
    sd_report <- sdreport_sdmTMB(object$tmb_obj, getJointPrecision = TRUE)
  } else {
    sd_report <- object$sd_report
  }
  project_historical_draws(lp, sd_report, nsim,
    sample_fe = req$sample_fe, sample_historical_re = TRUE)
}

# Matrix of draws (rows x draws) of `sims_var` from the per-draw reports `r`.
predict_draws <- function(r, req, object, tmb_data, nd) {
  sims_var <- req$sims_var
  pred_row_family_id <- tmb_data$proj_family_id + 1L
  if (sims_var == "est") {
    draw_pars <- attr(r, "draw_pars")
    out <- lapply(seq_along(r), function(i) {
      req_i <- req
      req_i[names(draw_pars[[i]])] <- draw_pars[[i]]
      predict_est(r[[i]], req_i, tmb_data, scale = req$type, model = req$model)
    })
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
        "The TMB model may differ between versions, so using this model may fail",
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
