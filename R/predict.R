# @importFrom stats predict
# @rdname predict
# @export

#' Predict from an sdmTMB model
#'
#' Make predictions from an \pkg{sdmTMB} model; can predict on the original or
#' new data.
#'
#' @param object A model fitted with [sdmTMB()].
#' @param newdata A data frame to make predictions on. This should be a data
#'   frame with the same predictor columns as in the fitted data and a time
#'   column (if this is a spatiotemporal model) with the same name as in the
#'   fitted data.
#' @param type Should the `est` column be in link (default) or response space?
#' @param se_fit Should standard errors on predictions be calculated? Warning:
#'   can be slow for large datasets or high-resolution projections when random
#'   fields are included. For faster uncertainty estimation, either use
#'   `re_form = NA` to exclude random fields or use the `nsim` argument to
#'   simulate from the joint precision matrix.
#' @param return_tmb_object Logical. If `TRUE`, will include the TMB object in a
#'   list format output. Necessary for the [get_index()] or [get_cog()]
#'   functions.
#' @param re_form `NULL` to include all spatial/spatiotemporal random fields in
#'   predictions. `~0` or `NA` for population-level predictions (predictions
#'   from fixed effects only, marginalizing over random fields). Often used with
#'   `se_fit = TRUE` to visualize marginal effects. Does not affect
#'   [get_index()] calculations.
#' @param re_form_iid `NULL` to specify including all random intercepts in the
#'   predictions. `~0` or `NA` for population-level predictions. No other
#'   options (e.g., some but not all random intercepts) are implemented yet.
#'   Only affects predictions with `newdata`. This *does* affects [get_index()].
#' @param nsim If `> 0`, simulate from the joint precision matrix with `nsim`
#'   draws. Returns a matrix of `nrow(newdata)` by `nsim` with each column
#'   representing one draw of the linear predictor (in link space). Simulating
#'   from the joint precision matrix accounts for uncertainty in both fixed and
#'   random effects. Use this to derive uncertainty on predictions (e.g.,
#'   `apply(x, 1, sd)`) or propagate uncertainty to derived quantities. This is
#'   the fastest way to characterize spatial uncertainty with sdmTMB.
#' @param sims_var Experimental: Which TMB reported variable from the model
#'   should be extracted from the joint precision matrix simulation draws?
#'   Defaults to link-space predictions. Options include: `"omega_s"`,
#'   `"zeta_s"`, `"epsilon_st"`, and `"est_rf"` (as described below).
#'   Other options will be passed verbatim.
#' @param mcmc_samples See `extract_mcmc()` in the
#'   \href{https://github.com/sdmTMB/sdmTMBextra}{sdmTMBextra} package for
#'   more details and the
#'   \href{https://sdmTMB.github.io/sdmTMB/articles/bayesian.html}{Bayesian vignette}.
#'   If specified, the predict function will return a matrix of a similar form
#'   as if `nsim > 0` but representing Bayesian posterior samples from the Stan
#'   model.
#' @param model Which component to predict from delta/hurdle models when `nsim >
#'   0` or `mcmc_samples` is supplied. `NA` (default) returns the combined
#'   prediction from both components; `1` returns the binomial component only; `2`
#'   returns the positive component only. Predictions are on the link or response
#'   scale depending on `type`. For regular predictions (without simulation),
#'   both components are returned. See the [delta-model
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html).
#' @param offset A numeric vector of optional offset values. If left at default
#'   `NULL`, the offset is implicitly left at 0.
#' @param return_tmb_report Logical: return the output from the TMB
#'   report? For regular prediction, this is all the reported variables
#'   at the MLE parameter values. For `nsim > 0` or when `mcmc_samples`
#'   is supplied, this is a list where each element is a sample and the
#'   contents of each element is the output of the report for that sample.
#' @param return_tmb_data Logical: return formatted data for TMB? Used
#'   internally.
#' @param ... Not implemented.
#'
#' @return
#' If `return_tmb_object = FALSE` (and `nsim = 0` and `mcmc_samples = NULL`):
#'
#' A data frame:
#' * `est`: Estimate in link space (everything included)
#' * `est_non_rf`: Estimate from everything except random fields (fixed effects, random intercepts, time-varying effects, etc.)
#' * `est_rf`: Estimate from all random fields combined
#' * `omega_s`: Spatial random field (models consistent spatial patterns)
#' * `zeta_s`: Spatially varying coefficient field (models how effects vary across space)
#' * `epsilon_st`: Spatiotemporal random field (models spatial patterns that vary over time)
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
#' * Rows represent predictions with one row per row of `newdata`
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
#' qqnorm(predictions$resids);abline(a = 0, b = 1)
#'
#' # Predictions onto new data --------------------------------------------
#'
#' qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
#' predictions <- predict(m, newdata = qcs_grid_2011)
#'
#' \donttest{
#' # A short function for plotting our predictions:
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
#' # element. If time isn't also a fixed effect then it doesn't matter what you pick:
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
#'   spatiotemporal = "AR1", # using an AR1 to have something to forecast with
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
  nsim = 0,
  sims_var = "est",
  model = c(NA, 1, 2),
  offset = NULL,
  mcmc_samples = NULL,
  return_tmb_object = FALSE,
  return_tmb_report = FALSE,
  return_tmb_data = FALSE,
  ...) {

  if ("version" %in% names(object)) {
    check_sdmTMB_version(object$version)
  } else {
    cli_abort(c("This looks like a very old version of a model fit.",
        "Re-fit the model before predicting with it."))
  }
  if (!"xy_cols" %in% names(object$spde)) {
    cli_warn(c("It looks like this model was fit with make_spde().",
    "Using `xy_cols`, but future versions of sdmTMB may not be compatible with this.",
    "Please replace make_spde() with make_mesh()."))
  } else {
    xy_cols <- object$spde$xy_cols
  }

  if (object$version < numeric_version("0.5.0.9001")) {
    cli_abort("This model was fit with an older version of sdmTMB before internal handling of `extra_time` was simplified. Please refit your model before predicting on it (or install version 0.5.0 or 0.5.0.9000).")
  }

  area <- 1
  sims <- nsim

  reinitialize(object)

  assert_that(model[[1]] %in% c(NA, 1, 2),
    msg = "`model` argument not valid; should be one of NA, 1, 2")
  if (missing(model)) {
    if (.has_delta_attr(object)) model <- attr(object, "delta_model_predict") # for ggpredict
  }
  model <- model[[1]]
  type <- match.arg(type)
  # FIXME parallel setup here?

  if (is.null(re_form) && isTRUE(se_fit)) {
    msg <- paste0("Prediction can be slow when `se_fit = TRUE` and random fields ",
      "are included (i.e., `re_form = NA`). Consider using the `nsim` argument ",
      "to take draws from the joint precision matrix and summarizing the standard ",
      "devation of those draws.")
    cli_inform(msg)
  }

  # places where we force newdata:
  nd_arg_was_null <- FALSE
  if (is.null(newdata)) {
    if (is_delta(object) || nsim > 0 || type == "response" || !is.null(mcmc_samples) || se_fit || !is.null(re_form) || !is.null(re_form_iid) || !is.null(offset) || isTRUE(object$family$delta)) {
      newdata <- object$data
      nd_arg_was_null <- TRUE # will be used to carry over the offset
    }
  }
  sys_calls <- unlist(lapply(sys.calls(), deparse)) # retrieve function that called this
  vr <- check_visreg(sys_calls)
  visreg_df <- vr$visreg_df
  if (visreg_df) {
    re_form <- vr$re_form
    se_fit <- vr$se_fit
  }

  # from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))
  pop_pred_iid <- (!is.null(re_form_iid) && ((re_form_iid == ~0) || identical(re_form_iid, NA)))

  exclude_RE <- if (pop_pred_iid) 1L else object$tmb_data$exclude_RE

  tmb_data <- object$tmb_data
  tmb_data$do_predict <- 1L
  no_spatial <- as.logical(object$tmb_data$no_spatial)

  if (!is.null(newdata)) {
    if (any(!xy_cols %in% names(newdata)) && isFALSE(pop_pred) && !no_spatial)
      cli_abort(c("`xy_cols` (the column names for the x and y coordinates) are not in `newdata`.",
          "Did you miss specifying the argument `xy_cols` to match your data?",
          "The newer `make_mesh()` (vs. `make_spde()`) takes care of this for you."))

    if (isFALSE(pop_pred) && !no_spatial) {
      xy_orig <- object$data[,xy_cols]
      xy_nd <- newdata[,xy_cols]
      all_outside <- function(x1, x2) {
        min(x1) > max(x2) || max(x1) < min(x2)
      }
      if (all_outside(xy_orig[,1], xy_nd[,1]) || all_outside(xy_orig[,2], xy_nd[,2])) {
        cli_warn(c("`newdata` prediction coordinates appear to be ouside the fitted coordinates.",
          "This will likely cause all your random field values to be returned as 0.",
          "Check your coordinates including any conversions between projections.",
          "If working with UTMs, are both in km or m?"))
      }
    }

    if (object$time == "_sdmTMB_time") newdata[[object$time]] <- 0L
    if (visreg_df) {
      if (!object$time %in% names(newdata)) {
        newdata[[object$time]] <- max(object$data[[object$time]], na.rm = TRUE)
      }
    }

    check_time_class(object, newdata)
    original_time <- object$time_lu$time_from_data
    new_data_time <- unique(newdata[[object$time]])

    if (!all(new_data_time %in% original_time))
      cli_abort(c("Some new time elements were found in `newdata`. ",
        "If you would like to predict on new time elements,",
        "see the `extra_time` argument in `?sdmTMB`.")
      )

    # If making population predictions (with standard errors), we don't need
    # to worry about space, so fill in dummy values if the user hasn't made any:
    fake_spatial_added <- FALSE
    if (pop_pred) {
      for (i in c(1, 2)) {
        if (!xy_cols[[i]] %in% names(newdata)) {
          suppressWarnings({
            newdata[[xy_cols[[i]]]] <- mean(object$data[[xy_cols[[i]]]], na.rm = TRUE)
            fake_spatial_added <- TRUE
          })
        }
      }
    }

    if (sum(is.na(new_data_time)) > 0)
      cli_abort(c("There is at least one NA value in the time column.",
        "Please remove it."))

    newdata$sdm_orig_id <- seq(1L, nrow(newdata))

    if (!no_spatial) {
      if (requireNamespace("dplyr", quietly = TRUE)) { # faster
        unique_newdata <- dplyr::distinct(newdata[, xy_cols, drop = FALSE])
      } else {
        unique_newdata <- unique(newdata[, xy_cols, drop = FALSE])
      }
      unique_newdata[["sdm_spatial_id"]] <- seq(1, nrow(unique_newdata)) - 1L

      if (requireNamespace("dplyr", quietly = TRUE)) { # much faster
        newdata <- dplyr::left_join(newdata, unique_newdata, by = xy_cols)
      } else {
        newdata <- base::merge(newdata, unique_newdata, by = xy_cols,
          all.x = TRUE, all.y = FALSE)
        newdata <- newdata[order(newdata$sdm_orig_id),, drop = FALSE]
      }
      proj_mesh <- fmesher::fm_basis(object$spde$mesh,
        loc = as.matrix(unique_newdata[, xy_cols, drop = FALSE]))
    } else {
      proj_mesh <- object$spde$A_st # fake
      if (!all(object$spde$xy_cols %in% names(newdata))) {
        newdata[[xy_cols[1]]] <- NA_real_ # fake
        newdata[[xy_cols[2]]] <- NA_real_ # fake
      }
      newdata[["sdm_spatial_id"]] <- rep(0L, nrow(newdata)) # fake
    }

    if (length(object$formula) == 1L) {
      # this formula has breakpt() etc. in it:
      thresh <- list(check_and_parse_thresh_params(object$formula[[1]], newdata))
      formula <- list(thresh[[1]]$formula) # this one does not
    } else {
      thresh <- list(check_and_parse_thresh_params(object$formula[[1]], newdata),
        check_and_parse_thresh_params(object$formula[[2]], newdata))
      formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
    }

    nd <- newdata
    response <- get_response(object$formula[[1]])
    sdmTMB_fake_response <- FALSE
    if (!response %in% names(nd)) {
      nd[[response]] <- 0 # fake for model.matrix
      sdmTMB_fake_response <- TRUE
    }

    if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE

    # FIXME check if random slopes and intercepts are the same in both linear predictors?
    # parse random intercept/slope sparse model matrices on new data:
    Zt_list <- list()

    if (sum(object$tmb_data$n_re_groups) > 0 && isFALSE(pop_pred_iid)) {
      for (ii in seq_len(length(formula))) {
        xx <- parse_formula(remove_s_and_t2(object$smoothers$formula_no_sm), nd)
        # factor level checks:
        RE_names <- xx$barnames
        for (i in seq_along(RE_names)) {
          assert_that(is.factor(newdata[[RE_names[i]]]),
            msg = sprintf("Random effect group column `%s` in newdata is not a factor.", RE_names[i]))
          levels_fit <- levels(object$data[[RE_names[i]]])
          levels_nd <- levels(newdata[[RE_names[i]]])
          if (sum(!levels_nd %in% levels_fit)) {
            msg <- paste0("Extra levels found in random intercept factor levels for `", RE_names[i],
              "`. Please remove them.")
            cli_abort(msg)
          }
        }

        # now do with a joint data frame to ensure factor levels match
        common_cols <- intersect(colnames(object$data), colnames(nd))
        joint_df <- rbind(object$data[,common_cols,drop=FALSE], nd[,common_cols,drop=FALSE])
        xx <- parse_formula(object$smoothers$formula_no_sm, joint_df)
        # drop the original data:
        Zt <- xx$re_cov_terms$Zt[,seq(nrow(object$data) + 1, nrow(object$data) + nrow(nd))]
        Zt_list[[ii]] <- Zt
      }
    }

    # deal with prediction IID random intercepts:
    # RE_names <- object$split_formula[[1]]$barnames # TODO DELTA HARDCODED TO 1 here; fine for now

    ## not checking so that not all factors need to be in prediction:
    # fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
    # proj_RE_indexes <- vapply(RE_names, function(x) as.integer(nd[[x]]) - 1L, rep(1L, nrow(nd)))

    # if (isFALSE(pop_pred_iid)) {
    #   for (i in seq_along(RE_names)) {
    #     # checking newdata random intercept columns are factors
    #     assert_that(is.factor(newdata[[RE_names[i]]]),
    #                 msg = sprintf("Random effect group column `%s` in newdata is not a factor.", RE_names[i]))
    #     levels_fit <- levels(object$data[[RE_names[i]]])
    #     levels_nd <- levels(newdata[[RE_names[i]]])
    #     if (sum(!levels_nd %in% levels_fit)) {
    #       msg <- paste0("Extra levels found in random intercept factor levels for `", RE_names[i],
    #         "`. Please remove them.")
    #       cli_abort(msg)
    #     }
    #   }
    # }

    proj_X_ij <- list()
    for (i in seq_along(object$formula)) {
      f2 <- remove_s_and_t2(object$split_formula[[i]]$form_no_bars)#object$smoothers$formula_no_bars_no_sm
      tt <- stats::terms(f2)
      attr(tt, "predvars") <- attr(object$terms[[i]], "predvars")
      Terms <- stats::delete.response(tt)
      mf <- model.frame(Terms, newdata, xlev = object$xlevels[[i]])
      proj_X_ij[[i]] <- model.matrix(Terms, mf, contrasts.arg = object$contrasts[[i]])
    }

    # TODO DELTA hardcoded to 1:
    sm <- parse_smoothers(object$smoothers$formula_no_bars, data = object$data,
      newdata = nd, basis_prev = object$smoothers$basis_out)

    if (!is.null(object$time_varying))
      proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
    else
      proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy

    if (length(area) != nrow(proj_X_ij[[1]]) && length(area) != 1L) {
      cli_abort("`area` should be of the same length as `nrow(newdata)` or of length 1.")
    }

    # newdata, null offset in predict, and non-null in fit #372
    if (isFALSE(nd_arg_was_null) && is.null(offset) && !all(object$offset == 0)) {
      msg <- c(
        "Fitted object contains an offset but the offset is `NULL` in `predict.sdmTMB()` and `newdata` were supplied.",
        "Prediction will proceed assuming the offset vector is 0 in the prediction.",
        "Specify an offset vector in `predict.sdmTMB()` to override this.")
      cli_inform(msg)
    }

    if (!is.null(offset)) {
      if (nrow(proj_X_ij[[1]]) != length(offset))
        cli_abort("Prediction offset vector does not equal number of rows in prediction dataset.")
    }
    tmb_data$proj_offset_i <- if (!is.null(offset)) offset else rep(0, nrow(proj_X_ij[[1]]))
    # if (nd_arg_was_null) tmb_data$proj_offset_i <- tmb_data$offset_i
    tmb_data$proj_X_threshold <- thresh[[1]]$X_threshold # TODO DELTA HARDCODED TO 1
    tmb_data$area_i <- if (length(area) == 1L) rep(area, nrow(proj_X_ij[[1]])) else area
    tmb_data$proj_mesh <- proj_mesh
    tmb_data$proj_X_ij <- proj_X_ij
    tmb_data$proj_X_rw_ik <- proj_X_rw_ik
    # tmb_data$proj_RE_indexes <- proj_RE_indexes

    tmb_data$Zt_list_proj <- Zt_list
    time_lu <- object$time_lu
    tmb_data$proj_year <- time_lu$year_i[match(nd[[object$time]], time_lu$time_from_data)] # was make_year_i(nd[[object$time]])
    tmb_data$proj_time_include <- as.integer(time_lu$time_from_data %in% nd[[object$time]])
    tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
    tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
    tmb_data$calc_se <- as.integer(se_fit)
    tmb_data$pop_pred <- as.integer(pop_pred)
    tmb_data$exclude_RE <- exclude_RE
    tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
    tmb_data$proj_Zs <- sm$Zs
    tmb_data$proj_Xs <- sm$Xs

    # SVC:
    if (!is.null(object$spatial_varying)) {
      # recreate original data SVC formula stuff:
      z_i_orig <- model.matrix(object$spatial_varying_formula, object$data)
      svc_contrasts <- attr(z_i_orig, which = "contrasts")
      ttsv <- stats::terms(object$spatial_varying_formula)
      mfsv <- model.frame(ttsv, object$data)
      mtsv <- attr(mfsv, "terms")
      xlevelssv <- stats::.getXlevels(mtsv, mfsv)
      # apply it to prediction data:
      mfsv_new <- model.frame(ttsv, newdata, xlev = xlevelssv)
      z_i <- model.matrix(ttsv, mfsv_new, contrasts.arg = svc_contrasts)
      .int <- grep("(Intercept)", colnames(z_i))
      if (sum(.int) > 0) z_i <- z_i[,-.int,drop=FALSE]
    } else {
      z_i <- matrix(0, nrow(newdata), 0L)
    }
    tmb_data$proj_z_i <- z_i

    epsilon_covariate <- rep(0, length(unique(newdata[[object$time]])))
    if (tmb_data$est_epsilon_model) {
      # covariate vector dimensioned by number of time steps
      time_steps <- unique(newdata[[object$time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- newdata[newdata[[object$time]] == time_steps[i],
            object$epsilon_predictor, drop = TRUE][[1]]
      }
    }
    tmb_data$epsilon_predictor <- epsilon_covariate

    if (return_tmb_data) {
      return(tmb_data)
    }

    new_tmb_obj <- TMB::MakeADFun(
      data = tmb_data,
      profile = object$control$profile,
      parameters = get_pars(object),
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      silent = TRUE
    )

    old_par <- object$model$par
    # need to initialize the new TMB object once:
    new_tmb_obj$fn(old_par)

    if (sims > 0 && is.null(mcmc_samples)) {
      if (!"jointPrecision" %in% names(object$sd_report) && !has_no_random_effects(object)) {
        message("Rerunning TMB::sdreport() with `getJointPrecision = TRUE`.")
        sd_report <- TMB::sdreport(object$tmb_obj, getJointPrecision = TRUE)
      } else {
        sd_report <- object$sd_report
      }
      if (has_no_random_effects(object)) {
        t_draws <- t(mvtnorm::rmvnorm(n = sims, mean = sd_report$par.fixed,
          sigma = sd_report$cov.fixed))
        row.names(t_draws) <- NULL
      } else {
        t_draws <- rmvnorm_prec(mu = new_tmb_obj$env$last.par.best,
          tmb_sd = sd_report, n_sims = sims)
      }
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(mcmc_samples)) {
      t_draws <- mcmc_samples
      if (nsim > 0) {
        if (nsim > ncol(t_draws)) {
          cli_abort("`nsim` must be <= number of MCMC samples.")
        } else {
          t_draws <- t_draws[,seq(ncol(t_draws) - nsim + 1, ncol(t_draws)), drop = FALSE]
        }
      }
      r <- apply(t_draws, 2L, new_tmb_obj$report)
    }
    if (!is.null(mcmc_samples) || sims > 0) {
      if (return_tmb_report) return(r)
      .var <-  switch(sims_var,
        "est" = "proj_eta",
        "est_rf" = "proj_rf",
        "omega_s" = "proj_omega_s_A",
        "zeta_s" = "proj_zeta_s_A",
        "epsilon_st" = "proj_epsilon_st_A_vec",
        sims_var)
      out <- lapply(r, `[[`, .var)

      predtype <- as.integer(model[[1]])
      if (isTRUE(object$family$delta) && sims_var == "est") {
        if (predtype %in% c(1L, NA)) {
          out1 <- lapply(out, function(x) x[, 1L, drop = TRUE])
          out1 <- do.call("cbind", out1)
        }
        if (predtype %in% c(2L, NA)) {
          out2 <- lapply(out, function(x) x[, 2L, drop = TRUE])
          out2 <- do.call("cbind", out2)
        }
        if (is.na(predtype)) {
          out <- object$family[[1]]$linkinv(out1) *
            object$family[[2]]$linkinv(out2)
          if (type != "response") out <- object$family[[2]]$linkfun(out) # transform combined back into link space
        } else if (predtype == 1L) {
          out <- out1
          if (type == "response") out <- object$family[[1]]$linkinv(out)
        } else if (predtype == 2L) {
          out <- out2
          if (type == "response") out <- object$family[[2]]$linkinv(out)
        } else {
          cli_abort("`model` type not valid.")
        }
      } else { # not a delta model OR not sims_var = "est":

        if (isTRUE(object$family$delta) && sims_var != "est" && is.na(model[[1]])) {
          cli_warn("`model` argument was left as NA; defaulting to 1st model component.")
          model <- 1L
        } else {
          model <- as.integer(model)
        }
        if (!isTRUE(object$family$delta)) {
          model <- 1L
        }
        if (length(dim(out[[1]])) == 2L) {
          out <- lapply(out, function(.x) .x[,model])
          out <- do.call("cbind", out)
        } else if (length(dim(out[[1]])) == 3L) {
          xx <- list()
          for (i in seq_len(dim(out[[1]])[2])) {
            xx[[i]] <- lapply(out, function(.x) .x[,i,model])
            xx[[i]] <- do.call("cbind", xx[[i]])
          }
          out <- xx
          if (sims_var == "zeta_s") names(out) <- object$spatial_varying
          if (length(out) == 1L) out <- out[[1]]
        } else {
          cli_abort("Too many dimensions returned from model. Try `return_tmb_report = TRUE` and parse the output yourself.")
        }

        if (type == "response") out <- object$family$linkinv(out)
      }

      if (sims_var == "est") {
        rownames(out) <- nd[[object$time]] # for use in index calcs
        attr(out, "time") <- object$time
        if (type == "response"){
          attr(out, "link") <- "response"
        } else {
          if(isTRUE(object$family$delta)){
            if (is.na(predtype)) {
              attr(out, "link") <- object$family[[2]]$link
            } else if (predtype == 1L) {
              attr(out, "link") <- object$family[[1]]$link
            } else if (predtype == 2L) {
              attr(out, "link") <- object$family[[2]]$link
            } else {
              cli_abort("`model` type not valid.")
            }
          } else {
            attr(out, "link") <- object$family$link
          }
        }
      }

      return(out)
    }

    lp <- new_tmb_obj$env$last.par.best
    r <- new_tmb_obj$report(lp)
    if (return_tmb_report) return(r)

    if (isFALSE(pop_pred)) {
      if (isTRUE(object$family$delta)) {
        nd$est1 <- r$proj_eta[,1]
        nd$est2 <- r$proj_eta[,2]
        nd$est_non_rf1 <- r$proj_fe[,1]
        nd$est_non_rf2 <- r$proj_fe[,2]
        nd$est_rf1 <- r$proj_rf[,1]
        nd$est_rf2 <- r$proj_rf[,2]
        nd$omega_s1 <- r$proj_omega_s_A[,1]
        nd$omega_s2 <- r$proj_omega_s_A[,2]
        for (z in seq_len(dim(r$proj_zeta_s_A)[2])) { # SVC:
          nd[[paste0("zeta_s_", object$spatial_varying[z], "1")]] <- r$proj_zeta_s_A[,z,1]
          nd[[paste0("zeta_s_", object$spatial_varying[z], "2")]] <- r$proj_zeta_s_A[,z,2]
        }
        nd$epsilon_st1 <- r$proj_epsilon_st_A_vec[,1]
        nd$epsilon_st2 <- r$proj_epsilon_st_A_vec[,2]
        if (type == "response" && !se_fit) {
          nd$est1 <- object$family[[1]]$linkinv(nd$est1)
          nd$est2 <- object$family[[2]]$linkinv(nd$est2)
          if (object$tmb_data$poisson_link_delta) {
            .n <- nd$est1 # expected group density (already exp())
            .p <- 1 - exp(-.n) # expected encounter rate
            .w <- nd$est2 # expected biomass per group (already exp())
            .r <- (.n * .w) / .p # (n * w)/p # positive expectation
            nd$est1 <- .p # expected encounter rate
            nd$est2 <- .r # positive expectation
            nd$est <- .n * .w # expected combined value
          } else {
            nd$est <- nd$est1 * nd$est2
          }
        }
      } else {
        nd$est <- r$proj_eta[,1]
        nd$est_non_rf <- r$proj_fe[,1]
        nd$est_rf <- r$proj_rf[,1]
        nd$omega_s <- r$proj_omega_s_A[,1]
        for (z in seq_len(dim(r$proj_zeta_s_A)[2])) { # SVC:
          nd[[paste0("zeta_s_", object$spatial_varying[z])]] <- r$proj_zeta_s_A[,z,1]
        }
        nd$epsilon_st <- r$proj_epsilon_st_A_vec[,1]
        if (type == "response") {
          nd$est <- object$family$linkinv(nd$est)
        }
      }
    }

    nd$sdm_spatial_id <- NULL
    nd$sdm_orig_id <- NULL

    obj <- new_tmb_obj

    if ("visreg_model" %in% names(object)) {
      model <- object$visreg_model
    } else {
      if (visreg_df)
        model <- 1L
    }

    if (se_fit) {
      sr <- TMB::sdreport(new_tmb_obj, bias.correct = FALSE)
      sr_est_rep <- as.list(sr, "Estimate", report = TRUE)
      sr_se_rep <- as.list(sr, "Std. Error", report = TRUE)
      if (pop_pred) {
        proj_eta <- sr_est_rep[["proj_fe"]]
        se <- sr_se_rep[["proj_fe"]]
      } else {
        proj_eta <- sr_est_rep[["proj_eta"]]
        se <- sr_se_rep[["proj_eta"]]
      }

      # For delta models with model = NA, use combined predictions
      if (isTRUE(object$family$delta) && is.na(model) && !pop_pred &&
          "proj_eta_delta" %in% names(sr_est_rep)) {
        nd$est <- sr_est_rep[["proj_eta_delta"]]
        nd$est_se <- sr_se_rep[["proj_eta_delta"]]
      } else {
        if (is.na(model)) model_temp <- 1L else model_temp <- model
        proj_eta <- proj_eta[,model_temp,drop=TRUE]
        se <- se[,model_temp,drop=TRUE]
        nd$est <- proj_eta
        nd$est_se <- se
      }
    }
    if (type == "response" && se_fit) {
      est_name <- if (isTRUE(object$family$delta)) "'est1' and 'est2'" else "'est'"
      msg <- paste0("predict(..., type = 'response', se_fit = TRUE) detected; ",
        "returning the prediction ", est_name, " in link space because the standard errors ",
        "are calculated in link space.")
      cli_warn(msg)
      type <- "link"
    }

    if (pop_pred) {
      if (isTRUE(object$family$delta)) {
        if (type == "response") {
          nd$est1 <- object$family[[1]]$linkinv(r$proj_fe[,1])
          nd$est2 <- object$family[[2]]$linkinv(r$proj_fe[,2])
          if (object$tmb_data$poisson_link_delta) {
            .n <- nd$est1 # expected group density (already exp())
            .p <- 1 - exp(-.n) # expected encounter rate
            .w <- nd$est2 # expected biomass per group (already exp())
            .r <- (.n * .w) / .p # (n * w)/p # positive expectation
            nd$est1 <- .p # expected encounter rate
            nd$est2 <- .r # positive expectation
            nd$est <- .n * .w # expected combined value
          } else {
            nd$est <- nd$est1 * nd$est2
          }
        } else {
          nd$est1 <- r$proj_fe[,1]
          nd$est2 <- r$proj_fe[,2]
          if (is.na(model)) {
            p1 <- object$family[[1]]$linkinv(r$proj_fe[,1])
            p2 <- object$family[[2]]$linkinv(r$proj_fe[,2])
            if (object$tmb_data$poisson_link_delta) {
              nd$est <- nd$est1 + nd$est2
            } else {
              nd$est <- object$family[[2]]$linkfun(p1 * p2)
            }
            if (se_fit) {
              nd$est <- sr_est_rep$proj_rf_delta
              nd$est_se <- sr_se_rep$proj_rf_delta
            }
          }
        }
      } else {
        if (type == "response") {
          nd$est <- object$family$linkinv(r$proj_fe[,1])
        } else {
          nd$est <- r$proj_fe[,1]
        }
      }
    }

    if (pop_pred && visreg_df) {
      nd$est <- r$proj_fe[,model,drop=TRUE] # FIXME re_form_iid??
    }

    orig_dat <- object$tmb_data$y_i
    if (model == 2L && nrow(nd) == nrow(orig_dat) && visreg_df) {
      nd <- nd[!is.na(orig_dat[,2]),,drop=FALSE] # drop NAs from delta positive component
    }

    if ("sdmTMB_fake_year" %in% names(nd)) {
      nd <- nd[!nd$sdmTMB_fake_year,,drop=FALSE]
      nd$sdmTMB_fake_year <- NULL
    }
    if (fake_spatial_added) {
      for (i in 1:2) nd[[xy_cols[[i]]]] <- NULL
    }
    if (sdmTMB_fake_response) {
      nd[[response]] <- NULL
    }

  } else { # We are not dealing with new data:
    if (se_fit) {
      cli_warn(paste0("Standard errors have not been implemented yet unless you ",
        "supply `newdata`. In the meantime you could supply your original data frame ",
        "to the `newdata` argument."))
    }
    nd <- object$data
    lp <- object$tmb_obj$env$last.par.best
    # object$tmb_obj$fn(lp) # call once to update internal structures?
    r <- object$tmb_obj$report(lp)

    nd$est <- r$eta_i[,1] # DELTA FIXME
    # The following is not an error,
    # IID and RW effects are baked into fixed effects for `newdata` in above code:
    nd$est_non_rf <- r$eta_fixed_i[,1] + r$eta_rw_i[,1] + r$eta_iid_re_i[,1] # DELTA FIXME
    nd$est_rf <- r$omega_s_A[,1] + r$epsilon_st_A_vec[,1] # DELTA FIXME
    nd$omega_s <- r$omega_s_A[,1]# DELTA FIXME
    for (z in seq_len(dim(r$zeta_s_A)[2])) { # SVC: # DELTA FIXME
      nd[[paste0("zeta_s_", object$spatial_varying[z])]] <- r$zeta_s_A[,z,1]
    }
    nd$epsilon_st <- r$epsilon_st_A_vec[,1]# DELTA FIXME
    obj <- object
  }

  # clean up:
  if (!object$tmb_data$include_spatial[1]) {
    nd$omega_s1 <- NULL
    nd$omega_s <- NULL
  }
  if (isTRUE(object$family$delta)) {
    if (!object$tmb_data$include_spatial[2]) {
      nd$omega_s2 <- NULL
    }
  }
  if (as.logical(object$tmb_data$spatial_only)[1]) {
    nd$epsilon_st1 <- NULL
    nd$epsilon_st <- NULL
  }
  if (isTRUE(object$family$delta)) {
    if (as.logical(object$tmb_data$spatial_only)[2]) {
      nd$epsilon_st2 <- NULL
    }
  }
  if (!object$tmb_data$spatial_covariate) {
    nd$zeta_s1 <- NULL
    nd$zeta_s1 <- NULL
    nd$zeta_s <- NULL
  }

  if (no_spatial) nd[,xy_cols] <- NULL
  nd[["_sdmTMB_time"]] <- NULL
  if (no_spatial) nd[["est_rf"]] <- NULL
  if (no_spatial) nd[["est_non_rf"]] <- NULL
  row.names(nd) <- NULL

  if (return_tmb_object) {
    return(list(data = nd, report = r, obj = obj, fit_obj = object, pred_tmb_data = tmb_data))
  } else {
    if (visreg_df) {
      # for visreg & related, return consistent objects with lm(), gam() etc.
      if (isTRUE(se_fit)) {
        return(list(fit = nd$est, se.fit = nd$est_se))
      } else {
        return(nd$est)
      }
    } else {
      return(nd) # data frame by default
    }
  }
}

# https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

remove_9000 <- function(x) {
  as.package_version(paste0(
    strsplit(as.character(x), ".", fixed = TRUE)[[1]][1:3],
    collapse = "."
  ))
}

check_sdmTMB_version <- function(version) {
  if (remove_9000(utils::packageVersion("sdmTMB")) >
    remove_9000(version)) {
    msg <- paste0(
      "The installed version of sdmTMB is newer than the version ",
      "that was used to fit this model. It is possible new parameters ",
      "have been added to the TMB model since you fit this model and ",
      "that prediction will fail. We recommend you fit and predict ",
      "from an sdmTMB model with the same version."
      )
    cli_warn(msg)
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

check_visreg <- function(sys_calls) {
  visreg_df <- FALSE
  re_form <- NULL
  se_fit <- FALSE
  if (any(grepl("setupV", substr(sys_calls, 1, 7)))) {
    visreg_df <- TRUE
    re_form <- NA
    if (any(sys_calls == "residuals(fit)")) visreg_df <- FALSE
    # turn on standard error if in a function call
    indx <- which(substr(sys_calls, 1, 10) == "visregPred")
    if (length(indx) > 0 && any(unlist(strsplit(sys_calls[indx], ",")) == " se.fit = TRUE"))
      se_fit <- TRUE
  }
  named_list(visreg_df, se_fit, re_form)
}
