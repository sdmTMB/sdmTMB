#' Extract a relative biomass/abundance index, center of gravity, effective
#' area occupied, or weighted average
#'
#' @param obj A model fitted with [sdmTMB()]. For backwards compatibility,
#'   output from [predict.sdmTMB()] with `return_tmb_object = TRUE` is also
#'   accepted.
#' @param newdata New data (e.g., a prediction grid by year) to pass to
#'   [predict.sdmTMB()]. Not used when `obj` is legacy prediction output.
#' @param bias_correct Should bias correction be implemented via
#'   [TMB::sdreport()]? Bias correction accounts for the non-linear
#'   transformation of random effects when calculating the index. Recommended to
#'   be `TRUE` for final analyses, but can be set to `FALSE` for faster
#'   calculation while experimenting with models. See Thorson and Kristensen
#'   (2016) in the References.
#' @param level The confidence level.
#' @param area Grid cell area for area weighting the index. Can be: (1) a
#'   numeric vector of length `nrow(newdata)` with area for each grid cell, (2)
#'   a single numeric value to apply to all grid cells, or (3) a character value
#'   giving the column name in `newdata` containing areas. See Details for
#'   non-spatial uses of `area` as an integration multiplier.
#' @param offset An optional numeric offset vector with one value per row of
#'   `newdata`.
#' @param silent Logical. Suppress progress messages?
#' @param derived_link Optional override for the inverse link used when
#'   calculating derived quantities such as the index. By default, the fitted
#'   family link is used. Currently supported for non-delta `binomial()` and
#'   `betabinomial()` models fit with `link = "cloglog"`.
#' @param predict_args A named list of less commonly used arguments to pass to
#'   [predict.sdmTMB()]. `newdata` and `offset` should be supplied directly.
#' @param ... Passed to [TMB::sdreport()].
#'
#' @details
#' More generally, `area` is the multiplier used to integrate response-scale
#' predictions. For binomial or beta-binomial models fit to proportions with
#' `weights` specifying the number of trials, predictions are expected
#' proportions per trial. In that case, `area` can be used as a standard number
#' of trials (e.g., hooks per longline set) to obtain an expected-count index.
#' The original fitting `weights` are not automatically reused for index
#' standardization; supply the desired standardization multiplier through
#' `area`.
#'
#' @seealso [get_index_sims()]
#' @return
#' For `get_index()`:
#' A data frame with columns for time, estimate (area-weighted total abundance
#' or biomass), lower and upper confidence intervals, log estimate, and standard
#' error of the log estimate.
#'
#' For `get_cog()`:
#' A data frame with columns for time, estimate (center of gravity: the
#' abundance-weighted mean x and y coordinates), lower and upper confidence
#' intervals, and standard error of center of gravity coordinates.
#'
#' For `get_eao()`:
#' A data frame with columns for time, estimate (effective area occupied: the
#' area required if the population was spread evenly at the arithmetic mean
#' density), lower and upper confidence intervals, log EAO, and standard error
#' of the log EAO estimates.
#'
#' For `get_weighted_average()`:
#' A data frame with columns for time, estimate (weighted average of the
#' provided vector, weighted by predicted density), lower and upper confidence
#' intervals, and standard error of the estimates.
#'
#' @references
#'
#' Geostatistical model-based indices of abundance
#' (along with many newer papers):
#'
#' Shelton, A.O., Thorson, J.T., Ward, E.J., and Feist, B.E. 2014. Spatial
#' semiparametric models improve estimates of species abundance and
#' distribution. Canadian Journal of Fisheries and Aquatic Sciences 71(11):
#' 1655--1666. \doi{10.1139/cjfas-2013-0508}
#'
#' Thorson, J.T., Shelton, A.O., Ward, E.J., and Skaug, H.J. 2015.
#' Geostatistical delta-generalized linear mixed models improve precision for
#' estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci.
#' 72(5): 1297–1310. \doi{10.1093/icesjms/fsu243}
#'
#' Geostatistical model-based center of gravity:
#'
#' Thorson, J.T., Pinsky, M.L., and Ward, E.J. 2016. Model-based inference for
#' estimating shifts in species distribution, area occupied and centre of
#' gravity. Methods Ecol Evol 7(8): 990–1002. \doi{10.1111/2041-210X.12567}
#'
#' Geostatistical model-based effective area occupied:
#'
#' Thorson, J.T., Rindorf, A., Gao, J., Hanselman, D.H., and Winker, H. 2016.
#' Density-dependent changes in effective area occupied for
#' sea-bottom-associated marine fishes. Proceedings of the Royal Society B:
#' Biological Sciences 283(1840): 20161853.
#' \doi{10.1098/rspb.2016.1853}
#'
#' Bias correction:
#'
#' Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method for
#' bias correction in statistical models using random effects, with spatial and
#' population dynamics examples. Fisheries Research 175: 66–74.
#' \doi{10.1016/j.fishres.2015.11.016}
#'
#' @examples
#' \donttest{
#' if (ggplot2_installed()) {
#' library(ggplot2)
#'
#' # use a small number of knots for this example to make it fast:
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 60)
#'
#' # fit a spatiotemporal model:
#' m <- sdmTMB(
#'  data = pcod,
#'  formula = density ~ 0 + as.factor(year),
#'  time = "year", mesh = mesh, family = tweedie(link = "log")
#' )
#'
#' # prepare a prediction grid:
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#'
#' # biomass index:
#' ind <- get_index(m, newdata = nd, bias_correct = TRUE)
#' ind
#' ggplot(ind, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
#'   ylim(0, NA)
#'
#' # do that in 2 chunks
#' # only necessary for very large grids to save memory
#' # will be slower but save memory
#' # note the first argument is the model fit object:
#' ind <- get_index_split(m, newdata = nd, nsplit = 2, bias_correct = TRUE)
#'
#' # center of gravity:
#' cog <- get_cog(m, newdata = nd, format = "wide")
#' cog
#' ggplot(cog, aes(est_x, est_y, colour = year)) +
#'   geom_point() +
#'   geom_linerange(aes(xmin = lwr_x, xmax = upr_x)) +
#'   geom_linerange(aes(ymin = lwr_y, ymax = upr_y)) +
#'   scale_colour_viridis_c()
#'
#' # effective area occupied:
#' eao <- get_eao(m, newdata = nd)
#' eao
#' ggplot(eao, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
#'   ylim(0, NA)
#'
#' # weighted average (e.g., depth-weighted by biomass):
#' wa <- get_weighted_average(m, newdata = nd, vector = nd$depth)
#' wa
#' ggplot(wa, aes(year, est)) + geom_line() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
#' }
#' }
#' @export
get_index <- function(obj, newdata = NULL, bias_correct = TRUE, level = 0.95,
  area = 1, offset = NULL, silent = TRUE, derived_link = NULL,
  predict_args = list(), ...)  {
  if (is.logical(newdata)) {
    if (length(newdata) != 1L || is.na(newdata) || !missing(bias_correct)) {
      cli_abort("A logical `newdata` is ambiguous; name `newdata` and `bias_correct` explicitly.")
    }
    cli_warn("Interpreting the second positional argument as `bias_correct`; use `bias_correct =` explicitly.")
    bias_correct <- newdata
    newdata <- NULL
  }
  area_missing <- missing(area)
  obj <- .prepare_index_input(obj, newdata, offset, predict_args)
  d <- get_generic(obj, value_name = "link_total",
    bias_correct = bias_correct, level = level, trans = exp, area = area,
    derived_link = derived_link, area_missing = area_missing, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d$type <- "index"
  d
}

.validate_derived_link <- function(fit_obj, derived_link) {
  if (is.null(derived_link)) {
    return(NULL)
  }
  if (!is.character(derived_link) || length(derived_link) != 1L || is.na(derived_link)) {
    cli_abort("`derived_link` must be a single character string naming a link.")
  }
  if (!derived_link %in% names(.valid_link)) {
    cli_abort(c(
      "`derived_link` is not valid.",
      "i" = "Choose one of {.val {names(.valid_link)}}."
    ))
  }
  family_spec <- .object_family_spec(fit_obj, caller = "`get_index()`")
  if (.family_spec_is_multi_family(family_spec)) {
    cli_abort("`derived_link` is not currently supported for multi-family models.")
  }
  family <- family_spec$family
  if (.family_spec_has_two_components(family_spec)) {
    cli_abort("`derived_link` is not currently supported for delta or hurdle families.")
  }
  if (!family$family %in% c("binomial", "betabinomial")) {
    cli_abort("`derived_link` is currently only supported for binomial and betabinomial models.")
  }
  if (!identical(family$link, "cloglog")) {
    cli_abort("`derived_link` is currently only supported when the fitted family uses `link = 'cloglog'`.")
  }
  derived_link
}

.resolve_link_pred <- function(fit_obj, tmb_data, derived_link = NULL) {
  if (is.null(tmb_data$link_pred)) {
    if (!is.null(tmb_data$link)) {
      tmb_data$link_pred <- tmb_data$link
    } else {
      family <- .object_family_spec(fit_obj, caller = "`get_index()`")$family
      tmb_data$link_pred <- unname(.valid_link[family$link])
    }
  }
  if (!is.null(derived_link)) {
    derived_link <- .validate_derived_link(fit_obj, derived_link)
    tmb_data$link_pred[] <- unname(.valid_link[derived_link])
  }
  tmb_data$link_pred
}

chunk_time <- function(x, chunks) {
  assert_that(is.numeric(chunks))
  assert_that(chunks > 0)
  chunks <- as.integer(chunks)
  if (chunks == 1L) {
    list(x)
  } else {
    return(split(x, cut(seq_along(x), chunks, labels = FALSE)))
  }
}

.prepare_index_input <- function(obj, newdata = NULL, offset = NULL,
  predict_args = list()) {
  if (!is.list(predict_args)) {
    cli_abort("`predict_args` must be a list.")
  }
  if (length(predict_args)) {
    if (is.null(names(predict_args)) || any(!nzchar(names(predict_args)))) {
      cli_abort("All elements of `predict_args` must be named.")
    }
    reserved <- c("object", "newdata", "offset", "return_tmb_object",
      "return_tmb_report", "return_tmb_data", "nsim", "mcmc_samples")
    invalid <- intersect(names(predict_args), reserved)
    if (length(invalid)) {
      cli_abort(c(
        "Reserved arguments cannot be supplied in `predict_args`.",
        "i" = "Supply {.arg {invalid}} directly or remove them."
      ))
    }
  }

  if (inherits(obj, "sdmTMB")) {
    if (is.null(newdata)) {
      if (!is.null(offset)) {
        cli_abort("`offset` can only be supplied with `newdata`.")
      }
      if (length(predict_args)) {
        cli_abort("`predict_args` can only be supplied with `newdata`.")
      }
      return(obj)
    }
    if (!is.data.frame(newdata)) {
      cli_abort("`newdata` must be a data frame.")
    }
    if (!is.null(offset) &&
        (!is.numeric(offset) || length(offset) != nrow(newdata))) {
      cli_abort("`offset` must be a numeric vector with one value per row of `newdata`.")
    }
    args <- c(list(
      object = obj,
      newdata = newdata,
      offset = offset,
      return_tmb_data = TRUE
    ), predict_args)
    tmb_data <- do.call(predict.sdmTMB, args)
    return(list(data = newdata, fit_obj = obj, pred_tmb_data = tmb_data))
  }

  if (!is.null(newdata)) {
    cli_abort("`newdata` cannot be supplied when `obj` is prediction output.")
  }
  if (!is.null(offset)) {
    cli_abort(c(
      "`offset` cannot be added to existing prediction output.",
      "i" = "Use `get_index(fit, newdata = ..., offset = ...)` instead."
    ))
  }
  if (length(predict_args)) {
    cli_abort("`predict_args` cannot be supplied when `obj` is prediction output.")
  }
  if (!is.list(obj) || is.null(obj$fit_obj) ||
      is.null(obj$pred_tmb_data) || is.null(obj$data)) {
    cli_abort(c(
      "`obj` must be an sdmTMB fit or legacy prediction output.",
      "i" = "The preferred form is `get_index(fit, newdata = data)`.",
      "i" = "Legacy prediction output requires `return_tmb_object = TRUE`."
    ))
  }
  obj
}

#' @rdname get_index
#' @param nsplit The number of splits to do the calculation in. For memory
#'   intensive operations (large grids and/or models), it can be helpful to
#'   do the prediction, area integration, and bias correction on subsets of
#'   time slices (e.g., years) instead of all at once. If `nsplit > 1`, this
#'   will usually be slower but with reduced memory use.
#' @export
get_index_split <- function(
    obj, newdata, bias_correct = FALSE, nsplit = 1,
    level = 0.95, area = 1, offset = NULL, silent = FALSE, predict_args = list(),
    derived_link = NULL, ...) {
  if (!inherits(obj, "sdmTMB")) {
    cli_abort("get_index_split() is meant to be run on a fitted object from sdmTMB() and not a prediction object as in get_index().")
  }
  if (!is.list(predict_args)) {
    cli_abort("`predict_args` must be a list.")
  }

  times <- sort(obj$time_lu$time_from_data)
  time_chunks <- chunk_time(times, nsplit)

  if ("offset" %in% names(predict_args)) {
    if (!is.null(offset)) {
      cli_abort("Supply `offset` directly or in `predict_args`, not both.")
    }
    cli_warn("Supplying `offset` in `predict_args` is deprecated; use the top-level `offset` argument.")
    offset <- predict_args$offset
    predict_args$offset <- NULL
  }
  if (is.null(offset)) {
    offset <- rep(0, nrow(newdata))
  }
  if (!is.numeric(offset) || length(offset) != nrow(newdata)) {
    cli_abort("`offset` must be a numeric vector with one value per row of `newdata`.")
  }
  if (is.character(area)) {
    if (length(area) != 1L || !area %in% names(newdata)) {
      cli_abort("A character `area` must name a column in `newdata`.")
    }
    area <- newdata[[area]]
  }
  if (length(area) == 1L) {
    area <- rep(area, nrow(newdata))
  }
  if (length(area) != nrow(newdata)) {
    cli_abort("`area` must have length one or one value per row of `newdata`.")
  }

  msg <- paste0("Calculating index in ", nsplit, " chunks")
  if (!silent) cli::cli_progress_bar(msg, total = length(time_chunks))
  index_list <- list()
  for (i in seq_along(time_chunks)) {
    if (!silent) {
      cli::cli_progress_update(set = i, total = length(time_chunks), force = TRUE)
    }
    this_chunk_i <- newdata[[obj$time]] %in% time_chunks[[i]]
    nd <- newdata[this_chunk_i, , drop = FALSE]

    index_list[[i]] <-
      get_index(
        obj,
        newdata = nd,
        bias_correct = bias_correct,
        level = level,
        area = area[this_chunk_i],
        offset = offset[this_chunk_i],
        silent = TRUE,
        derived_link = derived_link,
        predict_args = predict_args,
        ...
      )
  }
  if (!silent) cli::cli_progress_done()
  do.call(rbind, index_list)
}

#' @rdname get_index
#' @param format Long or wide.
#' @export
get_cog <- function(obj, newdata = NULL, bias_correct = FALSE, level = 0.95,
  format = c("long", "wide"), area = 1, offset = NULL, silent = TRUE,
  derived_link = NULL, predict_args = list(), ...)  {

  if (is.logical(newdata)) {
    if (length(newdata) != 1L || is.na(newdata) || !missing(bias_correct)) {
      cli_abort("A logical `newdata` is ambiguous; name `newdata` and `bias_correct` explicitly.")
    }
    cli_warn("Interpreting the second positional argument as `bias_correct`; use `bias_correct =` explicitly.")
    bias_correct <- newdata
    newdata <- NULL
  }
  area_missing <- missing(area)
  obj <- .prepare_index_input(obj, newdata, offset, predict_args)

  is_fit_obj <- inherits(obj, "sdmTMB")
  fit_obj <- if (is_fit_obj) obj else obj$fit_obj
  pred_tmb_data <- if (is_fit_obj) obj$tmb_data else obj$pred_tmb_data
  xy_cols <- fit_obj$spde$xy_cols
  if (is.null(xy_cols)) {
    cli_abort("`get_cog()` requires x/y coordinates and isn't available for areal models.")
  }
  # for a bare `do_index = TRUE` fit, `obj$data` is the *observation* data, so
  # coordinates must come from the stored prediction TMB data instead
  if (!is_fit_obj && all(xy_cols %in% names(obj$data))) {
    x_vec <- obj$data[[xy_cols[[1]]]]
    y_vec <- obj$data[[xy_cols[[2]]]]
  } else if (!is.null(pred_tmb_data$proj_lon) &&
             !is.null(pred_tmb_data$proj_lat)) {
    x_vec <- pred_tmb_data$proj_lon
    y_vec <- pred_tmb_data$proj_lat
  } else {
    cli_abort("Prediction data must include the x/y columns used for the model.")
  }
  # x and y share one objective function, sdreport, and bias correction
  d_xy <- get_generic(obj, value_name = "weighted_avg",
    bias_correct = bias_correct, level = level, trans = I, area = area,
    vector = cbind(x_vec, y_vec), derived_link = derived_link,
    area_missing = area_missing, ...)
  d_x <- d_xy[[1]]
  d_y <- d_xy[[2]]
  d_x <- d_x[, names(d_x) != "trans_est", drop = FALSE]
  d_y <- d_y[, names(d_y) != "trans_est", drop = FALSE]
  d_x$coord <- "X"
  d_y$coord <- "Y"
  d <- rbind(d_x, d_y)
  format <- match.arg(format)
  if (format == "wide") {
    x <- d[d$coord == "X", c("est", "lwr", "upr", "se"),drop=FALSE]
    y <- d[d$coord == "Y", c("est", "lwr", "upr", "se"),drop=FALSE]
    names(x) <- paste0(names(x), "_", "x")
    names(y) <- paste0(names(y), "_", "y")
    d <- cbind(d[d$coord == "X", fit_obj$time, drop=FALSE], cbind(x, y))
  }
  d$type <- "cog"
  d
}

#' @rdname get_index
#' @param vector A numeric vector of the same length as the prediction data,
#'   containing the values to be averaged (e.g., depth, temperature).
#' @export
get_weighted_average <- function(obj, newdata = NULL, vector, bias_correct = FALSE,
  level = 0.95, area = 1, offset = NULL, silent = TRUE, derived_link = NULL,
  predict_args = list(), ...)  {

  if (!is.data.frame(newdata) && missing(vector) && !is.null(newdata)) {
    cli_warn("Interpreting the second positional argument as `vector`; use `vector =` explicitly.")
    vector <- newdata
    newdata <- NULL
  }
  area_missing <- missing(area)
  obj <- .prepare_index_input(obj, newdata, offset, predict_args)

  d <- get_generic(obj, value_name = "weighted_avg",
    bias_correct = bias_correct, level = level, trans = I, area = area,
    vector = vector, derived_link = derived_link, area_missing = area_missing, ...)
  d <- d[, names(d) != "trans_est", drop = FALSE]
  d$type <- "weighted_average"
  d
}

#' @rdname get_index
#' @export
get_eao <- function(obj,
  newdata = NULL,
  bias_correct = FALSE,
  level = 0.95,
  area = 1,
  offset = NULL,
  silent = TRUE,
  derived_link = NULL,
  predict_args = list(),
  ...
)  {

  if (is.logical(newdata)) {
    if (length(newdata) != 1L || is.na(newdata) || !missing(bias_correct)) {
      cli_abort("A logical `newdata` is ambiguous; name `newdata` and `bias_correct` explicitly.")
    }
    cli_warn("Interpreting the second positional argument as `bias_correct`; use `bias_correct =` explicitly.")
    bias_correct <- newdata
    newdata <- NULL
  }
  area_missing <- missing(area)
  obj <- .prepare_index_input(obj, newdata, offset, predict_args)

  d <- get_generic(obj, value_name = c("log_eao"),
    bias_correct = bias_correct, level = level, trans = exp, area = area,
    derived_link = derived_link, area_missing = area_missing, ...)
  names(d)[names(d) == "trans_est"] <- "log_est"
  d$type <- "eoa"
  d
}

# Estimates and SEs of an index objective's reports, as from
# `summary(<sdreport>, "report")`. `args` are make_sdmTMB_adfun() arguments.
index_report <- function(fit, args, par, ...) {
  if (...length() == 0L) {
    out <- joint_precision_report(fit, args, par)
    if (!is.null(out)) return(out)
  }
  new_obj <- do.call(make_sdmTMB_adfun, args)
  summary(index_sdreport(fit, new_obj, par, ...), "report")
}

# sdreport()'s delta-method covariance of the reports is J Q^-1 J', with J
# their Jacobian in all (fixed and random) parameters at the fitted mode and Q
# the fit's joint precision. Projection rows don't enter the likelihood, so Q
# applies; this skips re-optimizing the random effects and computing their
# marginal variances. NULL unless it verifiably applies.
joint_precision_report <- function(fit, args, par) {
  sr <- fit$sd_report
  Q <- sr$jointPrecision
  random <- fit$tmb_obj$env$random
  if (!is.null(args$profile) || is.null(Q) ||
      !length(random) || !isTRUE(sr$pdHess) ||
      !identical(unname(sr$par.fixed), unname(par))) {
    return(NULL)
  }
  args$random <- NULL
  obj <- do.call(make_sdmTMB_adfun, c(args, ADreport = TRUE))
  x <- obj$par
  if (!identical(names(x), colnames(Q)) ||
      length(x) != length(par) + length(sr$par.random)) {
    return(NULL)
  }
  x[random] <- sr$par.random
  x[-random] <- par
  J <- obj$gr(x)
  QiJt <- tryCatch(as.matrix(Matrix::solve(Q, t(J))), error = function(e) NULL)
  if (is.null(QiJt)) return(NULL)
  se <- sqrt(rowSums(J * t(QiJt)))
  if (any(!is.finite(se))) return(NULL)
  est <- obj$fn(x)
  out <- cbind(Estimate = as.numeric(est), `Std. Error` = se)
  rownames(out) <- names(est)
  out
}

# `sdreport()` for an index objective. Projection rows don't enter the
# likelihood, so at the fitted parameters the fixed-effect Hessian is the fit's
# and needn't be recomputed (about two gradient evaluations per fixed effect).
index_sdreport <- function(fit, new_obj, par,
  hessian.fixed = fit_hessian_fixed(fit, new_obj, par), ...) {
  sdreport_sdmTMB(new_obj, par.fixed = par, hessian.fixed = hessian.fixed, ...)
}

# The fit's fixed-effect Hessian, or NULL (recompute) unless it verifiably
# applies to `new_obj` at `par`.
fit_hessian_fixed <- function(fit, new_obj, par) {
  sr <- fit$sd_report
  if (!is.null(fit$control$profile) || is.null(sr) || !isTRUE(sr$pdHess) ||
      !identical(unname(sr$par.fixed), unname(par)) ||
      !isTRUE(all.equal(as.numeric(new_obj$fn(par)), fit$model$objective,
        tolerance = 1e-8))) {
    return(NULL)
  }
  H <- tryCatch(solve(sr$cov.fixed), error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) NULL else H
}

get_generic <- function(obj, value_name, bias_correct = FALSE, level = 0.95,
  trans = I, area = 1, vector = NULL, silent = TRUE, derived_link = NULL,
  area_missing = FALSE, ...) {

  # if offset is a character vector, use the value in the dataframe
  if (is.character(area)) {
    area <- obj$data[[area]]
  }

  reinitialize(obj$fit_obj)

  is_fit_obj <- inherits(obj, "sdmTMB")
  # `do_index = TRUE` fits precompute only the index totals
  # (`calc_index_totals`) with the fit-time `area`; other derived quantities
  # and explicit `area` overrides must rebuild the objective function or the
  # requested values are missing from (or wrong in) the stored sdreport
  use_precomputed <- is_fit_obj &&
    isTRUE(obj$do_index) &&
    value_name[[1]] == "link_total" &&
    is.null(derived_link) &&
    isTRUE(area_missing)
  rebuild_from_fit <- is_fit_obj &&
    value_name[[1]] %in% c("link_total", "weighted_avg", "log_eao") &&
    !use_precomputed
  # Only these reports need standard errors from a rebuilt objective
  adreport <- c(value_name,
    switch(value_name[[1]], link_total = "total", log_eao = "eao"))

  if (!use_precomputed && !rebuild_from_fit) {
    if (is.null(obj$pred_tmb_data$proj_X_ij) ||
        is.null(obj$pred_tmb_data$proj_time_include)) {
      cli_abort(c(
        "Prediction data needed for index calculation are missing.",
        "i" = "Use `get_index(fit, newdata = data)` with the current sdmTMB version."
      ))
    }

    nr1 <- nrow(obj$data)
    nr2 <- nrow(obj$pred_tmb_data$proj_X_ij[[1]])
    if (nr1 != nr2) {
      cli_abort(c("Predicted data appears to be modified after prediction",
        "i" = "Please filter `newdata` before predicting."))
    }

    if (bias_correct && obj$fit_obj$control$parallel > 1) {
      cli_warn("Bias correction can be slower with multiple cores; using 1 core.")
      obj$fit_obj$control$parallel <- 1L
    }

    assert_that(!is.null(area))
    if (length(area) > 1L) {
      n_fakend <- if (!is.null(obj$fake_nd)) nrow(obj$fake_nd) else 0L
      area <- c(area, rep(1, n_fakend)) # pad area with any extra time
    }
    if (length(area) != nrow(obj$pred_tmb_data$proj_X_ij[[1]]) && length(area) != 1L) {
      cli_abort("`area` should be of the same length as `nrow(newdata)` or of length 1.")
    }
    if (length(area) == 1L)
      area <- rep(area, nrow(obj$pred_tmb_data$proj_X_ij[[1]]))

    tmb_data <- obj$pred_tmb_data
    if (is.null(tmb_data$proj_time_include)) {
      cli_abort("Missing `proj_time_include` in prediction data. Please recreate the index input with the current sdmTMB version.")
    }
    tmb_data$link_pred <- .resolve_link_pred(obj$fit_obj, tmb_data, derived_link)
    tmb_data$area_i <- area
    if (value_name[1] == "link_total")
      tmb_data$calc_index_totals <- 1L
    if (value_name[1] == "log_eao")
      tmb_data$calc_eao <- 1L
    if (value_name[1] == "weighted_avg") {
      if (is.null(vector)) {
        cli_abort("A vector must be provided for weighted average calculation.")
      }
      if (NROW(vector) != nrow(obj$pred_tmb_data$proj_X_ij[[1]])) {
        cli_abort("`vector` should be of the same length as `nrow(newdata)`.")
      }
      tmb_data$proj_vector <- vector
      tmb_data$calc_weighted_avg <- 1L
    }

    pars <- get_pars(obj$fit_obj)

    eps_name <- "eps_index" # FIXME break out into function; add for COG?
    pars[[eps_name]] <- numeric(0)

    args <- list(
      data = tmb_data,
      parameters = pars,
      profile = obj$fit_obj$control$profile,
      map = obj$fit_obj$tmb_map,
      random = obj$fit_obj$tmb_random,
      backend = backend_sdmTMB(obj$fit_obj),
      silent = silent,
      adreport = adreport
    )
    # bias correction is done below
    ssr <- index_report(obj$fit_obj, args, obj$fit_obj$model$par, ...)
  } else if (rebuild_from_fit) {
    reinitialize(obj)
    if (bias_correct && obj$control$parallel > 1) {
      cli_warn("Bias correction can be slower with multiple cores; using 1 core.")
      obj$control$parallel <- 1L
    }

    assert_that(!is.null(area))
    if (isTRUE(area_missing)) {
      area <- obj$tmb_data$area_i
    }
    if (length(area) != length(obj$tmb_data$area_i) && length(area) != 1L) {
      cli_abort("`area` should be of the same length as the original index prediction data or of length 1.")
    }
    tmb_data <- obj$tmb_data
    if (is.null(tmb_data$proj_time_include)) {
      cli_abort("Missing `proj_time_include` in fitted data. Please refit or re-run with the current sdmTMB version.")
    }
    if (length(area) == 1L) {
      tmb_data$area_i <- rep(area, length(tmb_data$area_i))
    } else {
      tmb_data$area_i <- area
    }
    tmb_data$link_pred <- .resolve_link_pred(obj, tmb_data, derived_link)
    if (value_name[1] == "link_total")
      tmb_data$calc_index_totals <- 1L
    if (value_name[1] == "log_eao")
      tmb_data$calc_eao <- 1L
    if (value_name[1] == "weighted_avg") {
      if (is.null(vector)) {
        cli_abort("A vector must be provided for weighted average calculation.")
      }
      if (NROW(vector) != nrow(tmb_data$proj_X_ij[[1]])) {
        cli_abort("`vector` should be of the same length as the original index prediction data.")
      }
      tmb_data$proj_vector <- vector
      tmb_data$calc_weighted_avg <- 1L
    }

    pars <- get_pars(obj)
    eps_name <- "eps_index"
    pars[[eps_name]] <- numeric(0)

    args <- list(
      data = tmb_data,
      parameters = pars,
      profile = obj$control$profile,
      map = obj$tmb_map,
      random = obj$tmb_random,
      backend = backend_sdmTMB(obj),
      silent = silent,
      adreport = adreport
    )
    ssr <- index_report(obj, args, obj$model$par, ...)
    obj <- list(fit_obj = obj)
  } else {
    # already done in sdmTMB(do_index = TRUE)
    ssr <- summary(obj$sd_report, "report")
    pars <- get_pars(obj)
    tmb_data <- obj$tmb_data
    if (is.null(tmb_data$proj_time_include)) {
      cli_abort("Missing `proj_time_include` in fitted data. Please refit or re-run with the current sdmTMB version.")
    }
    tmb_data$link_pred <- .resolve_link_pred(obj, tmb_data, derived_link)
    obj <- list(fit_obj = obj) # to match regular format
    eps_name <- "eps_index"
  }
  if (bias_correct && value_name[[1]] %in% c("link_total", "weighted_avg", "log_eao")) {
    # extract and modify parameters
    .n <- sum(row.names(ssr) == switch(value_name[[1]],
      link_total = "total", weighted_avg = "weighted_avg", log_eao = "eao"))
    pars[[eps_name]] <- rep(0, .n)
    new_values <- rep(0, .n)
    names(new_values) <- rep(eps_name, length(new_values))
    fixed <- c(obj$fit_obj$model$par, new_values)
    new_obj2 <- make_sdmTMB_adfun(
      data = tmb_data,
      parameters = pars,
      map = obj$fit_obj$tmb_map,
      profile = obj$fit_obj$control$profile,
      random = obj$fit_obj$tmb_random,
      backend = backend_sdmTMB(obj$fit_obj),
      silent = silent,
      # Ratio quantities (weighted_avg, eao) tag their sums in the template;
      # only the internal inner optimizer uses those tags (lowrank = TRUE)
      # and they are much faster there. The external optimizer is faster for
      # the index total.
      intern = value_name[[1]] != "link_total",
      inner.control = list(sparse = TRUE, lowrank = TRUE, trace = FALSE)
    )
    gradient <- new_obj2$gr(fixed)
    corrected_vals <- gradient[names(fixed) == eps_name]
  } else {
    if (value_name[[1]] == "link_total" || value_name[[1]] == "weighted_avg" || value_name[[1]] == "log_eao")
      cli_inform(c("Bias correction is turned off.", "
        It is recommended to turn this on for final inference."))
  }
  log_total <- ssr[row.names(ssr) %in% value_name, , drop = FALSE]
  row.names(log_total) <- NULL
  d <- as.data.frame(log_total)
  names(d) <- c("trans_est", "se")
  if (bias_correct) {
    if (value_name[[1]] == "weighted_avg") {
      d$trans_est <- corrected_vals
    } else {
      d$trans_est <- log(corrected_vals)
    }
    d$est <- corrected_vals
  } else {
    d$est <- as.numeric(trans(d$trans_est))
  }
  d$lwr <- as.numeric(trans(d$trans_est + stats::qnorm((1-level)/2) * d$se))
  d$upr <- as.numeric(trans(d$trans_est + stats::qnorm(1-(1-level)/2) * d$se))

  # also grab natural space SE:
  if ("link_total" %in% value_name) {
    .total <- ssr[row.names(ssr) %in% "total", , drop = FALSE]
    d$se_natural <- as.numeric(.total[,2])
  }

  # A matrix `vector` (e.g., x and y for COG) stacks one time series per
  # column; return one data frame per column
  if (value_name[[1]] == "weighted_avg" && NCOL(vector) > 1L) {
    chunk <- rep(seq_len(NCOL(vector)), each = nrow(d) / NCOL(vector))
    return(lapply(split(d, chunk), finish_generic, obj = obj,
      tmb_data = tmb_data, value_name = value_name))
  }
  finish_generic(d, obj, tmb_data, value_name)
}

# Match rows of derived quantities to time steps and drop unpredicted ones
finish_generic <- function(d, obj, tmb_data, value_name) {
  time_name <- obj$fit_obj$time
  time_include <- NULL
  if (!is.null(tmb_data) && "proj_time_include" %in% names(tmb_data)) {
    time_include <- as.integer(tmb_data$proj_time_include)
  }
  lu <- obj$fit_obj$time_lu
  if (!is.null(time_include) && length(time_include) == nrow(d)) {
    tt <- lu$time_from_data[time_include != 0]
    d <- d[time_include != 0, ,drop=FALSE]
    keep <- !is.na(d$est)
    d <- d[keep, ,drop=FALSE]
    tt <- tt[keep]
  } else {
    if ("pred_tmb_data" %in% names(obj)) { # standard case
      ii <- sort(unique(obj$pred_tmb_data$proj_year))
    } else { # fit with do_index = TRUE
      ii <- sort(unique(obj$fit_obj$tmb_data$proj_year))
    }
    d <- d[d$est != 0, ,drop=FALSE] # these were not predicted on
    d <- d[!is.na(d$est), ,drop=FALSE] # these were not predicted on
    tt <- lu$time_from_data[match(ii, lu$year_i)]
  }
  if (nrow(d) == 0L) {
    msg <- c(
      "There were no results returned by TMB.",
      "It's possible TMB ran out of memory.",
      "You could try a computer with more RAM or see the function `get_index_split()` with `nsplit > 1`",
      "which lets you split the TMB sdreport and bias correction into chunks."
    )
    cli_abort(msg)
  }
  d[[time_name]] <- tt

  # remove padded extra time fake data:
  if (!is.null(obj$fake_nd)) {
    d <- d[!d[[obj$fit_obj$time]] %in% obj$fake_nd[[obj$fit_obj$time]], ,drop = FALSE]
  }
  if ("do_index_time_missing_from_nd" %in% names(obj$fit_obj)) {
    d <- d[!d[[obj$fit_obj$time]] %in% obj$fit_obj$do_index_time_missing_from_nd, ,drop = FALSE]
  }
  if (!"link_total" %in% value_name) {
    d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se'), drop = FALSE]
  } else {
    d[,c(time_name, 'est', 'lwr', 'upr', 'trans_est', 'se', 'se_natural'), drop = FALSE]
  }
}
