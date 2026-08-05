#' Project from an \pkg{sdmTMB} model using simulation
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' @description Project forward in time from an \pkg{sdmTMB} model using a
#' simulation approach for computational efficiency.
#' This can be helpful for calculating predictive intervals for long
#' projections where including those time elements in `extra_time` during model
#' estimation can be slow.
#'
#' @description Inspiration for this approach comes from the \pkg{VAST} function
#' `project_model()`.
#'
#' @param object A fitted model from [sdmTMB()].
#' @param newdata A new data frame to predict on. It should contain both the
#'   historical time elements and any new time elements to project to. Time
#'   values must be numeric and align with the constant cadence in the fitted
#'   model. Missing intermediate future time elements are simulated.
#' @param nsim Number of projection simulations.
#' @param uncertainty How to sample uncertainty for the fitted parameters:
#'   `"both"` for the joint fixed and random effect precision matrix,
#'   `"random"` for the random effect precision matrix (holding the fixed
#'   effects at their MLE), or `"none"` for neither.
#' @param silent Logical. Suppress progress messages?
#' @param sims_var Optional single element to extract from the \pkg{TMB}
#'   simulation report. If `NULL`, the default, the usual named projection list
#'   is returned. If supplied, the simulation dimension is last. See also
#'   `return_tmb_report`.
#' @param sim_re A vector of `0`s and `1`s representing which random effects to
#'   simulate in the projection. Generally, leave this untouched. Order is:
#'   spatial fields, spatiotemporal fields, spatially varying coefficient
#'   fields, random intercepts, time-varying coefficients, smoothers.
#'   The default is to simulate spatiotemporal fields and time-varying
#'   coefficients, when present.
#' @param return_tmb_report Return the \pkg{TMB} report from `simulate()`? This
#'   lets you parse out whatever elements you want from the simulation including
#'   grabbing multiple elements from one set of simulations. See examples.
#' @param ... Passed to [predict.sdmTMB()].
#'
#' @references `project_model()` in the \pkg{VAST} package.
#' @author
#' J.T. Thorson wrote the original version in the \pkg{VAST} package.
#' S.C. Anderson wrote this version inspired by the \pkg{VAST} version with
#' help from A.J. Allyn.
#' @importFrom cli cli_abort cli_inform cli_warn
#'
#' @return
#' Default: a list with elements `est` and `epsilon_st` (if spatiotemporal
#' effects are present). Each list element includes a matrix with rows
#' corresponding to rows in `newdata` and `nsim` columns. For delta models, the
#' components are `est1`, `est2`, `epsilon_st`, and `epsilon_st2` for the 1st
#' and 2nd linear predictors. In all cases, these returned values are in *link*
#' space.
#'
#' If `sims_var` is supplied, an array containing that report element, with the
#' simulation dimension last. If `return_tmb_report = TRUE`, a list of
#' \pkg{TMB} reports returned by `simulate()`. Run `names()` on an element of the
#' output to see the available `sims_var` options.
#' @export
#'
#' @examples
#' \donttest{
#' if (ggplot2_installed()) {
#' library(ggplot2)
#'
#' mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 25)
#' historical_years <- 2004:2022
#' to_project <- 10
#' future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
#' all_years <- c(historical_years, future_years)
#' proj_grid <- replicate_df(wcvi_grid, "year", all_years)
#'
#' # we could fit our model like this, but for long projections, this becomes slow:
#' if (FALSE) {
#'   fit <- sdmTMB(
#'     catch_weight ~ 1,
#'     time = "year",
#'     offset = log(dogfish$area_swept),
#'     extra_time = all_years, #< note that all years here
#'     spatial = "on",
#'     spatiotemporal = "ar1",
#'     data = dogfish,
#'     mesh = mesh,
#'     family = tweedie(link = "log")
#'   )
#' }
#'
#' # instead, we could fit our model like this and then take simulation draws
#' # from the projection time period:
#' fit2 <- sdmTMB(
#'   catch_weight ~ 1,
#'   time = "year",
#'   offset = log(dogfish$area_swept),
#'   extra_time = historical_years, #< does *not* include projection years
#'   spatial = "on",
#'   spatiotemporal = "ar1",
#'   data = dogfish,
#'   mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#'
#' # we will only use 20 `nsim` so this example runs quickly
#' # you will likely want many more (> 200) in practice so the result
#' # is relatively stable
#'
#' set.seed(1)
#' out <- project(fit2, newdata = proj_grid, nsim = 20)
#' names(out)
#' est_mean <- apply(out$est, 1, mean) # summarize however you'd like
#' est_se <- apply(out$est, 1, sd)
#'
#' # visualize:
#' proj_grid$est_mean <- est_mean
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = est_mean)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   coord_fixed() +
#'   scale_fill_viridis_c() +
#'   ggtitle("Projection simulation (mean)")
#'
#' # visualize the spatiotemporal random fields:
#' proj_grid$eps_mean <- apply(out$epsilon_st, 1, mean)
#' proj_grid$eps_se <- apply(out$epsilon_st, 1, sd)
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = eps_mean)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   scale_fill_gradient2() +
#'   coord_fixed() +
#'   ggtitle("Projection simulation\n(spatiotemporal fields)")
#'
#' ggplot(subset(proj_grid, year > 2021), aes(X, Y, fill = eps_se)) +
#'   geom_raster() +
#'   facet_wrap(~year) +
#'   scale_fill_viridis_c() +
#'   coord_fixed() +
#'   ggtitle("Projection simulation\n(spatiotemporal fields standard error)")
#' }
#' }
project <- function(
    object,
    newdata,
    nsim = 1,
    uncertainty = c("both", "random", "none"),
    silent = FALSE,
    sims_var = NULL,
    sim_re = c(0, 1, 0, 0, 1, 0),
    return_tmb_report = FALSE,
    ...) {
  assert_that(inherits(object, "sdmTMB"))
  assert_that(is.data.frame(newdata))
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim != floor(nsim)) {
    cli_abort("`nsim` must be one finite, positive whole number.")
  }
  if (nsim > .Machine$integer.max) {
    cli_abort("`nsim` must be less than or equal to `.Machine$integer.max`.")
  }
  nsim <- as.integer(nsim)
  for (x in c("return_tmb_report", "silent")) {
    value <- get(x)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      cli_abort("`{x}` must be one non-missing logical value.")
    }
  }
  if (!is.null(sims_var) &&
      (!is.character(sims_var) || length(sims_var) != 1L ||
        is.na(sims_var) || !nzchar(sims_var))) {
    cli_abort("`sims_var` must be `NULL` or one non-empty character value.")
  }
  if ((!is.numeric(sim_re) && !is.logical(sim_re)) || length(sim_re) != 6L ||
      anyNA(sim_re) || !all(sim_re %in% c(0, 1))) {
    cli_abort("`sim_re` must contain exactly six non-missing 0 or 1 values.")
  }
  sim_re <- as.integer(sim_re)

  reinitialize(object)

  if (object$time == "_sdmTMB_time")
    cli_abort("Please refit the sdmTMB model with the 'time' argument specified.")

  if (!object$time %in% names(newdata)) {
    cli_abort("The time column `{object$time}` is missing from `newdata`.")
  }
  time_extension <- project_time_extension(object$time_lu, newdata[[object$time]])
  nproj <- length(time_extension$time_from_data)
  newdata[[object$time]] <- time_extension$canonical_time

  if (all(object$spatiotemporal == "off") && is.null(object$time_varying)) {
    cli_inform("No spatiotemporal or time-varying structures found. Proceeding with projection anyways.")
  }

  uncertainty <- match.arg(uncertainty)
  ee <- object$tmb_obj$env
  lpb <- ee$last.par.best

  if (uncertainty == "both") {
    if (has_no_random_effects(object)) {
      msg <- c("This model has no random effects.",
        "Sampling only from the fixed effects.")
      cli_inform(msg)
      lp <- t(mvtnorm::rmvnorm(n = nsim, mean = lpb, sigma = object$sd_report$cov.fixed))
    } else {
      lp <- rmvnorm_prec(lpb, object$sd_report, nsim)
    }
  } else if (uncertainty == "random") {
    if (has_no_random_effects(object)) {
      msg <- c("This model has no random effects.",
        "Choose a different type of uncertainty.")
      cli_abort(msg)
    }
    lp <- lpb %o% rep(1, nsim)
    mc <- ee$MC(keep = TRUE, n = nsim, antithetic = FALSE)
    lp[ee$random, ] <- attr(mc, "samples")
  } else { ## 'none'
    lp <- lpb %o% rep(1, nsim)
  }

  ## extend time keeping elements of sdmTMB object
  max_year_i <- max(object$time_lu$year_i)
  new_year_i <- max_year_i + seq_len(nproj)
  new_time_from_data <- time_extension$time_from_data
  object$extra_time <- c(object$extra_time, new_time_from_data)
  object$time_lu$sim_projected <- FALSE
  object$time_lu <- rbind(
    object$time_lu,
    data.frame(
      year_i = new_year_i, time_from_data = new_time_from_data,
      extra_time = TRUE, sim_projected = TRUE
    )
  )

  ## generate prediction TMB data list
  p <- predict(object, newdata = newdata, return_tmb_data = TRUE, ...)

  ## move data elements over
  p <- move_proj_to_tmbdat(p, object, newdata)

  ## extend time
  p$n_t <- nrow(object$time_lu)
  p$simulate_t <- as.integer(object$time_lu$sim_projected)
  ## sim random effects? order: omega, epsilon, zeta, IID, RW, smoothers
  p$sim_re <- sim_re

  ## parameters: add zeros as needed to all time-based parameters
  pars <- get_pars(object)
  n_m <- if (is_delta(object)) 2L else 1L
  n_s <- dim(pars$epsilon_st)[1]

  new_eps <- array(0, c(n_s, nproj, n_m))
  n_time_varying <- dim(pars$b_rw_t)[2]
  new_b_rw_t <- array(0, c(nproj, n_time_varying, n_m))
  pars$epsilon_st <- abind::abind(pars$epsilon_st, new_eps, along = 2)
  pars$b_rw_t <- abind::abind(pars$b_rw_t, new_b_rw_t, along = 1)

  map <- object$tmb_map
  if ("b_rw_t" %in% names(map)) {
    if (any(!is.na(map$b_rw_t))) {
      cli_abort("Function not set up yet for non-NA mapping of `b_rw_t`.")
    }
    map$b_rw_t <- factor(rep(NA, length(as.numeric(pars$b_rw_t))))
  }

  delta <- is_delta(object)
  n_active_st <- sum(object$spatiotemporal != "off")
    map$epsilon_st <- array(
      seq_len(length(pars$epsilon_st)),
      dim = dim(pars$epsilon_st)
    )
    for (i in which(object$spatiotemporal == "off")) {
      map$epsilon_st[, , i] <- NA
      new_eps[, , i] <- NA
    }
    map$epsilon_st <- as.factor(map$epsilon_st)
    new_eps <- rep(0, length(new_eps[!is.na(new_eps)]))

  ## rebuild TMB object
  if (!silent) cli::cli_inform("Rebuilding TMB object with TMB::MakeADFun()")
  obj <- TMB::MakeADFun(
    data = p,
    profile = object$control$profile,
    parameters = pars,
    map = map,
    random = object$tmb_random,
    DLL = "sdmTMB",
    silent = TRUE
  )

  ## do simulations
  if (!silent) cli::cli_progress_bar("Simulating projections", total = nsim)
  ret <- list()
  for (i in seq_len(nsim)) {
    if (!silent) cli::cli_progress_update()
    lpx <- lp[, i, drop = TRUE]
    if (!is.null(object$time_varying)) { ## pad time-varying random effects
      lpx <- insert_pars(
        lpx, "b_rw_t", .n = length(as.vector(new_b_rw_t)),
        n_groups = n_time_varying * n_m
      )
    }
    if (length(as.vector(new_eps))) { ## pad spatiotemporal random effects
      lpx <- insert_pars(
        lpx, "epsilon_st", .n = length(as.vector(new_eps)),
        n_groups = n_active_st
      )
    }
    ret[[i]] <- obj$simulate(par = lpx)
  }
  if (!silent) cli::cli_progress_done()
  if (return_tmb_report) {
    return(ret)
  }
  if (!is.null(sims_var)) {
    return(extract_project_sims(ret, sims_var))
  }

  out <- list()
  if (delta) {
    element_names <- c("est1", "est2", "epsilon_st1", "epsilon_st2")
    element_internal <- c("eta_i", "eta_i", "epsilon_st_A_vec", "epsilon_st_A_vec")
    linear_predictor <- c(1L, 2L, 1L, 2L)
  } else {
    element_names <- c("est", "epsilon_st")
    element_internal <- c("eta_i", "epsilon_st_A_vec")
    linear_predictor <- c(1L, 1L)
  }
  for (i in seq_along(element_names)) {
    eni <- element_names[i]
    out[[eni]] <- lapply(ret, \(x) x[[element_internal[i]]][, linear_predictor[i]])
    out[[eni]] <- do.call(cbind, out[[eni]])
  }
  if (all("off" == object$spatiotemporal)) {
    out$epsilon_st1 <- out$epsilon_st2 <- out$epsilon_st <- NULL
  }
  out
}

insert_pars <- function(par, nm, .n, n_groups = 1L) {
  rn <- names(par)
  index <- which(rn == nm)
  if (!length(index)) cli_abort("Parameter `{nm}` was not found.")
  if (!identical(index, seq.int(min(index), max(index)))) {
    cli_abort("Parameter `{nm}` is not stored in one contiguous block.")
  }
  first <- min(index)
  last <- max(index)
  target <- par[index]
  if (length(.n) != 1L || is.na(.n) || !is.finite(.n) || .n < 0 ||
      .n != floor(.n)) {
    cli_abort("The padding length must be one non-negative whole number.")
  }
  if (length(n_groups) != 1L || is.na(n_groups) || !is.finite(n_groups) ||
      n_groups < 1L || n_groups != floor(n_groups)) {
    cli_abort("`n_groups` must be one positive whole number.")
  }
  n_groups <- as.integer(n_groups)
  fill <- rep(0, .n)
  names(fill) <- rep(nm, length(fill))
  if (length(target) %% n_groups || length(fill) %% n_groups) {
    cli_abort("Parameter and padding lengths must divide evenly among parameter groups.")
  }
  target_per_group <- length(target) / n_groups
  fill_per_group <- length(fill) / n_groups
  replacement <- unlist(lapply(seq_len(n_groups), function(i) {
    target_index <- (i - 1L) * target_per_group + seq_len(target_per_group)
    fill_index <- (i - 1L) * fill_per_group + seq_len(fill_per_group)
    c(target[target_index], fill[fill_index], use.names = TRUE)
  }), use.names = FALSE)
  names(replacement) <- rep(nm, length(replacement))
  before <- if (first > 1L) par[seq_len(first - 1L)] else par[0L]
  after <- if (last < length(par)) par[seq.int(last + 1L, length(par))] else par[0L]
  c(before, replacement, after, use.names = TRUE)
}

project_time_extension <- function(time_lu, new_time) {
  fitted_time <- time_lu$time_from_data
  if (!is.numeric(fitted_time) || !is.numeric(new_time)) {
    cli_abort("Projection requires a numeric time column.")
  }
  if (!length(fitted_time) || anyNA(fitted_time) || any(!is.finite(fitted_time))) {
    cli_abort("The fitted time values must contain finite, non-missing values.")
  }
  if (!length(new_time) || anyNA(new_time) || any(!is.finite(new_time))) {
    cli_abort("The projection time column must contain finite, non-missing values.")
  }

  tolerance <- sqrt(.Machine$double.eps) *
    max(1, abs(c(fitted_time, new_time)))
  collapse_time <- function(x) {
    x <- sort(x)
    x[c(TRUE, diff(x) > tolerance)]
  }
  match_time <- function(x, table) {
    vapply(x, function(value) {
      distance <- abs(table - value)
      index <- which.min(distance)
      if (distance[[index]] <= tolerance) index else NA_integer_
    }, integer(1L))
  }

  fitted_time <- collapse_time(fitted_time)
  requested_time <- collapse_time(new_time)
  last_fitted <- max(fitted_time)
  fitted_match <- match_time(requested_time, fitted_time)
  unknown_historical <- requested_time[
    is.na(fitted_match) & requested_time < last_fitted - tolerance
  ]
  if (length(unknown_historical)) {
    cli_abort(c(
      "`newdata` contains time values that are neither fitted nor future values.",
      "i" = "Unknown values: {paste(unknown_historical, collapse = ', ')}."
    ))
  }
  future_time <- requested_time[
    is.na(fitted_match) & requested_time > last_fitted + tolerance
  ]
  if (!length(future_time)) {
    cli_abort(c(
      "No new time elements for projection were found in `newdata`.",
      "i" = "Use `predict(..., nsim = ...)` to simulate without projecting forward."
    ))
  }
  if (length(fitted_time) > 1L) {
    spacing <- diff(fitted_time)
    cadence <- spacing[[1L]]
    if (cadence <= 0 || any(abs(spacing - cadence) > tolerance)) {
      cli_abort(c(
        "The fitted time values do not have a constant cadence.",
        "i" = "Refit with missing time slices supplied through `extra_time`."
      ))
    }
  } else {
    if (length(future_time) < 2L) {
      cli_abort(c(
        "A model with only one fitted time value needs at least two future time values to infer the cadence.",
        "i" = "Supply an additional fitted or `extra_time` time slice."
      ))
    }
    future_spacing <- diff(future_time)
    cadence <- future_spacing[[1L]]
    if (cadence <= 0 || any(abs(future_spacing - cadence) > tolerance)) {
      cli_abort(c(
        "The future projection values do not have a constant cadence.",
        "i" = "Supply an additional fitted or `extra_time` time slice."
      ))
    }
  }
  steps <- (future_time - last_fitted) / cadence
  rounded_steps <- round(steps)
  if (any(abs(future_time - (last_fitted + rounded_steps * cadence)) > tolerance)) {
    cli_abort(c(
      "Future time values cannot be aligned with a constant cadence.",
      "i" = "Supply an additional fitted or `extra_time` time slice to define the cadence."
    ))
  }
  nproj <- max(as.integer(rounded_steps))
  generated_time <- last_fitted + cadence * seq_len(nproj)
  full_time <- collapse_time(c(fitted_time, generated_time))
  canonical_match <- match_time(new_time, full_time)
  if (anyNA(canonical_match)) {
    cli_abort("Could not match all projection time values to the fitted or generated time grid.")
  }
  list(
    time_from_data = generated_time,
    canonical_time = full_time[canonical_match],
    cadence = cadence
  )
}

extract_project_sims <- function(reports, sims_var) {
  available <- names(reports[[1L]])
  if (!sims_var %in% available) {
    cli_abort(c(
      "`sims_var = \"{sims_var}\"` was not found in the TMB simulation report.",
      "i" = "Available elements include: {paste(available, collapse = ', ')}."
    ))
  }
  values <- lapply(reports, `[[`, sims_var)
  first <- values[[1L]]
  first_dim <- dim(first)
  value_length <- length(first)
  if (!all(vapply(values, length, integer(1L)) == value_length) ||
      !all(vapply(values, function(x) identical(dim(x), first_dim), logical(1L)))) {
    cli_abort("TMB report element `{sims_var}` changed dimensions between simulations.")
  }
  out_dim <- c(if (is.null(first_dim)) value_length else first_dim, length(values))
  out <- array(unlist(values, use.names = FALSE), dim = out_dim)
  first_dimnames <- dimnames(first)
  if (is.null(first_dim)) first_dimnames <- list(names(first))
  if (!is.null(first_dimnames)) dimnames(out) <- c(first_dimnames, list(NULL))
  out
}

move_proj_to_tmbdat <- function(x, object, newdata, called_by_simulate = FALSE, size = NULL) {
  x$do_predict <- 0L
  x$year_i <- x$proj_year
  ## x$A_st <- x$proj_mesh
  ## .cpp uses unique locations in projection but not in fitting:
  xy_cols <- object$spde$xy_cols
  proj_mesh <- fmesher::fm_basis(object$spde$mesh, loc = as.matrix(newdata[, xy_cols, drop = FALSE]))
  x$A_st <- proj_mesh
  x$A_spatial_index <- seq_len(dim(proj_mesh)[1]) - 1L
  x$X_threshold <- x$proj_X_threshold
  x$X_ij <- x$proj_X_ij
  x$X_rw_ik <- x$proj_X_rw_ik
  x$z_i <- x$proj_z_i
  x$Zs <- x$proj_Zs
  x$Xs <- x$proj_Xs
  x$Zt_list <- x$Zt_list_proj
  x$offset_i <- x$proj_offset_i
  n_m <- length(x$X_ij) ## n linear predictor [m]odels
  x$y_i <- matrix(NA, ncol = n_m, nrow = nrow(x$proj_X_ij[[1]])) # fake
  x$weights_i <- rep(1, nrow(x$y_i)) # fake: FIXME: bring in?
  x$area_i <- rep(1, nrow(x$y_i)) # fake FIXME: bring in for index??

  if (!called_by_simulate) {
    unique_size <- unique(x$size)
    if (length(unique_size) != 1L) {
      cli_abort("This function hasn't been set up to work with binomial size specified yet.")
    }
    x$size <- rep(1, nrow(x$y_i)) # FIXME: bring in?
  }
  if (!is.null(size)) {
    if (length(size) != nrow(x$y_i)) cli_abort("size argument doesn't match rows of newdata")
    x$size <- size
  }
  if (!all(x$size == 1L)) {
    if (is.null(size)) {
      msg <- c(
        "It looks like the size/trials were specified for a binomial family with the weights argument.",
        "When simulating with newdata all sizes/trials will be assumed to be 1 unless you specify a vector via the `size` argument."
      )
      cli_inform(msg)
      x$size <- rep(1, nrow(x$y_i))
    }
  } else { # all 1s, but newdata; expand as needed silently
    x$size <- rep(1, nrow(x$y_i))
  }
  x$do_predict <- 0L

  # nullify large data objects that are no longer needed:
  x$proj_X_ij <- list(matrix(0, ncol = 1, nrow = 1))
  x$proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy
  x$proj_mesh <- Matrix::Matrix(c(0, 0, 2:0), 3, 5) # dummy
  x$proj_Zs <- list()
  x$proj_Xs <- matrix(nrow = 0L, ncol = 0L)
  x$proj_lon <- 0
  x$proj_lat <- 0

  x
}
