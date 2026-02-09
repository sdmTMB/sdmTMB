#' Project from an \pkg{sdmTMB} model using simulation
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' @description The function enables projecting forward in time from an
#' \pkg{sdmTMB} model using a simulation approach for computational efficiency.
#' This can be helpful for calculating predictive intervals for long
#' projections where including those time elements in `extra_time` during model
#' estimation can be slow.
#'
#' @description Inspiration for this approach comes from the \pkg{VAST} function
#' `project_model()`.
#'
#' @param object A fitted model from [sdmTMB()].
#' @param newdata A new data frame to predict on. Should contain both historical
#'   and any new time elements to predict on.
#' @param nsim Number of simulations.
#' @param uncertainty How to sample uncertainty for the fitted parameters:
#'   `"both"` for the joint fixed and random effect precision matrix,
#'   `"random"` for the random effect precision matrix (holding the fixed
#'   effects at their MLE), or `"none"` for neither.
#' @param silent Silent?
#' @param sims_var Element to extract from the \pkg{TMB} report. Also see
#'   `return_tmb_report`.
#' @param sim_re A vector of `0`s and `1`s representing which random effects to
#'   simulate in the projection. Generally, leave this untouched. Order is:
#'   spatial fields, spatiotemporal fields, spatially varying coefficient
#'   fields, random intercepts, time-varying coefficients, smoothers.
#'   The default is to simulate spatiotemporal fields and time-varying
#'   coefficients, if present.
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
#' If `return_tmb_report = TRUE`, a list of \pkg{TMB} reports from `simulate()`.
#' Run `names()` on the output to see the options.
#' @export
#'
#' @examplesIf ggplot2_installed()
#' \donttest{
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
project <- function(
    object,
    newdata,
    nsim = 1,
    uncertainty = c("both", "random", "none"),
    silent = FALSE,
    sims_var = "eta_i",
    sim_re = c(0, 1, 0, 0, 1, 0),
    return_tmb_report = FALSE,
    ...) {
  assert_that(inherits(object, "sdmTMB"))
  assert_that(length(nsim) == 1L)
  assert_that(is.numeric(nsim))
  assert_that(nsim >= 1)
  # assert_that(length(nproj) == 1L)
  # assert_that(is.numeric(nproj))
  # assert_that(nproj >= 1)
  assert_that(is.data.frame(newdata))
  assert_that(length(sims_var) == 1L)
  assert_that(is.logical(return_tmb_report))
  assert_that(is.logical(silent))
  assert_that(length(sim_re) == 6L)
  assert_that(all(as.integer(sim_re) %in% c(0L, 1L)))

  reinitialize(object)

  if (object$time == "_sdmTMB_time")
    cli_abort("Please refit the sdmTMB model with the 'time' argument specified.")

  ti <- sort(unique(newdata[[object$time]]))
  fitted_time <- object$time_lu$time_from_data
  nproj <- max(ti) - max(fitted_time)
  if (nproj == 0L) {
    msg <- c(
      "No new time elements for projection found in 'newdata'.",
      "Please supply new time elements in 'newdata'.",
      "If you wish to simulate from predictions on new data without supplying new time elements, please use sdmTMB::predict(..., nsim = ...)."
    )
    cli_abort(msg)
  }

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
  new_year_i <- seq(max_year_i + 1, max_year_i + nproj)
  max_time <- max(object$time_lu$time_from_data)
  new_time_from_data <- seq(max_time + 1, max_time + nproj)
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
  n_m <- if (.is_regular_delta_model(object)) 2L else 1L
  n_s <- dim(pars$epsilon_st)[1]

  new_eps <- array(0, c(n_s, nproj, n_m))
  new_b_rw_t <- array(0, c(nproj, 1, n_m))
  pars$epsilon_st <- abind::abind(pars$epsilon_st, new_eps, along = 2)
  pars$b_rw_t <- abind::abind(pars$b_rw_t, new_b_rw_t, along = 1)

  ## FIXME: drop dims to 0 to avoid mapping complexity
  map <- object$tmb_map
  if ("b_rw_t" %in% names(map)) {
    if (any(!is.na(map$b_rw_t))) {
      cli_abort("Function not set up yet for non-NA mapping of `b_rw_t`.")
    }
    map$b_rw_t <- factor(rep(NA, length(as.numeric(pars$b_rw_t))))
  }

  delta <- .is_regular_delta_model(object)
  ## epsilon_st is always in map??
  # if (delta && "off") {
  # if (delta && "off" %in% object$spatiotemporal) {
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
  # }

  ## rebuild TMB object
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
      lpx <- insert_pars(lpx, "b_rw_t", .n = length(as.vector(new_b_rw_t)), delta_model = delta)
    }
    if (length(as.vector(new_eps))) { ## pad spatiotemporal random effects
      lpx <- insert_pars(lpx, "epsilon_st", .n = length(as.vector(new_eps)), delta_model = delta)
    }
    ret[[i]] <- obj$simulate(par = lpx)
  }
  cli::cli_progress_done()
  if (return_tmb_report) {
    return(ret)
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

insert_pars <- function(par, nm, .n, delta_model = FALSE) {
  rn <- names(par)
  first <- min(which(rn == nm))
  last <- max(which(rn == nm))
  npar <- sum(rn == nm)
  if (delta_model) { ## must inject 1st linear predictor mid-way through
    mid <- npar / 2 ## guaranteed to be even
    fill <- rep(0, .n)
    names(fill) <- rep(nm, length(fill))
    fill1 <- fill[seq(1, length(fill) / 2)] ## split fill in 2
    fill2 <- fill[seq(length(fill) / 2 + 1, length(fill))]
    ret <- c(
      par[seq(1, first - 1)], ## up to our param of interest
      par[seq(first, first + mid)], ## 1st half of param of interest
      fill1,
      par[seq(first + mid + 1, last)], ## 2nd half of param of interest
      fill2,
      use.names = TRUE
    )
  } else {
    fill <- rep(0, .n)
    names(fill) <- rep(nm, length(fill))
    ret <- c(par[seq(last)], fill, use.names = TRUE)
  }
  if (rn[length(rn)] != nm) { ## not at end; append on the rest of par
    ret <- c(ret, par[seq(last + 1, length(par))], use.names = TRUE)
    ret
  }
  ret
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
  if (isTRUE(as.logical(x$has_dispersion_model))) {
    if (is.null(x$proj_Xdisp_ij)) {
      cli_abort("Internal error: missing dispersion projection matrix in prediction data.")
    }
    x$Xdisp_ij <- x$proj_Xdisp_ij
  }
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
  x$proj_Xdisp_ij <- matrix(0, ncol = 1, nrow = 1)
  x$proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy
  x$proj_mesh <- Matrix::Matrix(c(0, 0, 2:0), 3, 5) # dummy
  x$proj_Zs <- list()
  x$proj_Xs <- matrix(nrow = 0L, ncol = 0L)
  x$proj_lon <- 0
  x$proj_lat <- 0

  x
}
