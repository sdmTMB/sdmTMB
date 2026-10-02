#' Run extra optimization on an already fitted object
#'
#' @param object An object from [sdmTMB()].
#' @param nlminb_loops How many extra times to run [stats::nlminb()]
#'   optimization. Sometimes restarting the optimizer at the previous best
#'   values aids convergence.
#' @param newton_loops How many extra Newton optimization loops to try with
#'   [stats::optimHess()]. Sometimes aids convergence.
#'
#' @return An updated model fit of class `sdmTMB`. `object` is left unchanged.
#' @export
#'
#' @examples
#' # Run extra optimization steps to help convergence:
#' # (Not typically needed)
#' fit <- sdmTMB(density ~ 0 + poly(depth, 2) + as.factor(year),
#'   data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie())
#' fit_1 <- run_extra_optimization(fit, newton_loops = 1)
#' max(fit$gradients)
#' max(fit_1$gradients)
run_extra_optimization <- function(object,
  nlminb_loops = 0,
  newton_loops = 1) {
  new_obj <- object
  # a fresh objective so optimizing doesn't modify `object`'s TMB environment
  new_obj$tmb_obj <- remake_tmb_obj(object)
  tmb_opt <- run_nlminb_loops(
    nlminb_loops = nlminb_loops, opt = object$model, obj = new_obj$tmb_obj,
    lower = object$lower, upper = object$upper,
    control = object$nlminb_control, silent = FALSE,
    suppress_warnings = isTRUE(object$control$suppress_nlminb_warnings)
  )
  tmb_opt <- run_newton_loops(
    newton_loops = newton_loops, opt = tmb_opt, obj = new_obj$tmb_obj,
    silent = FALSE, lower = object$lower, upper = object$upper
  )
  new_obj$model <- tmb_opt
  new_obj$sd_report <- sdreport_sdmTMB(new_obj$tmb_obj,
    getJointPrecision = "jointPrecision" %in% names(object$sd_report))
  conv <- get_convergence_diagnostics(new_obj$sd_report)
  new_obj$gradients <- conv$final_grads
  new_obj$bad_eig <- conv$bad_eig
  new_obj$pos_def_hessian <- new_obj$sd_report$pdHess
  lpb <- new_obj$tmb_obj$env$last.par.best
  new_obj$parlist <- new_obj$tmb_obj$env$parList(par = lpb)
  new_obj$last.par.best <- lpb
  new_obj
}

# Rebuild a fit's TMB objective and evaluate it at the fitted parameters. The
# result shares no environment with `object$tmb_obj`.
remake_tmb_obj <- function(object) {
  if (is.null(object$parlist)) {
    cli_abort("This model was fit with an sdmTMB version that is too old; please refit it.")
  }
  obj <- make_sdmTMB_adfun(
    data = object$tmb_data,
    parameters = object$parlist,
    map = object$tmb_map,
    random = object$tmb_random,
    backend = backend_sdmTMB(object),
    profile = object$control$profile
  )
  if (isTRUE(object$control$normalize) &&
      isTRUE(object$tmb_data$normalize_in_r == 1L)) {
    obj <- TMB::normalize(obj, flag = "flag", value = 0)
  }
  obj$env$beSilent()
  obj$fn(object$model$par) # restores last.par.best etc.
  obj
}

# `suppress_nlminb_warnings`: nlminb() warns when it probes parameters where the
# objective is `NA`/`NaN`; it treats these as `Inf` and steps back.
maybe_suppress_warnings <- function(suppress) {
  if (isTRUE(suppress)) suppressWarnings else identity
}

# First of two phases: fit the fixed effects with random fields off (and
# censored Poisson as Poisson for stability), returning starting values for
# the full model. Falls back to `tmb_params` if this phase fails.
fit_first_phase <- function(tmb_data, tmb_params, tmb_map, profile, backend,
                            lower, upper, mesh, nlminb_control, silent,
                            suppress_warnings = FALSE) {
  tmb_data$no_spatial <- 1L
  tmb_data$include_spatial <- integer(ncol(tmb_data$component_active)) # per component
  censored <- tmb_data$component_active == 1L &
    tmb_data$family_code == .valid_family[["censored_poisson"]]
  tmb_data$family_code[censored] <- as.integer(.valid_family[["poisson"]])
  obj <- make_sdmTMB_adfun(data = tmb_data, parameters = tmb_params,
    profile = profile, map = tmb_map, backend = backend, silent = silent)
  lim <- set_limits(obj, lower = lower, upper = upper, mesh = mesh,
    spatial_model = tmb_data$spatial_model, silent = TRUE)
  opt <- tryCatch(
    maybe_suppress_warnings(suppress_warnings)(stats::nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr,
      lower = lim$lower, upper = lim$upper, control = nlminb_control)),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective) || !all(is.finite(opt$par))) {
    if (!silent) cli_inform("First optimization phase failed; using default starting values")
    return(tmb_params)
  }
  out <- obj$env$parList(opt$par) # no random effects in this phase
  # phase-one threshold estimates often cause optimization problems
  out$b_threshold <- tmb_params$b_threshold
  out
}

# Restart nlminb() from `opt$par` up to `nlminb_loops` times, stopping if a
# restart fails to improve the objective.
run_nlminb_loops <- function(nlminb_loops, opt, obj, lower, upper, control,
                             silent = TRUE, suppress_warnings = FALSE) {
  if (nlminb_loops < 1 || !length(opt$par)) return(opt)
  if (!silent) cli_inform("running extra nlminb optimization")
  quietly <- maybe_suppress_warnings(suppress_warnings)
  for (i in seq_len(nlminb_loops)) {
    new_opt <- tryCatch(
      quietly(stats::nlminb(
        start = opt$par, objective = obj$fn, gradient = obj$gr,
        control = control, lower = lower, upper = upper
      )),
      error = function(e) NULL
    )
    if (is.null(new_opt) || !is.finite(new_opt$objective) ||
        new_opt$objective > opt$objective) {
      if (!silent) cli_inform("extra nlminb optimization did not improve the objective; retaining previous parameters")
      obj$fn(opt$par) # reset the objective's state to the retained parameters
      break
    }
    new_opt$iterations <- new_opt$iterations + opt$iterations
    new_opt$evaluations <- new_opt$evaluations + opt$evaluations
    opt <- new_opt
  }
  opt
}

run_newton_loops <- function(newton_loops, opt, obj, silent = TRUE,
                             lower = NULL, upper = NULL) {
  if (newton_loops < 1 || !length(opt$par)) return(opt)
  inform <- function(...) if (!silent) cli_inform(c(...))
  inform("attempting to improve convergence with Newton update(s)")
  for (i in seq_len(newton_loops)) {
    g <- as.numeric(obj$gr(opt$par))
    if (!all(is.finite(g))) {
      inform("non-finite gradient; skipping Newton update(s)")
      break
    }
    if (max(abs(g)) < 1e-9) {
      inform("maximum absolute gradient is already < 1e-9;",
        "skipping any remaining Newton updates for speed")
      break
    }
    step <- tryCatch(
      solve(stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr), g),
      error = function(e) NULL
    )
    if (is.null(step) || !all(is.finite(step))) {
      inform("Hessian could not be computed or inverted; skipping Newton update(s)")
      break
    }
    new_opt <- try_newton_step(opt, obj, step, g, lower, upper)
    if (is.null(new_opt)) {
      inform("retaining parameters from before Newton update",
        "and skipping further Newton updates")
      break
    }
    inform("accepting parameters from Newton update")
    opt <- new_opt
  }
  # optimHess() and rejected steps leave the objective at other parameters
  obj$fn(opt$par)
  sync_last_par_best(obj, opt$par)
  opt
}

# TMB and RTMB keep the conditional random-effect mode in `last.par`, but only
# update `last.par.best` when the objective strictly decreases. A Newton step
# that is accepted at an equal objective would otherwise leave `last.par.best`
# (used by `sdreport()`, `predict()`, etc.) at the previous fixed effects. Keep
# the reported fixed effects and conditional mode in sync after a Newton update.
sync_last_par_best <- function(obj, par) {
  env <- obj$env
  last_par <- env$last.par
  if (is.null(env$random)) {
    env$last.par.best <- par
  } else {
    last_par[-env$random] <- par
    env$last.par.best <- last_par
  }
  invisible(NULL)
}

# Take the Newton step, halving it until it stays within the limits, does not
# increase the objective beyond numerical tolerance, and improves the
# gradient. Returns `NULL` if no step length works.
try_newton_step <- function(opt, obj, step, gradient, lower, upper,
                            max_halvings = 5L) {
  for (k in 0:max_halvings) {
    new_par <- opt$par - step / 2^k
    within_bounds <-
      (is.null(lower) || all(new_par >= lower)) &&
      (is.null(upper) || all(new_par <= upper))
    if (!within_bounds) next
    new_objective <- as.numeric(obj$fn(new_par))
    new_gradient <- tryCatch(as.numeric(obj$gr(new_par)), error = function(e) NULL)
    objective_tolerance <- 100 * .Machine$double.eps *
      max(1, abs(opt$objective))
    gradient_improved <- !is.null(new_gradient) &&
      all(is.finite(new_gradient)) &&
      max(abs(new_gradient)) < max(abs(gradient))
    objective_improved <- is.finite(new_objective) &&
      new_objective <= opt$objective + objective_tolerance
    if (gradient_improved && objective_improved) {
      opt$par <- new_par
      opt$objective <- new_objective
      return(opt)
    }
  }
  NULL
}
