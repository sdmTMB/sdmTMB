#' Run extra optimization on an already fitted object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @param object An object from [sdmTMB()].
#' @param nlminb_loops How many extra times to run [stats::nlminb()]
#'   optimization. Sometimes restarting the optimizer at the previous best
#'   values aids convergence.
#' @param newton_loops How many extra Newton optimization loops to try with
#'   [stats::optimHess()]. Sometimes aids convergence.
#'
#' @return An updated model fit of class `sdmTMB`.
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
  # FIXME: DRY; use this function in sdmTMB()

  new_obj <- object
  old_par <- object$model$par
  new_obj$tmb_obj$fn(old_par) # initialize the TMB object
  tmb_opt <- new_obj$model
  for (i in seq_len(nlminb_loops)) {
    temp <- new_obj$model[c("iterations", "evaluations")]
    tmb_opt <- stats::nlminb(
      start = tmb_opt$par,
      objective = new_obj$tmb_obj$fn,
      gradient = new_obj$tmb_obj$gr,
      control = object$nlminb_control,
      lower = object$lower,
      upper = object$upper
    )
    tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
    tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
  }

  if (any(!is.infinite(object$upper)) || any(!is.infinite(object$lower))) {
    if (newton_loops > 0) {
      cli_inform("Upper or lower limits were set. `stats::optimHess()` will ignore these limits. Set `control = sdmTMBcontrol(newton_loops = 0)` to avoid the `stats::optimHess()` optimization if desired.")
    }
  }

  tmb_opt <- run_newton_loops(newton_loops = newton_loops, opt = tmb_opt, obj = new_obj$tmb_obj, silent = FALSE)
  new_obj$model <- tmb_opt
  new_obj$sd_report <- TMB::sdreport(new_obj$tmb_obj,
    getJointPrecision = "jointPrecision" %in% names(object$sd_report))
  conv <- get_convergence_diagnostics(new_obj$sd_report)
  new_obj$gradients <- conv$final_grads
  new_obj$bad_eig <- conv$bad_eig
  new_obj
}

run_newton_loops <- function(newton_loops, opt, obj, silent = TRUE) {
  if (newton_loops > 0) {
    if (!silent) cli_inform("attempting to improve convergence with Newton update(s)")
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(obj$gr(opt$par))
      if (max(abs(g)) < 1e-9) {
        if (!silent) cli_inform(c("maximum absolute gradient is already < 1e-9;",
          "skipping any remaining Newton updates for speed"))
        break
      }
      h <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr)
      new_par <- tryCatch(
        opt$par - solve(h, g),
        error = function(e) {
          if (!silent) cli_inform("Hessian is singular; skipping Newton update(s)")
          NULL
        }
      )
      if (is.null(new_par)) break
      new_objective <- obj$fn(new_par) # also updates obj$env$last.par and obj$env$last.par.best!
      if (new_objective < opt$objective) {
        if (!silent) cli_inform("accepting parameters from Newton update")
        opt$par <- new_par
        opt$objective <- new_objective
      } else {
        if (!silent) cli_inform(c("retaining parameters from before Newton update",
          "and skipping further Newton updates"))
        break
      }
    }
  }
  opt
}
