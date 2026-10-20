# Precision check for the RTMB censored beta-binomial likelihood. Rows with
# `cens_direct == 0` sum the interval or its complement, whichever is shorter
# (see rtmb_dcensbetabinom()). The complement loses accuracy when the interval
# probability is small, so after fitting, rows that disagree with the direct
# sum are switched to it and the model is refit.

# Initial `cens_direct` flags: 1 sums the interval directly. Rows with zero
# likelihood weight, such as held-out cross-validation rows, are summed
# directly because they may fit poorly and are only evaluated for scoring.
.censored_direct_init <- function(weights, n, method) {
  if (identical(method, "direct")) return(rep(1L, n))
  if (is.null(weights)) return(rep(0L, n))
  as.integer(weights %in% 0)
}

# Rows not yet summed directly whose observation likelihood at the estimate
# differs from the direct sum by more than `tol`, or isn't finite.
check_censored_betabinomial <- function(obj, data, tol = 1e-8) {
  check <- data$cens_direct == 0L
  if (!any(check)) return(integer(0))
  par <- obj$env$parList(par = obj$env$last.par.best)
  shorter <- rtmb_report_values(data, par)$jnll_obs
  data$cens_direct[] <- 1L
  direct <- rtmb_report_values(data, par)$jnll_obs
  which(check & !(abs(shorter - direct) <= tol) %in% TRUE)
}

# Switch rows that fail the check to the direct sum and refit from the current
# estimates until the check passes. Returns the objective, optimization
# result, and data with the final `cens_direct` flags.
refit_censored_betabinomial <- function(obj, opt, data, map, random, lim,
                                        control, profile, silent,
                                        suppress_warnings = FALSE) {
  repeat {
    rows <- check_censored_betabinomial(obj, data)
    if (!length(rows)) break
    if (!silent) {
      cli_inform("Refitting with the direct censored beta-binomial sum for {length(rows)} row{?s} that failed the precision check.")
    }
    data$cens_direct[rows] <- 1L
    obj <- make_sdmTMB_adfun(data,
      parameters = obj$env$parList(par = obj$env$last.par.best),
      map = map, random = random, profile = profile, backend = "rtmb",
      silent = silent)
    opt <- maybe_suppress_warnings(suppress_warnings)(stats::nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr,
      lower = lim$lower, upper = lim$upper, control = control))
  }
  list(obj = obj, opt = opt, data = data)
}
