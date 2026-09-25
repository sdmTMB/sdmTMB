# Deterministic, distinct values for every estimated parameter cell, so
# misplaced components or array axes change the objective. Cells fixed by
# `map` keep their values, as in a real fit.
rtmb_test_parameters <- function(parameters, map = list(), scale = 0.2) {
  for (name in names(parameters)) {
    n <- length(parameters[[name]])
    if (!n) next
    free <- if (is.null(map[[name]])) rep(TRUE, n) else !is.na(map[[name]])
    values <- seq(-scale, scale, length.out = n)
    parameters[[name]][free] <- values[free]
  }
  parameters
}

# A deterministic perturbation of a parameter vector, small enough to stay
# within each parameter's domain.
rtmb_perturb <- function(par) par + 0.05 * sin(seq_along(par))

# ADREPORT order differs between backends, so `sdreport()` summaries are
# compared by name. `split()` keeps the order within each name and the number
# of rows per name, so swapped elements or repeated reports still fail.
sdreport_by_name <- function(sd, columns = "Estimate") {
  lapply(stats::setNames(columns, columns),
    function(column) split(unname(sd[, column]), rownames(sd)))
}

# Compare the `sdreport()` estimates and standard errors of two fits.
expect_sdreport_matches <- function(actual, expected, tolerance = 1e-5,
                                    info = NULL) {
  columns <- c("Estimate", "Std. Error")
  for (type in c("fixed", "report")) {
    expect_equal(
      sdreport_by_name(summary(actual$sd_report, type), columns),
      sdreport_by_name(summary(expected$sd_report, type), columns),
      tolerance = tolerance, info = paste(info, type))
  }
}

# Compare the RTMB and C++ objectives: value and gradient at the given and a
# perturbed point of the same tapes, every report's value and dimension, and
# `sdreport()` estimates. With random effects, the joint objectives are also
# compared, so the latent arrays are exercised rather than optimized away.
expect_rtmb_matches_tmb <- function(data, parameters, map, random,
                                    info = NULL, sdreport = TRUE) {
  objectives <- function(random) list(
    cpp = make_sdmTMB_adfun(data, parameters, map, random, backend = "tmb"),
    rtmb = make_sdmTMB_adfun(data, parameters, map, random, backend = "rtmb"))
  expect_same_objective <- function(obj, par, label) {
    expect_equal(obj$rtmb$fn(par), obj$cpp$fn(par), tolerance = 1e-7,
      info = label)
    expect_equal(obj$rtmb$gr(par), obj$cpp$gr(par), ignore_attr = TRUE,
      tolerance = 1e-6, info = label)
  }
  obj <- objectives(random)
  cpp <- obj$cpp
  rt <- obj$rtmb
  start <- cpp$par
  expect_same_objective(obj, start, info)
  expected <- cpp$report()
  actual <- rt$report()
  expect_setequal(names(actual), names(expected))
  for (name in intersect(names(expected), names(actual))) {
    expect_equal(as.vector(actual[[name]]), as.vector(expected[[name]]),
      tolerance = 1e-6, info = paste(info, name))
    expect_equal(dim(actual[[name]]), dim(expected[[name]]),
      info = paste(info, name))
  }
  if (sdreport) {
    sd_cpp <- suppressWarnings(summary(sdreport_sdmTMB(cpp), "report"))
    sd_rt <- suppressWarnings(summary(sdreport_sdmTMB(rt), "report"))
    expect_equal(sdreport_by_name(sd_rt), sdreport_by_name(sd_cpp),
      tolerance = 1e-6, info = info)
  }
  expect_same_objective(obj, rtmb_perturb(start), paste(info, "perturbed"))
  if (length(random)) {
    joint <- objectives(NULL)
    for (par in list(joint$cpp$par, rtmb_perturb(joint$cpp$par))) {
      expect_same_objective(joint, par, paste(info, "joint"))
    }
  }
  invisible(obj)
}

# Check fitted rows and, unless `project = FALSE`, projected rows.
expect_rtmb_fit_data_matches <- function(fit, newdata, info,
                                         project = TRUE, edit = identity) {
  p <- edit(rtmb_test_parameters(fit$tmb_params, fit$tmb_map))
  expect_rtmb_matches_tmb(fit$tmb_data, p, fit$tmb_map, fit$tmb_random,
    info = info)
  if (project) {
    projection <- predict(fit, newdata = newdata, offset = newdata$off,
      return_tmb_data = TRUE)
    expect_rtmb_matches_tmb(projection, p, fit$tmb_map, fit$tmb_random,
      info = paste(info, "projection"))
  }
}
