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

# Frozen objective values. The C++ TMB template was the oracle for the RTMB
# objective; its `fn()` and `gr()` at deterministic parameters (plus a few
# values outside the likelihood: projections, deviance residuals, derived
# indices) are stored as `deparse` snapshots in `_snaps/`, and RTMB is compared
# against them. They were recorded from the C++ template before it was removed
# (sdmTMB 1.2.0 still has both backends). After an intended change, rerun the
# test and accept the new RTMB values with `testthat::snapshot_review()` or
# `testthat::snapshot_accept()`, and say why in the commit.
frozen_tolerance <- 1e-8

expect_frozen_snapshot <- function(key, actual, tolerance = frozen_tolerance) {
  # 10 significant digits resolve differences well below the tolerance.
  actual <- rapply(actual, signif, classes = "numeric", how = "replace",
    digits = 10)
  expect_snapshot_value(c(list(key = key), actual), style = "deparse",
    tolerance = tolerance)
}

# Short vectors are stored whole. Long ones are reduced to their length, sum,
# absolute sum, and a position-weighted sum, so swapped elements still fail.
frozen_summary <- function(x) {
  if (!is.numeric(x)) return(NULL)
  x <- as.vector(x)
  if (length(x) <= 10L) return(x)
  c(length(x), sum(x), sum(abs(x)), sum(sin(seq_along(x)) * x))
}

# ADREPORT order differed between the backends, so `sdreport()` summaries are
# summarized by name, keeping the order within each name.
frozen_sdreport <- function(obj) {
  sd <- suppressWarnings(summary(RTMB::sdreport(obj), "report"))
  sd <- sd[!rownames(sd) %in% c("log_range_min", "log_range_max"), ,
    drop = FALSE]
  lapply(split(seq_len(nrow(sd)), rownames(sd)), function(i)
    list(est = frozen_summary(sd[i, "Estimate"]),
      se = frozen_summary(sd[i, "Std. Error"])))
}

# RTMB computes deviance residuals only on request, off the tape.
frozen_devresid <- function(obj, data) {
  rtmb_report_values(data, obj$env$parList(), deviance = TRUE)$devresid
}

frozen_fn_gr <- function(obj, par) {
  list(fn = as.numeric(obj$fn(par)), gr = frozen_summary(as.vector(obj$gr(par))))
}

# `fn()` and `gr()` at the start, the named `reports` (there and, with random
# effects, at its inner mode), and optionally the `sdreport()` estimates and
# SEs. With random effects and `joint = TRUE`, the joint objective is included
# too, so the latent arrays are exercised rather than optimized away.
# `objective = FALSE` keeps only the reports and `sdreport()` values, e.g., for
# projection data, whose likelihood is the fitted data's.
frozen_objective_values <- function(data, parameters, map = list(),
                                    random = NULL, reports = character(0),
                                    sdreport = FALSE, joint = TRUE,
                                    objective = TRUE) {
  make <- function(random) make_sdmTMB_adfun(data, parameters, map, random)
  obj <- make(random)
  # Evaluated even when not stored: reports use the inner mode at the start.
  fn_gr <- frozen_fn_gr(obj, obj$par)
  out <- if (objective) fn_gr else list()
  # The TMB template had no NB1 deviance residuals (family code 10).
  if (any(data$family_code == 10L)) reports <- setdiff(reports, "devresid")
  if (length(reports)) {
    rep <- obj$report()
    if ("devresid" %in% reports) rep$devresid <- frozen_devresid(obj, data)
    out$report <- lapply(stats::setNames(reports, reports), function(name)
      frozen_summary(rep[[name]]))
  }
  if (sdreport) out$sdreport <- frozen_sdreport(obj)
  if (objective && joint && length(random)) {
    j <- make(NULL)
    out$joint <- frozen_fn_gr(j, j$par)
  }
  list(obj = obj, values = out)
}

# Compare an RTMB objective with frozen values (see
# `frozen_objective_values()`). Returns the RTMB object.
expect_frozen_objective <- function(key, data, parameters, map = list(),
                                    random = NULL, reports = character(0),
                                    sdreport = FALSE, joint = TRUE,
                                    objective = TRUE) {
  x <- frozen_objective_values(data, parameters, map, random, reports,
    sdreport, joint, objective)
  expect_frozen_snapshot(key, x$values)
  invisible(x$obj)
}

# Check fitted rows and, unless `project = FALSE`, the projected linear
# predictor, against frozen values keyed by `paste(prefix, info)`.
expect_frozen_fit_data <- function(prefix, fit, newdata, info,
                                   project = TRUE, edit = identity) {
  key <- paste(prefix, info)
  p <- edit(rtmb_test_parameters(fit$tmb_params, fit$tmb_map))
  expect_frozen_objective(key, fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  if (project) {
    projection <- predict(fit, newdata = newdata, offset = newdata$off,
      return_tmb_data = TRUE)
    expect_frozen_objective(paste(key, "projection"), projection, p,
      fit$tmb_map, fit$tmb_random, reports = "proj_eta", objective = FALSE)
  }
}

# Compare other deterministic values (a numeric vector or a list of them)
# with frozen ones.
expect_frozen_values <- function(key, values, tolerance = frozen_tolerance) {
  summarize <- function(x) {
    if (is.list(x)) lapply(x, summarize) else frozen_summary(x)
  }
  expect_frozen_snapshot(key, list(values = summarize(values)),
    tolerance = tolerance)
}
