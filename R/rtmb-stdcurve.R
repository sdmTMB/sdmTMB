# eDNA standard curves (the `stdcurve()` family) --------------------------
#
# A `stdcurve()` model has two data sets that share parameters:
#
# - field samples (`data`): qPCR cycle thresholds (Ct; 0 = no detection)
#   whose latent log DNA concentration is the usual sdmTMB linear predictor
#   (a log link, so the response-scale mean is concentration);
# - laboratory standards (`stdcurve(standards = )`): Ct values of reactions
#   with known DNA concentration.
#
# Each PCR plate j has four coefficients, the columns of the random-effect
# matrix `std_xi` (n_plate x 4):
#
#   1. Ct intercept    \  positive Ct ~ Normal(xi1 + xi2 * x, phi)
#   2. Ct slope        /
#   3. detection intercept \  logit P(Ct > 0) = xi3 + xi4 * f(x)
#   4. detection slope     /
#
# where x is log concentration (known for standards, the linear predictor for
# field samples) and f() is exp() (`detection = "concentration"`) or the
# identity (`detection = "log_concentration"`). Each column is
# Normal(std_mu[k], exp(ln_std_sigma[k])) across plates, independently (no
# correlation among a plate's coefficients). `phi` (`ln_phi`) is shared by
# standards and field samples.
#
# This is RTMB only: the R side (`.stdcurve_tmb_data()`) adds a `stdcurve`
# list to the C++ data list; `rtmb_prepare()` passes it through; the
# standards and plate random effects enter the objective in
# `rtmb_stdcurve_nll()`; and the field likelihood is the `stdcurve` entry in
# `rtmb_obs_families`, with each row's plate coefficients in its state.

.stdcurve_coef_names <- c("ct_intercept", "ct_slope", "det_intercept",
  "det_slope")

# Data for the standards and the plate index of each field row. `data` is the
# analysis data after rows with missing values are dropped.
.stdcurve_tmb_data <- function(family, data, multi_family = FALSE,
                               dispformula = FALSE) {
  if (multi_family) {
    cli_abort("`stdcurve()` is not supported in multi-family models.")
  }
  if (dispformula) {
    cli_abort("`dispformula` is not supported with `stdcurve()`.")
  }
  std <- family$standards
  plate_col <- family$plate
  plates <- sort(unique(std$plate))
  plate_i <- .stdcurve_plate_index(data, plate_col, plates)
  if (length(plates) < 2L) {
    cli_abort("`stdcurve()` needs standards from at least 2 plates to estimate plate variation.")
  }
  list(
    detection_log = as.integer(family$detection == "log_concentration"),
    plates = plates,
    plate_i = plate_i,
    std_plate = match(std$plate, plates),
    std_x = log(std$conc),
    std_ct = std$ct
  )
}

# Plate index (into `plates`, the plates with standards) of each row of
# `data`, from its `plate_col` column. `what` names the data in messages.
.stdcurve_plate_index <- function(data, plate_col, plates, what = "data") {
  if (!plate_col %in% names(data)) {
    cli_abort("`{plate_col}` (the `stdcurve()` plate column) must be a column in `{what}`.")
  }
  plate <- as.character(data[[plate_col]])
  if (anyNA(plate)) {
    cli_abort("`{what}${plate_col}` must not contain missing values.")
  }
  unknown <- setdiff(unique(plate), plates)
  if (length(unknown)) {
    cli_abort(c("Every plate in `{what}` must have standards.",
      "x" = "No standards for plate{?s} {.val {unknown}}."))
  }
  match(plate, plates)
}

# Starting values from per-plate fits to the standards: Ct on log
# concentration (positives) and logistic detection. The plate parameters are
# added to `params`, and the intercept of the field linear predictor (`X`
# columns) starts where the mean curve puts the positive field Ct values `y`.
.stdcurve_start <- function(params, sc, y, X) {
  y <- as.matrix(y)[, 1L]
  if (any(y < 0, na.rm = TRUE)) {
    cli_abort("`stdcurve()` responses (Ct) must be >= 0, with 0 for non-detections.")
  }
  f <- function(x) if (sc$detection_log == 1L) x else exp(x)
  pos <- sc$std_ct > 0
  if (sum(pos) < 2L) {
    cli_abort("`stdcurve()` needs at least 2 positive (Ct > 0) standards.")
  }
  ct_fit <- stats::lm.fit(cbind(1, sc$std_x[pos]), sc$std_ct[pos])
  ct_coef <- unname(ct_fit$coefficients)
  per_plate <- vapply(seq_along(sc$plates), function(j) {
    k <- pos & sc$std_plate == j
    if (sum(k) < 2L || length(unique(sc$std_x[k])) < 2L) return(c(NA, NA))
    unname(stats::lm.fit(cbind(1, sc$std_x[k]), sc$std_ct[k])$coefficients)
  }, numeric(2))
  ct_sd <- apply(per_plate, 1L, stats::sd, na.rm = TRUE)
  bad <- !is.finite(ct_sd) | ct_sd <= 0
  ct_sd[bad] <- (0.1 * abs(ct_coef) + 0.01)[bad]
  det <- suppressWarnings(tryCatch(
    stats::glm.fit(cbind(1, f(sc$std_x)), as.numeric(pos),
      family = stats::binomial())$coefficients,
    error = function(e) c(0, 0)))
  det <- unname(det)
  det[!is.finite(det)] <- 0
  # Perfect separation is common (all high standards amplify); keep the
  # start finite and moderate.
  det <- pmin(pmax(det, -10), 10)
  # Matrix of plate coefficients (n_plate x 4), hyper means, and log SDs.
  mu <- c(ct_coef, det)
  params$std_xi <- matrix(mu, length(sc$plates), 4L, byrow = TRUE)
  params$std_mu <- mu
  params$ln_std_sigma <- log(c(ct_sd, 0.5, 0.5 * abs(det[2L]) + 0.01))
  # Invert the mean Ct curve for the latent log concentration
  y_pos <- !is.na(y) & y > 0
  intercept <- colnames(X) == "(Intercept)"
  if (any(y_pos) && any(intercept)) {
    params$b_j[intercept] <- mean((y[y_pos] - mu[[1]]) / mu[[2]])
  }
  params
}

# Plate coefficients for fitted rows `i`, as an n x 4 matrix.
rtmb_stdcurve_row_xi <- function(par, prepared, i) {
  par$std_xi[prepared$stdcurve$plate_i[i], , drop = FALSE]
}

# The plate curves at log concentration `x`, given each row's plate
# coefficients `xi`: the mean positive Ct and the logit detection probability.
# Plain arithmetic, so shared by the RTMB likelihood and the R-side methods.
.stdcurve_ct_mean <- function(xi, x) xi[, 1L] + xi[, 2L] * x

.stdcurve_logit_p <- function(xi, x, detection_log) {
  xi[, 3L] + xi[, 4L] * (if (detection_log == 1L) x else exp(x))
}

# Negative log density of the plate random effects and the standards, plus
# reports. `phi` is the natural-scale Ct SD.
rtmb_stdcurve_nll <- function(par, phi, prepared) {
  sc <- prepared$stdcurve
  xi <- par$std_xi
  std_sigma <- exp(par$ln_std_sigma)
  nll <- 0
  for (k in seq_len(4L)) {
    nll <- nll - sum(RTMB::dnorm(xi[, k], par$std_mu[[k]], std_sigma[[k]],
      log = TRUE))
  }
  xi_p <- xi[sc$std_plate, , drop = FALSE]
  present <- as.numeric(sc$std_ct > 0)
  logit_p <- .stdcurve_logit_p(xi_p, sc$std_x, sc$detection_log)
  ct_hat <- .stdcurve_ct_mean(xi_p, sc$std_x)
  nll <- nll - sum(RTMB::dbinom_robust(present, 1, logit_p, log = TRUE))
  # Positive standards only; data-dependent indexing, not an AD branch.
  pos <- which(present == 1)
  nll <- nll - sum(RTMB::dnorm(sc$std_ct[pos], ct_hat[pos], phi, log = TRUE))
  list(
    nll = nll,
    reports = list(std_xi = xi, std_mu = par$std_mu, std_sigma = std_sigma,
      std_ct_hat = ct_hat, std_logit_p = logit_p)
  )
}

# Field-sample observation family. With the log link, `mu` is DNA
# concentration; the curves are in log concentration, the linear predictor
# `s$eta`. `s$xi` holds each row's plate coefficients.
rtmb_obs_families$stdcurve <- list(
  logpdf = function(y, mu, s) {
    present <- as.numeric(y > 0)
    logit_p <- .stdcurve_logit_p(s$xi, s$eta, s$detection_log)
    ld <- RTMB::dbinom_robust(present, 1, logit_p, log = TRUE)
    # Non-detections have no Ct; the Gaussian term is finite and multiplied
    # by 0 for them, so no branching on AD values is needed.
    ld + present * RTMB::dnorm(y, .stdcurve_ct_mean(s$xi, s$eta), s$phi,
      log = TRUE)
  },
  mean = function(mu, s) {
    RTMB::plogis(.stdcurve_logit_p(s$xi, s$eta, s$detection_log)) *
      .stdcurve_ct_mean(s$xi, s$eta)
  },
  simulate = function(mu, s) {
    n <- length(mu)
    p <- stats::plogis(.stdcurve_logit_p(s$xi, s$eta, s$detection_log))
    stats::rbinom(n, 1, p) * stats::rnorm(n, .stdcurve_ct_mean(s$xi, s$eta),
      s$phi)
  }
)

# Expected Ct (detection probability times the mean positive Ct) of each
# fitted row, given log concentration `eta`, with plate coefficients from
# `theta` (by default at their modes). Used by `fitted()` and response
# residuals.
.stdcurve_expected_ct <- function(object, eta, theta = get_pars(object)) {
  sc <- object$tmb_data$stdcurve
  if (length(eta) != length(sc$plate_i)) {
    cli_abort("`stdcurve()` expected Ct values need one prediction per fitted row.")
  }
  xi <- theta$std_xi[sc$plate_i, , drop = FALSE]
  stats::plogis(.stdcurve_logit_p(xi, eta, sc$detection_log)) *
    .stdcurve_ct_mean(xi, eta)
}

#' Get eDNA standard curve estimates
#'
#' Extract the estimated PCR plate standard curves from a model fit with the
#' [stdcurve()] family.
#'
#' @param object A model fit with `family = stdcurve()`.
#' @param conf.level Confidence level for the intervals.
#'
#' @return A list with two data frames:
#' * `plates`: one row per plate and coefficient (`ct_intercept`, `ct_slope`,
#'   `det_intercept`, `det_slope`) with the estimate (conditional mode) and
#'   standard error.
#'
#' Parameters fixed with `sdmTMBcontrol(map = )` have `NA` standard errors.
#' * `hyper`: the mean (`mu`) and standard deviation (`sigma`) of each
#'   coefficient across plates, with confidence intervals.
#'
#' Positive Ct values have mean `ct_intercept + ct_slope * log(concentration)`;
#' the logit detection probability is
#' `det_intercept + det_slope * f(concentration)`, where `f()` is the identity
#' (`detection = "concentration"`) or `log()`
#' (`detection = "log_concentration"`).
#' @export
#' @examples
#' \dontrun{
#' fit <- sdmTMB(ct ~ 1, data = d, mesh = mesh,
#'   family = stdcurve(standards = standards))
#' get_stdcurve(fit)
#' }
get_stdcurve <- function(object, conf.level = 0.95) {
  assert_that(inherits(object, "sdmTMB"))
  sc <- object$tmb_data$stdcurve
  if (is.null(sc)) {
    cli_abort("`object` was not fit with the `stdcurve()` family.")
  }
  sdr <- object$sd_report
  if (is.null(sdr)) {
    cli_abort("Standard errors are unavailable; refit with `getsd = TRUE` in `sdmTMBcontrol()`.")
  }
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  n_p <- length(sc$plates)
  coefs <- .stdcurve_coef_names
  pars <- object$tmb_obj$env$parList(par = object$tmb_obj$env$last.par.best)
  # SEs from sdreport(); entries mapped off (fixed) get NA
  se <- function(name, what) {
    n <- length(pars[[name]])
    tab <- summary(sdr, what)
    vals <- unname(tab[rownames(tab) == name, 2L])
    m <- object$tmb_map[[name]]
    if (is.null(m)) return(if (length(vals) == n) vals else rep(NA_real_, n))
    lev <- as.integer(m)
    out <- rep(NA_real_, n)
    out[!is.na(lev)] <- vals[lev[!is.na(lev)]]
    out
  }
  plates <- data.frame(
    plate = rep(sc$plates, 4L),
    term = rep(coefs, each = n_p),
    estimate = c(pars$std_xi),
    std.error = se("std_xi", "random"),
    stringsAsFactors = FALSE
  )
  mu <- pars$std_mu
  mu_se <- se("std_mu", "fixed")
  ln_sigma <- pars$ln_std_sigma
  ln_sigma_se <- se("ln_std_sigma", "fixed")
  hyper <- data.frame(
    term = rep(coefs, 2L),
    parameter = rep(c("mu", "sigma"), each = 4L),
    estimate = c(mu, exp(ln_sigma)),
    std.error = c(mu_se, exp(ln_sigma) * ln_sigma_se),
    conf.low = c(mu - z * mu_se, exp(ln_sigma - z * ln_sigma_se)),
    conf.high = c(mu + z * mu_se, exp(ln_sigma + z * ln_sigma_se)),
    stringsAsFactors = FALSE
  )
  rownames(hyper) <- NULL
  list(plates = plates, hyper = hyper)
}
