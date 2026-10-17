# Priors by their `sdmTMBpriors()` names. `fit.R` passes most priors to C++
# as one unnamed vector, so rebuild the names from the same template. Each
# prior is NULL unless all of its hyperparameters are set.
rtmb_prior_inputs <- function(data) {
  template <- sdmTMBpriors()
  template$b <- template$sigma_V <- NULL
  template$custom <- template$custom_log_jacobian <- template$aniso <- NULL
  sizes <- lengths(template)
  stopifnot(sum(sizes) == length(data$priors))
  values <- split(data$priors,
    factor(rep(names(sizes), sizes), levels = names(sizes)))
  priors <- lapply(values, function(x) if (anyNA(x)) NULL else unname(x))
  if (data$priors_b_n > 0L) {
    priors$b <- list(index = data$priors_b_index + 1L,
      mean = data$priors_b_mean, Sigma = data$priors_b_Sigma)
  }
  priors$sigma_V <- data$priors_sigma_V
  priors$stan <- data$stan_flag == 1L
  priors$custom <- data$priors_custom
  priors$aniso <- rtmb_aniso_prior_inputs(data$priors_aniso, priors)
  priors
}

# PC prior for a Matern field; `prior` is `pc_matern()`'s
# c(range_gt, sigma_lt, range_prob, sigma_prob).
# `lambda_range` replaces the range rate under the anisotropy prior.
rtmb_pc_matern <- function(log_tau, log_kappa, prior, share_range = FALSE,
                           stan = FALSE,
                           lambda_range = -log(prior[[3L]]) * prior[[1L]]) {
  lambda_sigma <- -log(prior[[4L]]) / prior[[2L]]
  range <- sqrt(8) / exp(log_kappa)
  sigma <- exp(-log_tau - log_kappa) / sqrt(4 * pi)
  log_density <- log(lambda_sigma) - lambda_sigma * sigma
  if (!share_range) {
    log_density <- log_density + log(lambda_range) -
      2 * log(range) - lambda_range / range
  }
  # Jacobian from (log_tau, log_kappa) to (range, sigma) is range * sigma;
  # a shared range is counted in the spatial prior only.
  if (stan) {
    log_density <- log_density + log(sigma)
    if (!share_range) log_density <- log_density + log(range)
  }
  log_density
}

# Priors apply to every component with the same hyperparameters, as in the
# C++ template. With `stan`, Jacobians put priors on the transformed scale.
rtmb_prior_nll <- function(par, theta, prepared) {
  prior <- prepared$priors
  stan <- prior$stan
  # Normal prior `p = c(mean, sd)`; nothing when the prior or value is absent.
  normal <- function(x, p) {
    if (is.null(p) || !length(x)) return(0)
    -sum(RTMB::dnorm(x, p[[1L]], p[[2L]], log = TRUE))
  }
  threshold <- theta$threshold
  nll <- 0
  aniso <- prior$aniso
  for (m in seq_len(prepared$n_m)) {
    # Under the anisotropy prior, the range rate scales with f(|v|)
    if (!is.null(aniso)) {
      cosh_v <- rtmb_aniso_cosh(par$ln_H_input[, m])
      f_v <- pc_aniso_f(cosh_v)
      if (m %in% aniso$components) {
        nll <- nll - rtmb_pc_aniso(cosh_v, aniso$lambda_v)
      }
    }
    if (!is.null(prior$b)) {
      b <- if (m == 1L) par$b_j else par$b_j2
      nll <- nll - sum(RTMB::dmvnorm(b[prior$b$index] - prior$b$mean,
        Sigma = prior$b$Sigma, log = TRUE))
    }
    if (!is.null(prior$matern_s)) {
      nll <- nll - rtmb_pc_matern(par$ln_tau_O[[m]], par$ln_kappa[1L, m],
        prior$matern_s, stan = stan,
        lambda_range = if (is.null(aniso)) -log(prior$matern_s[[3L]]) *
          prior$matern_s[[1L]] else aniso$range_s * f_v)
    }
    if (!is.null(prior$matern_st)) {
      nll <- nll - rtmb_pc_matern(par$ln_tau_E[[m]], par$ln_kappa[2L, m],
        prior$matern_st, stan = stan,
        lambda_range = if (is.null(aniso)) -log(prior$matern_st[[3L]]) *
          prior$matern_st[[1L]] else aniso$range_st * f_v,
        # the range is shared and the spatial prior already includes it
        share_range = prepared$share_range[[m]] && !is.null(prior$matern_s))
    }
    nll <- nll + normal(theta$rho[m], prior$ar1_rho)
    if (stan && !is.null(prior$ar1_rho)) {
      raw <- par$ar1_phi[[m]]
      nll <- nll - (log(2) + raw - 2 * RTMB::logspace_add(0, raw))
    }
    nll <- nll +
      normal(threshold$s_slope[m], prior$threshold_breakpt_slope) +
      normal(threshold$s_cut[m], prior$threshold_breakpt_cut) +
      normal(threshold$s50[m], prior$threshold_logistic_s50) +
      normal(threshold$s95[m], prior$threshold_logistic_s95) +
      normal(threshold$s_max[m], prior$threshold_logistic_smax)
    # Jacobian for s95 = s50 + exp(b_threshold[2, ])
    if (stan && length(threshold$s95) &&
      !is.null(prior$threshold_logistic_s95)) {
      nll <- nll - par$b_threshold[2L, m]
    }
    for (k in seq_len(nrow(prior$sigma_V))) {
      if (!anyNA(prior$sigma_V[k, ])) {
        sigma_V <- theta$sigma_V[k, m]
        # column 3: 0 = gamma, 1 = lognormal (absent in older fits)
        nll <- nll - if (ncol(prior$sigma_V) > 2L && prior$sigma_V[k, 3L] == 1) {
          RTMB::dnorm(log(sigma_V), prior$sigma_V[k, 1L],
            prior$sigma_V[k, 2L], log = TRUE) - log(sigma_V)
        } else {
          RTMB::dgamma(sigma_V, shape = prior$sigma_V[k, 1L],
            scale = prior$sigma_V[k, 2L], log = TRUE)
        }
        if (stan) nll <- nll - log(sigma_V)
      }
    }
  }
  nll <- nll + normal(theta$phi, prior$phi)
  if (stan && !is.null(prior$phi)) nll <- nll - sum(par$ln_phi)
  nll
}

# Anisotropy PC prior ------------------------------------------------------

# Check the `pc_anisotropy()` prior and record which components' `H` it
# applies to: each distinct estimated column of `ln_H_input`, so a shared `H`
# is counted once. NULL without the prior.
aniso_prior_spec <- function(prior, backend, anisotropy, map, n_m) {
  if (is.null(prior) || is.na(prior[[1L]])) return(NULL)
  if (backend != "rtmb") {
    cli_abort("`pc_anisotropy()` needs `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  if (!isTRUE(anisotropy == 1L)) {
    cli_abort("`pc_anisotropy()` needs `anisotropy = TRUE` and a mesh without barriers.")
  }
  components <- seq_len(n_m)
  if (!is.null(map)) {
    cols <- t(matrix(as.integer(map), nrow = 2L))
    components <- which(!duplicated(cols) & rowSums(!is.na(cols)) > 0L)
  }
  list(ratio_lt = prior[[1L]], ratio_prob = prior[[2L]],
    components = components)
}

# Rates for the anisotropy prior from `aniso_prior_spec()` and the
# `pc_matern()` priors, following Theorem 5 of Llamazares-Elias et al. (2024):
# P(major / minor range ratio > ratio_lt) = ratio_prob and, marginally,
# P(range < range_gt) = range_prob. `range_s` and `range_st` times f(|v|) are
# the rates of the conditional priors on the ranges given the anisotropy.
rtmb_aniso_prior_inputs <- function(spec, priors) {
  if (is.null(spec)) return(NULL)
  f0 <- pc_aniso_f(1)
  lambda_v <- -log(spec$ratio_prob) /
    (pc_aniso_f((spec$ratio_lt + 1 / spec$ratio_lt) / 2) - f0)
  range_rate <- function(p) {
    if (is.null(p)) return(NULL)
    y <- lambda_v * f0
    w <- lambert_w0_log(y + log(y) - log(p[[3L]]))
    p[[1L]] * (w / f0 - lambda_v)
  }
  c(spec, list(lambda_v = lambda_v,
    range_s = range_rate(priors$matern_s),
    range_st = range_rate(priors$matern_st)))
}

# Principal branch of the Lambert W function of exp(log_z), by Newton's
# method on log(W), which is convex in log(W) and safe for large z.
lambert_w0_log <- function(log_z) {
  u <- if (log_z > 1) log(log_z) else 0
  for (i in seq_len(100L)) {
    step <- (u + exp(u) - log_z) / (1 + exp(u))
    u <- u - step
    if (abs(step) < 1e-12) break
  }
  exp(u)
}

# cosh(|v|) for `H = MakeH(x)`, which is tr(H) / 2 since det(H) = 1 and the
# eigenvalues of H are exp(+/- |v|).
rtmb_aniso_cosh <- function(x) {
  (exp(x[[1L]]) + (1 + x[[2L]]^2) * exp(-x[[1L]])) / 2
}

# The distance function f(|v|) = sqrt(pi / 3 * (3 * cosh(2|v|) + 1)) in
# terms of cosh(|v|).
pc_aniso_f <- function(cosh_v) sqrt(pi * (6 * cosh_v^2 - 2) / 3)

# Log PC prior density of the anisotropy as a density of `ln_H_input`: the
# density of `v`, lambda_v * f'(|v|) / (2 * pi * |v|) *
# exp(-lambda_v * (f(|v|) - f(0))), times the Jacobian |v| / sinh(|v|) from
# `ln_H_input` to `v`. It is smooth at isotropy and has its mode there.
rtmb_pc_aniso <- function(cosh_v, lambda_v) {
  log(lambda_v) + 0.5 * log(4 * pi) + log(cosh_v) -
    0.5 * log(2 * cosh_v^2 - 2 / 3) - log(2 * pi) -
    lambda_v * (pc_aniso_f(cosh_v) - pc_aniso_f(1))
}

# Custom priors ------------------------------------------------------------

# The `custom` and `custom_log_jacobian` functions from `sdmTMBpriors()` as
# they are stored in the data list, or NULL without them.
custom_prior_spec <- function(priors, backend) {
  if (is.null(priors$custom)) return(NULL)
  if (backend != "rtmb") {
    cli_abort("Custom priors need `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  list(density = priors$custom, log_jacobian = priors$custom_log_jacobian)
}

# Log density terms from one custom prior function, naming the function if
# it fails.
rtmb_custom_terms <- function(par, theta, custom, which) {
  arg <- c(density = "custom", log_jacobian = "custom_log_jacobian")[[which]]
  f <- custom[[which]]
  if (is.null(f)) return(NULL)
  tryCatch(f(par, theta), error = function(e) {
    cli_abort("The {.arg {arg}} prior function failed.", parent = e,
      call = NULL)
  })
}

# Summed log density of the custom priors, with the Jacobian only for
# `bayesian = TRUE`. Zero without custom priors.
rtmb_custom_log_density <- function(par, theta, prepared) {
  custom <- prepared$priors$custom
  if (is.null(custom)) return(0)
  out <- sum(rtmb_custom_terms(par, theta, custom, "density"))
  if (prepared$priors$stan) {
    out <- out + sum(rtmb_custom_terms(par, theta, custom, "log_jacobian"))
  }
  out
}

# Check that the custom prior functions return finite numeric values at the
# starting `parameters`, evaluated with plain numbers before taping. The
# Jacobian is only checked if it enters the objective (`bayesian = TRUE`).
rtmb_check_custom_priors <- function(parameters, prepared) {
  custom <- prepared$priors$custom
  if (is.null(custom)) return(invisible())
  theta <- rtmb_transform(parameters, prepared)
  used <- c("density", if (prepared$priors$stan) "log_jacobian")
  for (which in used) {
    if (is.null(custom[[which]])) next
    arg <- c(density = "custom", log_jacobian = "custom_log_jacobian")[[which]]
    x <- rtmb_custom_terms(parameters, theta, custom, which)
    if (!is.numeric(x) || !length(x)) {
      cli_abort("The {.arg {arg}} prior function must return numeric log densities.")
    }
    if (!all(is.finite(x))) {
      cli_abort("The {.arg {arg}} prior function returned non-finite values at the starting parameters.")
    }
  }
  invisible()
}
