# Priors by their `sdmTMBpriors()` names. `fit.R` passes most priors to C++
# as one unnamed vector, so rebuild the names from the same template. Each
# prior is NULL unless all of its hyperparameters are set.
rtmb_prior_inputs <- function(data) {
  template <- sdmTMBpriors()
  template$b <- template$sigma_V <- NULL
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
  priors
}

# PC prior for a Matern field; `prior` is `pc_matern()`'s
# c(range_gt, sigma_lt, range_prob, sigma_prob).
rtmb_pc_matern <- function(log_tau, log_kappa, prior, share_range = FALSE,
                           stan = FALSE) {
  lambda_range <- -log(prior[[3L]]) * prior[[1L]]
  lambda_sigma <- -log(prior[[4L]]) / prior[[2L]]
  range <- sqrt(8) / exp(log_kappa)
  sigma <- exp(-log_tau - log_kappa) / sqrt(4 * pi)
  log_density <- log(lambda_sigma) - lambda_sigma * sigma
  if (!share_range) {
    log_density <- log_density + log(lambda_range) -
      2 * log(range) - lambda_range / range
  }
  if (stan) {
    log_density <- log_density + log(sqrt(8)) - 2 * log(range) +
      log(sqrt(4 * pi)) + log_kappa
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
  for (m in seq_len(prepared$n_m)) {
    if (!is.null(prior$b)) {
      b <- if (m == 1L) par$b_j else par$b_j2
      nll <- nll - sum(RTMB::dmvnorm(b[prior$b$index] - prior$b$mean,
        Sigma = prior$b$Sigma, log = TRUE))
    }
    if (!is.null(prior$matern_s)) {
      nll <- nll - rtmb_pc_matern(par$ln_tau_O[[m]], par$ln_kappa[1L, m],
        prior$matern_s, stan = stan)
    }
    if (!is.null(prior$matern_st)) {
      nll <- nll - rtmb_pc_matern(par$ln_tau_E[[m]], par$ln_kappa[2L, m],
        prior$matern_st, share_range = prepared$share_range[[m]], stan = stan)
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
    for (k in seq_len(nrow(prior$sigma_V))) {
      if (!anyNA(prior$sigma_V[k, ])) {
        sigma_V <- theta$sigma_V[k, m]
        nll <- nll - RTMB::dgamma(sigma_V, shape = prior$sigma_V[k, 1L],
          scale = prior$sigma_V[k, 2L], log = TRUE)
        if (stan) nll <- nll - log(sigma_V)
      }
    }
  }
  nll <- nll + normal(theta$phi, prior$phi)
  if (stan && !is.null(prior$phi)) nll <- nll - sum(par$ln_phi)
  nll
}
