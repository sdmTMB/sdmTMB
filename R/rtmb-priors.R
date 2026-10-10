# Priors by their `sdmTMBpriors()` names. `fit.R` passes most priors to C++
# as one unnamed vector, so rebuild the names from the same template. Each
# prior is NULL unless all of its hyperparameters are set.
rtmb_prior_inputs <- function(data) {
  template <- sdmTMBpriors()
  template$b <- template$sigma_V <- NULL
  template$custom <- template$custom_log_jacobian <- NULL
  sizes <- lengths(template)
  # fits saved before `matern_svc` lack its (trailing) values
  values <- c(data$priors, rep(NA, sum(sizes) - length(data$priors)))
  stopifnot(sum(sizes) == length(values))
  values <- split(values,
    factor(rep(names(sizes), sizes), levels = names(sizes)))
  priors <- lapply(values, function(x) if (anyNA(x)) NULL else unname(x))
  if (data$priors_b_n > 0L) {
    priors$b <- list(index = data$priors_b_index + 1L,
      mean = data$priors_b_mean, Sigma = data$priors_b_Sigma)
  }
  priors$sigma_V <- data$priors_sigma_V
  priors$stan <- data$stan_flag == 1L
  priors$custom <- data$priors_custom
  priors
}

# PC prior for a Matern field; `prior` is `pc_matern()`'s
# c(range_gt, sigma_lt, range_prob, sigma_prob).
rtmb_pc_matern <- function(log_tau, log_kappa, prior, include_sigma = TRUE,
                           include_range = TRUE, stan = FALSE) {
  lambda_range <- -log(prior[[3L]]) * prior[[1L]]
  lambda_sigma <- -log(prior[[4L]]) / prior[[2L]]
  range <- sqrt(8) / exp(log_kappa)
  sigma <- exp(-log_tau - log_kappa) / sqrt(4 * pi)
  # The sigma part applies to estimated fields and the range part once per
  # shared range. Each part's Jacobian term, from (log_tau, log_kappa) to
  # (sigma, range), is log(sigma) or log(range).
  log_density <- 0
  if (include_sigma) {
    log_density <- log_density + log(lambda_sigma) - lambda_sigma * sigma
    if (stan) log_density <- log_density + log(sigma)
  }
  if (include_range) {
    log_density <- log_density + log(lambda_range) -
      2 * log(range) - lambda_range / range
    if (stan) log_density <- log_density + log(range)
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
        prior$matern_s, include_sigma = prepared$sigma_prior[1L, m],
        include_range = prepared$range_prior[1L, m], stan = stan)
    }
    if (!is.null(prior$matern_st)) {
      nll <- nll - rtmb_pc_matern(par$ln_tau_E[[m]], par$ln_kappa[2L, m],
        prior$matern_st, include_sigma = prepared$sigma_prior[2L, m],
        include_range = prepared$range_prior[2L, m], stan = stan)
    }
    if (!is.null(prior$matern_svc)) {
      for (z in seq_len(nrow(prepared$svc_kappa_row))) {
        nll <- nll - rtmb_pc_matern(par$ln_tau_Z[z, m],
          par$ln_kappa[prepared$svc_kappa_row[z, m], m], prior$matern_svc,
          include_sigma = prepared$sigma_prior[2L + z, m],
          include_range = prepared$range_prior[2L + z, m], stan = stan)
      }
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
