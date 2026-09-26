# The objective as a closure over `prepared` alone, built in its own frame so
# the AD tape doesn't also capture `data` and `parameters` from
# make_sdmTMB_adfun()'s frame.
rtmb_make_objective <- function(prepared) {
  force(prepared)
  function(par) rtmb_objective(par, prepared)
}

# RTMB objective: resolve parameters to their natural scale, evaluate latent
# effects, compute predictors from them, evaluate observations and priors,
# then register reports. `par` keeps the C++ parameter names and shapes.
rtmb_objective <- function(par, prepared) {
  # During simulate(), random parameters arrive as simulation references.
  # Record which, then work with their current values throughout.
  simulating <- names(par)[vapply(par, inherits, NA, "simref")]
  par <- lapply(par, rtmb_value)
  theta <- rtmb_transform(par, prepared)
  effects <- rtmb_latent_effects(par, theta, prepared, simulating)
  fitted <- rtmb_linear_predictors(par, theta, effects, prepared, prepared$fit)
  # Inactive family components have no predictor, as in the C++ template.
  active <- prepared$fit$active
  fitted$eta <- fitted$eta * active
  fitted$rw <- fitted$rw * active
  fitted$epsilon <- fitted$epsilon * active
  obs <- rtmb_observations(par, theta, prepared, fitted$eta)
  jnll <- effects$nll + sum(obs$jnll_obs) + rtmb_prior_nll(par, theta, prepared)
  sampling <- NULL
  if (!is.null(prepared$preferential)) {
    sampling <- rtmb_sampling(par, theta, effects, prepared)
    jnll <- jnll + sampling$nll
  }
  projected <- derived <- NULL
  if (!is.null(prepared$proj)) {
    projected <- rtmb_linear_predictors(par, theta, effects, prepared,
      prepared$proj)
    if (prepared$mixture) {
      projected$eta <- rtmb_mixture_eta(projected$eta, theta, prepared)
    }
    if (prepared$n_m > 1L) {
      projected$combined <- rtmb_combined_projection(projected, theta, prepared)
    }
    if (any(prepared$derived)) {
      derived <- rtmb_derived_indices(par, theta, prepared, projected)
      jnll <- jnll + derived$nll
    }
  }
  rtmb_report(par, theta, prepared, effects, fitted, obs, projected, derived,
    sampling, simulating)
  jnll
}

# Natural-scale parameters, resolved once per evaluation so that fields,
# observations, priors, and reports share one definition of each transform.
rtmb_transform <- function(par, prepared) {
  "[<-" <- RTMB::ADoverload("[<-")
  inputs <- prepared$precision
  n_m <- prepared$n_m
  # Maps the real line to a correlation in (-1, 1).
  correlation <- function(x) 2 * RTMB::plogis(x) - 1
  # Components without a spatial field report sigma_O = log_sigma_O = 0.
  log_sigma_O <- sigma_O <- matrix(0, 1L, n_m)
  log_sigma_Z <- sigma_Z <- matrix(0, nrow(par$ln_tau_Z), n_m)
  for (m in seq_len(n_m)) {
    ln_kappa <- par$ln_kappa[1L, m]
    if (prepared$include_spatial[[m]]) {
      log_sigma_O[1L, m] <- rtmb_log_field_sd(par$ln_tau_O[[m]], ln_kappa,
        inputs)
      sigma_O[1L, m] <- exp(log_sigma_O[1L, m])
    }
    if (prepared$svc) {
      log_sigma_Z[, m] <- rtmb_log_field_sd(par$ln_tau_Z[, m], ln_kappa,
        inputs)
      sigma_Z[, m] <- exp(log_sigma_Z[, m])
    }
  }
  # As in C++, sigma_V is zero without time-varying coefficients.
  sigma_V <- if (prepared$time_varying) exp(par$ln_tau_V) else
    matrix(0, nrow(par$ln_tau_V), n_m)
  rho_time <- correlation(par$rho_time_unscaled)
  dim(rho_time) <- dim(par$rho_time_unscaled)
  b <- par$b_threshold
  kappa <- exp(par$ln_kappa)
  xi <- !is.null(prepared$preferential) && prepared$preferential$xi
  list(
    # Random fields
    kappa = kappa,
    range = sqrt(8) / kappa,
    sigma_O = sigma_O,
    log_sigma_O = log_sigma_O,
    sigma_Z = sigma_Z,
    log_sigma_Z = log_sigma_Z,
    rho = correlation(par$ar1_phi),
    rho_sar = correlation(par$logit_rho_sar),
    alpha_car = RTMB::plogis(par$logit_rho_sar),
    H = if (prepared$anisotropy) {
      lapply(seq_len(n_m), function(m) rtmb_aniso_H(par$ln_H_input[, m]))
    },
    # Preferential-sampling field
    xi = if (xi) {
      log_sigma <- rtmb_log_field_sd(par$ln_tau_xi, par$ln_kappa_xi,
        prepared$preferential$precision)
      list(kappa = exp(par$ln_kappa_xi), log_sigma = log_sigma,
        sigma = exp(log_sigma), range = sqrt(8) / exp(par$ln_kappa_xi))
    },
    # Preferential-sampling temporal process SDs
    sigma_b_pref = if (!is.null(par$ln_sigma_b_pref)) exp(par$ln_sigma_b_pref),
    sigma_alpha_pref = if (!is.null(par$ln_sigma_alpha_pref)) {
      exp(par$ln_sigma_alpha_pref)
    },
    # Time-varying coefficients
    sigma_V = sigma_V,
    rho_time = rho_time,
    # Threshold covariate: breakpoint slope and cut, or logistic s50, s95,
    # and maximum, with s95 > s50.
    threshold = switch(prepared$threshold,
      none = list(),
      breakpt = list(s_slope = b[1L, ], s_cut = b[2L, ]),
      logistic = list(s50 = b[1L, ], s95 = b[1L, ] + exp(b[2L, ]),
        s_max = b[3L, ])),
    # Observation families
    phi = exp(par$ln_phi),
    tweedie_p = RTMB::plogis(par$thetaf) + 1,
    student_df = exp(par$ln_student_df) + 1,
    # The larger mixture component's mean is `mix_ratio` times the smaller.
    p_extreme = RTMB::plogis(par$logit_p_extreme),
    mix_ratio = exp(par$log_ratio_mix) + 1
  )
}
