# Latent effect densities for the RTMB objective. Each helper evaluates
# component `m` and returns its negative log density and the values after any
# requested simulation, so predictors are always computed from the simulated
# draw.

rtmb_value <- function(x) if (inherits(x, "simref")) x$orig else x

# One GMRF slice `x`: drawn first when simulating, then its negative log
# density.
rtmb_gmrf <- function(x, Q, scale, simulate) {
  if (simulate) {
    draw <- RTMB::simref(length(x))
    RTMB::dgmrf(draw, Q = Q, scale = scale, log = TRUE)
    x <- draw$value
  }
  if (inherits(Q, "dgeMatrix")) {
    # Infinite parameter draws can turn a sparse precision into a dense
    # matrix of NaNs through 0 * Inf. Its density is undefined, but reports
    # from the same parameter draw can still be evaluated.
    if (any(!is.finite(Q@x))) return(list(nll = Inf, value = x))
    Q <- methods::as(Q, "CsparseMatrix")
  }
  list(nll = -sum(RTMB::dgmrf(x, Q = Q, scale = scale, log = TRUE)),
    value = x)
}

# Log spatiotemporal SD by time step: constant, or log-linear with a slope.
rtmb_log_sigma_E <- function(par, prepared, m) {
  value <- rep(rtmb_log_field_sd(par$ln_tau_E[[m]], par$ln_kappa[2L, m],
    prepared$precision), prepared$n_t)
  if (prepared$epsilon_trend) {
    value <- value + par$b_epsilon[[m]] * prepared$epsilon_predictor
  }
  value
}

# Spatiotemporal field `epsilon_st` with IID, AR1, or RW time structure.
# Simulation replaces only `simulate_t` steps and conditions each draw on the
# preceding retained or simulated field.
rtmb_spatiotemporal_field <- function(epsilon_st, par, theta, prepared, Q,
                                      log_sigma_E, simulate, m) {
  inputs <- prepared$precision
  ar1 <- prepared$epsilon_ar1[[m]]
  rw <- prepared$epsilon_rw[[m]]
  rho <- theta$rho[[m]]
  nll <- 0
  ln_kappa <- par$ln_kappa[2L, m]
  scale <- rtmb_gmrf_scale(log_sigma_E, ln_kappa, inputs)
  previous_mean <- function(t) {
    if (t == 1L) return(0)
    if (ar1) rho * epsilon_st[, t - 1L, m]
    else if (rw) epsilon_st[, t - 1L, m]
    else 0
  }
  innovation_scale <- if (ar1) sqrt(1 - rho^2) else 1
  for (t in seq_len(prepared$n_t)) {
    step <- if (t > 1L) innovation_scale else 1
    draw <- simulate && t %in% prepared$simulate_t
    innovation <- rtmb_gmrf((epsilon_st[, t, m] - previous_mean(t)) / step,
      Q, scale[[t]], draw)
    if (draw) epsilon_st[, t, m] <- previous_mean(t) + step * innovation$value
    nll <- nll + innovation$nll
  }
  if (ar1) {
    nll <- nll + (prepared$n_t - 1L) * dim(epsilon_st)[1L] *
      log(innovation_scale)
  }
  list(nll = nll, value = epsilon_st)
}

# Time-varying coefficients `b_rw_t` as RW (flat first step), RW0, or AR1.
rtmb_time_varying <- function(b_rw_t, theta, prepared, simulate, m) {
  ar1 <- prepared$time_varying_type == "ar1"
  nll <- 0
  for (k in seq_len(dim(b_rw_t)[2L])) {
    sigma <- theta$sigma_V[k, m]
    rho <- if (ar1) theta$rho_time[k, m] else 0
    center <- function(t) {
      if (t == 1L) 0 else if (ar1) rho * b_rw_t[t - 1L, k, m]
      else b_rw_t[t - 1L, k, m]
    }
    step_sd <- function(t) {
      if (t > 1L && ar1) sigma * sqrt(1 - rho^2) else sigma
    }
    steps <- seq_len(prepared$n_t)
    if (prepared$time_varying_type == "rw") steps <- steps[-1L]
    if (simulate) {
      for (t in intersect(steps, prepared$simulate_t)) {
        b_rw_t[t, k, m] <- stats::rnorm(1L, center(t), step_sd(t))
      }
    }
    for (t in steps) {
      nll <- nll - RTMB::dnorm(b_rw_t[t, k, m], center(t), step_sd(t),
        log = TRUE)
    }
  }
  list(nll = nll, value = b_rw_t)
}

# IID random intercepts and slopes, correlated within a group when requested.
rtmb_iid_effects <- function(re_b_pars, par, groups, simulate, m) {
  nll <- 0
  for (group in groups) {
    sds <- exp(par$re_cov_pars[group$sd_rows, m])
    if (group$dimension == 1L) {
      index <- unlist(group$coefficients)
      if (simulate) re_b_pars[index, m] <- stats::rnorm(length(index), 0, sds)
      nll <- nll - sum(RTMB::dnorm(re_b_pars[index, m], 0, sds, log = TRUE))
    } else {
      correlation <- RTMB::unstructured(group$dimension)$corr(
        par$re_cov_pars[group$corr_rows, m])
      covariance <- correlation * (sds %o% sds)
      if (simulate) covariance_chol <- t(chol(covariance))
      for (index in group$coefficients) {
        if (simulate) {
          re_b_pars[index, m] <- as.vector(covariance_chol %*%
            stats::rnorm(group$dimension))
        }
        nll <- nll - sum(RTMB::dmvnorm(re_b_pars[index, m],
          Sigma = covariance, log = TRUE))
      }
    }
  }
  list(nll = nll, value = re_b_pars)
}

# Penalized smoother coefficients `b_smooth`.
rtmb_smooth_effects <- function(b_smooth, par, smooth_index, simulate, m) {
  nll <- 0
  for (s in seq_along(smooth_index)) {
    index <- smooth_index[[s]]
    sd <- exp(par$ln_smooth_sigma[s, m])
    if (simulate) b_smooth[index, m] <- stats::rnorm(length(index), 0, sd)
    nll <- nll - sum(RTMB::dnorm(b_smooth[index, m], 0, sd, log = TRUE))
  }
  list(nll = nll, value = b_smooth)
}

# Evaluate or simulate all latent effects for every component. Returns the
# effect arrays under their C++ parameter names and shapes, `log_sigma_E`
# (time by component), and `nll`, their summed negative log density.
# `simulating` names the parameters passed as simulation references; of
# those, only effects requested in `prepared$simulate_re` are drawn.
rtmb_latent_effects <- function(par, theta, prepared, simulating) {
  "[<-" <- RTMB::ADoverload("[<-")
  inputs <- prepared$precision
  n_m <- prepared$n_m
  # The preferential-sampling field `xi_s` has no simulation option yet, so
  # it is never drawn.
  simulate <- function(name) {
    name %in% simulating && isTRUE(prepared$simulate_re[name])
  }
  effects <- c(list(nll = 0), par[c("omega_s", "epsilon_st", "zeta_s",
    "b_rw_t", "re_b_pars", "b_smooth")],
    list(log_sigma_E = matrix(0, prepared$n_t, n_m)))
  add <- function(result, name) {
    effects$nll <<- effects$nll + result$nll
    effects[[name]] <<- result$value
  }
  # Add a GMRF slice's density and return the slice, drawn if simulating.
  gmrf <- function(x, Q, scale, name) {
    result <- rtmb_gmrf(x, Q, scale, simulate(name))
    effects$nll <<- effects$nll + result$nll
    result$value
  }
  for (m in seq_len(n_m)) {
    if (prepared$spatial[[m]] || prepared$temporal[[m]] ||
        prepared$svc_density[[m]]) {
      Q <- rtmb_precision(inputs, theta, 1L, m)
    }
    if (prepared$spatial[[m]]) {
      scale <- rtmb_gmrf_scale(theta$log_sigma_O[1L, m], par$ln_kappa[1L, m],
        inputs)
      effects$omega_s[, m] <- gmrf(effects$omega_s[, m], Q, scale, "omega_s")
    }
    if (prepared$svc_density[[m]]) {
      for (z in seq_len(dim(effects$zeta_s)[2L])) {
        scale <- rtmb_gmrf_scale(theta$log_sigma_Z[z, m], par$ln_kappa[1L, m],
          inputs)
        effects$zeta_s[, z, m] <- gmrf(effects$zeta_s[, z, m], Q, scale,
          "zeta_s")
      }
    }
    log_sigma_E <- rtmb_log_sigma_E(par, prepared, m)
    effects$log_sigma_E[, m] <- log_sigma_E
    if (prepared$temporal[[m]]) {
      shared <- prepared$share_range[[m]] || rtmb_areal(inputs)
      Q_st <- if (shared) Q else rtmb_precision(inputs, theta, 2L, m)
      add(rtmb_spatiotemporal_field(effects$epsilon_st, par, theta, prepared,
        Q_st, log_sigma_E, simulate("epsilon_st"), m), "epsilon_st")
    }
    if (prepared$iid_re[[m]]) {
      add(rtmb_iid_effects(effects$re_b_pars, par, prepared$re_groups[[m]],
        simulate("re_b_pars"), m), "re_b_pars")
    }
    if (prepared$time_varying) {
      add(rtmb_time_varying(effects$b_rw_t, theta, prepared,
        simulate("b_rw_t"), m), "b_rw_t")
    }
    if (prepared$smooths) {
      add(rtmb_smooth_effects(effects$b_smooth, par, prepared$smooth_index,
        simulate("b_smooth"), m), "b_smooth")
    }
  }
  # Sampling-only field of a preferential-sampling model: time invariant,
  # isotropic, and independent of the catch fields.
  if (!is.null(prepared$preferential) && prepared$preferential$xi) {
    inputs <- prepared$preferential$precision
    kappa <- theta$xi$kappa
    dim(kappa) <- c(1L, 1L)
    Q <- rtmb_precision(inputs, list(kappa = kappa), 1L, 1L)
    scale <- rtmb_gmrf_scale(theta$xi$log_sigma, par$ln_kappa_xi, inputs)
    effects$xi_s <- gmrf(par$xi_s, Q, scale, "xi_s")
  }
  effects
}
