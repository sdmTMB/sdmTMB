rtmb_inverse_link <- function(eta, link) {
  switch(link,
    identity = eta,
    log = exp(eta),
    logit = RTMB::plogis(eta),
    inverse = 1 / eta,
    cloglog = 1 - exp(-exp(eta)))
}

rtmb_link <- function(mu, link) {
  switch(link,
    identity = mu,
    log = log(mu),
    logit = RTMB::qlogis(mu),
    inverse = 1 / mu,
    cli::cli_abort("Link not implemented."))
}

# Logit of the inverse link without losing accuracy, as used by binomial
# likelihoods. log(exp(exp(eta)) - 1) is the cloglog case.
rtmb_logit_inverse_link <- function(eta, link) {
  switch(link,
    logit = eta,
    cloglog = RTMB::logspace_sub(exp(eta), 0),
    RTMB::qlogis(rtmb_inverse_link(eta, link)))
}

# Response mean of one component, as used for combined projections. The
# truncated negative binomials report the mean of the truncated distribution.
rtmb_component_mean <- function(eta, family, m, theta) {
  mu <- rtmb_inverse_link(eta, family$link[[m]])
  spec <- rtmb_obs_family(family$family[[m]])
  if (!is.null(spec$log_nzprob)) {
    mu <- mu / exp(spec$log_nzprob(mu, theta$phi[[family$phi]]))
  }
  mu
}

# Pick rows `i` from a row-aligned vector; scalars apply to every row.
rtmb_pick <- function(x, i) if (length(x) == 1L) x else x[i]

# Parameters and means of component `m` for fitted rows `i` of one family.
# For families with `logit_mu` (binomial), `mu` is on the logit scale.
rtmb_obs_state <- function(i, m, family, eta, par, theta, prepared, ln_phi_i) {
  fit <- prepared$fit
  name <- family$family[[m]]
  spec <- rtmb_obs_family(name)
  state <- list(eta = eta[i, m], ln_phi = 0, link = family$link[[m]],
    size = fit$size[i])
  if (!is.na(family$phi)) state$ln_phi <- par$ln_phi[[family$phi]]
  if (prepared$dispersion_model && m == prepared$n_m) {
    state$ln_phi <- ln_phi_i[i]
  }
  state$phi <- exp(state$ln_phi)
  if (!is.na(family$thetaf)) {
    state$tweedie_p <- theta$tweedie_p[[family$thetaf]]
  }
  if (!is.na(family$student_df)) {
    state$df <- theta$student_df[[family$student_df]]
  }
  if (!is.na(family$gengamma_Q)) {
    state$Q <- par$gengamma_Q[[family$gengamma_Q]]
  }
  if (name == "ordbeta") state$psi <- par$psi
  if (name == "censored_poisson") {
    state$upr <- fit$upr[i]
  }
  if (family$combine == "poisson_link") {
    # Poisson-link delta: component 1 models log numbers density, and
    # component 2 log weight, with the offset inside both means.
    offset <- fit$offset[i]
    state$log_one_minus_p <- -exp(offset + eta[i, 1L])
    state$log_p <- RTMB::logspace_sub(0, state$log_one_minus_p)
    state$mu <- if (m == 1L) exp(state$log_p)
      else exp(offset + eta[i, 1L] + eta[i, 2L] - state$log_p)
  } else if (isTRUE(spec$logit_mu)) {
    state$mu <- rtmb_logit_inverse_link(state$eta, state$link)
  } else {
    state$mu <- rtmb_inverse_link(state$eta, state$link)
  }
  if (isTRUE(spec$mixture)) {
    state$p_extreme <- theta$p_extreme
    state$mix_ratio <- theta$mix_ratio
  }
  state
}

# Evaluate or simulate every family component on its fitted rows. With
# `simulate_obs`, `OBS()` makes the response a simulation target; otherwise
# `y_i` reports response means instead.
rtmb_observations <- function(par, theta, prepared, eta) {
  "[<-" <- RTMB::ADoverload("[<-")
  fit <- prepared$fit
  n <- nrow(eta)
  out <- list()
  ln_phi_i <- NULL
  if (prepared$dispersion_model) {
    ln_phi_i <- rtmb_product(fit$Xdisp, par$b_disp_k)
    out$ln_phi_i <- ln_phi_i
    out$phi_i <- exp(ln_phi_i)
  }
  observing <- prepared$simulate_obs
  y_i <- if (observing) RTMB::OBS(fit$y) else fit$y
  simulating <- inherits(y_i, "simref")
  jnll_obs <- rep(0, n)
  devresid <- matrix(0, n, prepared$n_m)
  for (m in seq_len(prepared$n_m)) {
    for (group in prepared$obs_groups[[m]]) {
      family <- prepared$families[[group$f]]
      poisson_link <- family$combine == "poisson_link" && m == 1L
      spec <- rtmb_obs_family(family$family[[m]], poisson_link)
      state <- function(i) {
        rtmb_obs_state(i, m, family, eta, par, theta, prepared, ln_phi_i)
      }
      rows <- group$rows
      if (simulating) {
        # `[<-.simref` evaluates its indices in its own frame, so index once
        # and fill the resulting child reference.
        target <- y_i[rows, m]
        s <- state(rows)
        target[] <- spec$simulate(s$mu, s)
        next
      }
      if (!observing) {
        mu <- state(rows)$mu
        if (isTRUE(spec$logit_mu)) {
          mu <- RTMB::plogis(mu) * fit$size[rows]
        }
        y_i[rows, m] <- mu
        next
      }
      if (poisson_link) {
        # The C++ template records this residual for every row, treating a
        # missing response as zero.
        s <- state(rows)
        y <- fit$y[rows, m]
        devresid[rows, m] <- sqrt(-2 * spec$logpdf(ifelse(is.na(y), 0, y),
          s$mu, s))
      }
      i <- group$observed
      if (!length(i)) next
      s <- state(i)
      y <- fit$y[i, m]
      log_density <- spec$logpdf(y, s$mu, s)
      jnll_obs[i] <- jnll_obs[i] - fit$weights[i] * log_density
      if (!is.null(spec$deviance)) {
        # As in C++, some residuals can be NaN; keep them without R's
        # warnings.
        devresid[i, m] <- suppressWarnings(
          spec$deviance(y, s$mu, s, log_density))
      }
    }
  }
  c(out, list(y_i = y_i, jnll_obs = jnll_obs, devresid = devresid))
}
