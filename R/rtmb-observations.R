rtmb_mixture_families <- c("gamma_mix", "lognormal_mix", "nbinom2_mix")

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

# Log probability of a nonzero count for zero-truncated negative binomials.
rtmb_log_nzprob <- function(family, mu, phi) {
  switch(family,
    truncated_nbinom1 = RTMB::logspace_sub(0,
      -mu / phi * RTMB::logspace_add(0, log(phi))),
    truncated_nbinom2 = RTMB::logspace_sub(0,
      -phi * RTMB::logspace_add(0, log(mu) - log(phi))),
    0)
}

# Response mean of one component, as used for combined projections. The
# truncated negative binomials report the mean of the truncated distribution.
rtmb_component_mean <- function(eta, family, m, theta) {
  mu <- rtmb_inverse_link(eta, family$link[[m]])
  name <- family$family[[m]]
  if (name %in% c("truncated_nbinom1", "truncated_nbinom2")) {
    mu <- mu / exp(rtmb_log_nzprob(name, mu, theta$phi[[family$phi]]))
  }
  mu
}

# Pick rows `i` from a row-aligned vector; scalars apply to every row.
rtmb_pick <- function(x, i) if (length(x) == 1L) x else x[i]

# Parameters and means of component `m` for fitted rows `i` of one family.
# For binomial components `mu` is on the logit scale, as in the C++ template.
rtmb_obs_state <- function(i, m, family, eta, par, theta, prepared, ln_phi_i) {
  fit <- prepared$fit
  name <- family$family[[m]]
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
  } else if (name == "binomial") {
    state$mu <- rtmb_logit_inverse_link(state$eta, state$link)
  } else {
    state$mu <- rtmb_inverse_link(state$eta, state$link)
  }
  if (name %in% rtmb_mixture_families) {
    state$p_extreme <- theta$p_extreme
    state$mu_large <- state$mu * theta$mix_ratio
  }
  state
}

# Observation log densities. Branches on the observed response are static;
# each branch is evaluated only on its own rows.
rtmb_obs_log_density <- function(family, y, s, poisson_link) {
  lognormal <- function(mu) {
    RTMB::dnorm(log(y), log(mu) - s$phi^2 / 2, s$phi, log = TRUE) - log(y)
  }
  nbinom2 <- function(mu) {
    RTMB::dnbinom_robust(y, log(mu), 2 * log(mu) - s$ln_phi, log = TRUE)
  }
  nbinom1 <- function(mu) {
    RTMB::dnbinom_robust(y, log(mu), log(mu) + s$ln_phi, log = TRUE)
  }
  gamma <- function(mu) {
    RTMB::dgamma(y, shape = s$phi, scale = mu / s$phi, log = TRUE)
  }
  # Zero counts are impossible under zero truncation.
  truncated <- function(log_density) {
    log_density <- log_density - rtmb_log_nzprob(family, s$mu, s$phi)
    log_density[y < 0.001] <- -Inf
    log_density
  }
  mixture <- function(density) {
    RTMB::logspace_add(log(1 - s$p_extreme) + density(s$mu),
      log(s$p_extreme) + density(s$mu_large))
  }
  switch(family,
    gaussian = RTMB::dnorm(y, s$mu, s$phi, log = TRUE),
    binomial = if (poisson_link) {
      rtmb_ifelse_positive(y, s$log_p, s$log_one_minus_p)
    } else {
      RTMB::dbinom_robust(y, s$size, s$mu, log = TRUE)
    },
    tweedie = RTMB::dtweedie(y, s$mu, s$phi, s$tweedie_p, log = TRUE),
    poisson = RTMB::dpois(y, s$mu, log = TRUE),
    Gamma = gamma(s$mu),
    nbinom2 = nbinom2(s$mu),
    lognormal = lognormal(s$mu),
    student = RTMB::dt((y - s$mu) / s$phi, s$df, log = TRUE) - log(s$phi),
    Beta = RTMB::dbeta(y, s$mu * s$phi, (1 - s$mu) * s$phi, log = TRUE),
    truncated_nbinom2 = truncated(nbinom2(s$mu)),
    nbinom1 = nbinom1(s$mu),
    truncated_nbinom1 = truncated(nbinom1(s$mu)),
    censored_poisson = rtmb_dcenspois(y, s$mu, s$upr),
    gamma_mix = mixture(gamma),
    lognormal_mix = mixture(lognormal),
    nbinom2_mix = mixture(nbinom2),
    gengamma = rtmb_dgengamma(y, s$mu, s$phi, s$Q),
    betabinomial = rtmb_dbetabinom(y, s),
    ordbeta = rtmb_dordbeta(y, s))
}

# Choose between two AD vectors by an observed (data) condition.
rtmb_ifelse_positive <- function(y, positive, other) {
  out <- other
  out[y > 0] <- positive[y > 0]
  out
}

# Right-censored Poisson (`upr = NA`), interval-censored Poisson
# (`y <= count <= upr`), or an exact count (`upr == y`).
rtmb_dcenspois <- function(y, lambda, upr) {
  out <- lambda * 0
  exact <- !is.na(upr) & upr == y
  out[exact] <- RTMB::dpois(y[exact], lambda[exact], log = TRUE)
  if (any(!exact)) {
    U <- ifelse(is.na(upr), Inf, upr)[!exact]
    out[!exact] <- rtmb_censpois_logprob(y[!exact], U)(lambda[!exact])
  }
  out
}

# AD function of `lambda` returning log P(L <= Y <= U), Y ~ Poisson(lambda),
# with fixed bounds; `U = Inf` is right censoring. Mirrors the C++ atomic
# `censpois_logprob()`. `RTMB::ppois()` has no log or upper-tail AD method,
# so the value is computed with plain numbers at every evaluation and the
# derivative is supplied: with y = log P,
# dy/dlambda = exp(log p(L - 1) - y) - exp(log p(U) - y), omitting the terms
# for L = 0 or U = Inf. The rule uses AD functions, so higher-order
# derivatives work.
rtmb_censpois_logprob <- function(L, U) {
  lower <- L > 0
  upper <- is.finite(U)
  RTMB::ADjoint(
    function(lambda) rtmb_censpois_logprob_value(lambda, L, U),
    function(lambda, y, dy) {
      d <- lambda * 0
      d[lower] <- exp(RTMB::dpois(L[lower] - 1, lambda[lower], log = TRUE) -
        y[lower])
      d[upper] <- d[upper] - exp(RTMB::dpois(U[upper], lambda[upper],
        log = TRUE) - y[upper])
      dy * d
    })
}

# Plain-number log P(L <= Y <= U) using the tail that avoids cancellation.
rtmb_censpois_logprob_value <- function(lambda, L, U) {
  log1mexp <- function(x) ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
  log_diff <- function(a, b) a + log1mexp(b - a)
  log_cdf <- function(q, i, lower.tail = TRUE) {
    stats::ppois(q, lambda[i], lower.tail = lower.tail, log.p = TRUE)
  }
  out <- rep(-Inf, length(lambda))
  right <- is.infinite(U)
  out[right & L <= 0] <- 0
  i <- right & L > 0
  out[i] <- log_cdf(L[i] - 1, i, lower.tail = FALSE)
  i <- !right & L <= 0
  out[i] <- log_cdf(U[i], i)
  wide <- !right & L > 0 & U - L > 64
  i <- wide & lambda < L # mass above the interval: difference of upper tails
  out[i] <- log_diff(log_cdf(L[i] - 1, i, FALSE), log_cdf(U[i], i, FALSE))
  i <- wide & lambda > U # mass below: difference of lower CDFs
  out[i] <- log_diff(log_cdf(U[i], i), log_cdf(L[i] - 1, i))
  i <- wide & lambda >= L & lambda <= U # 1 minus both tails
  out[i] <- log1p(-(exp(log_cdf(L[i] - 1, i)) +
    exp(log_cdf(U[i], i, FALSE))))
  # Narrow intervals (or a failed difference): sum the PMF.
  for (j in which(!right & L > 0 & !is.finite(out))) {
    ld <- stats::dpois(L[j]:U[j], lambda[j], log = TRUE)
    m <- max(ld)
    out[j] <- if (m == -Inf) -Inf else m + log(sum(exp(ld - m)))
  }
  out
}

# Generalized gamma with Prentice (1974) parameterization and a mean
# parameterization, following the C++ `dgengamma()`.
rtmb_gengamma_log_theta <- function(mean, sigma, Q) {
  k <- Q^-2
  beta <- Q / sigma
  log(mean) - lgamma((k * beta + 1) / beta) + lgamma(k)
}

rtmb_dgengamma <- function(x, mean, sigma, Q) {
  k <- Q^-2
  mu <- rtmb_gengamma_log_theta(mean, sigma, Q) + log(k) / (Q / sigma)
  qw <- Q * (log(x) - mu) / sigma
  -log(sigma * x) + 0.5 * log(Q^2) * (1 - 2 * k) + k * (qw - exp(qw)) -
    lgamma(k)
}

# Beta-binomial on shape parameters mu * phi and (1 - mu) * phi.
rtmb_betabinom_shapes <- function(s) {
  logit_p <- rtmb_logit_inverse_link(s$eta, s$link)
  list(a = RTMB::plogis(logit_p) * s$phi, b = RTMB::plogis(-logit_p) * s$phi)
}

rtmb_dbetabinom <- function(y, s) {
  shape <- rtmb_betabinom_shapes(s)
  a <- shape$a
  b <- shape$b
  n <- s$size
  lgamma(n + 1) - lgamma(y + 1) - lgamma(n - y + 1) + lgamma(a + b) +
    lgamma(y + a) + lgamma(n - y + b) - lgamma(n + a + b) - lgamma(a) -
    lgamma(b)
}

# Ordered beta (Kubinec 2023) with logit-scale cutpoints psi[1] < psi[2].
rtmb_dordbeta <- function(y, s) {
  out <- s$eta * 0
  zero <- y == 0
  one <- y == 1
  mid <- !zero & !one
  eta <- s$eta
  out[zero] <- -RTMB::logspace_add(0, eta[zero] - s$psi[[1L]])
  out[one] <- -RTMB::logspace_add(0, -(eta[one] - s$psi[[2L]]))
  e <- eta[mid]
  log_F0 <- -RTMB::logspace_add(0, -(e - s$psi[[1L]]))
  log_F1 <- -RTMB::logspace_add(0, -(e - s$psi[[2L]]))
  mu <- s$mu[mid]
  phi <- rtmb_pick(s$phi, mid)
  out[mid] <- RTMB::logspace_sub(log_F0, log_F1) +
    RTMB::dbeta(y[mid], mu * phi, (1 - mu) * phi, log = TRUE)
  out
}

# Observation draws. Simulation evaluates the objective with plain numbers,
# so base R generators apply. Weights do not rescale draws.
rtmb_obs_simulate <- function(family, s) {
  n <- length(s$mu)
  mu <- s$mu
  phi <- s$phi
  lognormal <- function(mu) exp(stats::rnorm(n, log(mu) - phi^2 / 2, phi))
  # Zero-truncated negative binomial by inverting the upper tail: draw a
  # survival probability uniformly on (0, P(Y > 0)) in log space, which
  # stays finite when P(Y = 0) rounds to one.
  truncated_nbinom <- function(size) {
    log_nonzero <- stats::pnbinom(0, size = size, mu = mu,
      lower.tail = FALSE, log.p = TRUE)
    stats::qnbinom(log(stats::runif(n)) + log_nonzero, size = size, mu = mu,
      lower.tail = FALSE, log.p = TRUE)
  }
  mixture <- function(draw) {
    ifelse(stats::rbinom(n, 1, s$p_extreme) == 0, draw(mu), draw(s$mu_large))
  }
  switch(family,
    gaussian = stats::rnorm(n, mu, phi),
    binomial = stats::rbinom(n, s$size,
      if (is.null(s$log_p)) stats::plogis(mu) else mu),
    tweedie = {
      # Compound Poisson-gamma: a Poisson number of gamma summands.
      p <- s$tweedie_p
      count <- stats::rpois(n, mu^(2 - p) / (phi * (2 - p)))
      stats::rgamma(n, shape = count * (2 - p) / (p - 1),
        scale = phi * (p - 1) * mu^(p - 1))
    },
    poisson = stats::rpois(n, mu),
    Gamma = stats::rgamma(n, shape = phi, scale = mu / phi),
    nbinom2 = stats::rnbinom(n, size = phi, mu = mu),
    lognormal = lognormal(mu),
    student = mu + phi * stats::rt(n, s$df),
    Beta = stats::rbeta(n, mu * phi, (1 - mu) * phi),
    truncated_nbinom2 = truncated_nbinom(phi),
    nbinom1 = stats::rnbinom(n, size = mu / phi, mu = mu),
    truncated_nbinom1 = truncated_nbinom(mu / phi),
    censored_poisson = stats::rpois(n, mu),
    gamma_mix = mixture(function(m) stats::rgamma(n, phi, scale = m / phi)),
    lognormal_mix = mixture(lognormal),
    nbinom2_mix = mixture(function(m) stats::rnbinom(n, size = phi, mu = m)),
    gengamma = {
      Q <- s$Q
      w <- log(stats::rgamma(n, Q^-2, 1))
      exp(w / (Q / phi) + rtmb_gengamma_log_theta(mu, phi, Q))
    },
    betabinomial = {
      shape <- rtmb_betabinom_shapes(s)
      stats::rbinom(n, s$size, stats::rbeta(n, shape$a, shape$b))
    },
    ordbeta = {
      # One uniform selects the component: P(0) = p0, P(1) = p1.
      u <- stats::runif(n)
      p0 <- stats::plogis(s$psi[[1L]] - s$eta)
      p1 <- stats::plogis(s$eta - s$psi[[2L]])
      ifelse(u < p0, 0, ifelse(u > 1 - p1, 1, stats::rbeta(n, mu * phi,
        (1 - mu) * phi)))
    })
}

# Deviance residuals for observed rows; families without one report zero.
rtmb_obs_deviance <- function(family, y, s, log_density) {
  nb_deviance <- function(log_theta) {
    log_y <- log(y + 1e-10)
    saturated <- RTMB::dnbinom_robust(y, log_y, 2 * log_y - log_theta,
      log = TRUE)
    sign(y - s$mu) * sqrt(2 * (saturated - log_density))
  }
  switch(family,
    gaussian = y - s$mu,
    binomial = {
      # `y` successes out of `size` trials; `mu` is logit(p).
      n <- s$size
      xlogx <- function(x) ifelse(x > 0, x * log(x / n), 0) # data only
      log_p <- -RTMB::logspace_add(0, -s$mu)
      log_one_minus_p <- -RTMB::logspace_add(0, s$mu)
      deviance <- 2 * (xlogx(y) - y * log_p + xlogx(n - y) -
        (n - y) * log_one_minus_p)
      sign(y - n * exp(log_p)) * sqrt(deviance)
    },
    tweedie = {
      p <- s$tweedie_p
      deviance <- 2 * (y^(2 - p) / (1 - p) / (2 - p) -
        y * s$mu^(1 - p) / (1 - p) + s$mu^(2 - p) / (2 - p))
      sign(y - s$mu) * sqrt(deviance)
    },
    poisson = sign(y - s$mu) *
      sqrt(2 * (y * log((1e-10 + y) / s$mu) - (y - s$mu))),
    Gamma = sign(y - s$mu) * sqrt(2 * ((y - s$mu) / s$mu - log(y / s$mu))),
    nbinom2 = nb_deviance(s$ln_phi),
    nbinom1 = nb_deviance(log(s$mu) - s$ln_phi),
    lognormal = log(y) - (log(s$mu) - s$phi^2 / 2),
    student = {
      residual <- (y - s$mu) / s$phi
      sign(y - s$mu) * sqrt((s$df + 1) * log(1 + residual^2 / s$df))
    },
    rep(0, length(y)))
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
      name <- family$family[[m]]
      poisson_link <- family$combine == "poisson_link" && m == 1L
      state <- function(i) {
        rtmb_obs_state(i, m, family, eta, par, theta, prepared, ln_phi_i)
      }
      rows <- group$rows
      if (simulating) {
        # `[<-.simref` evaluates its indices in its own frame, so index once
        # and fill the resulting child reference.
        target <- y_i[rows, m]
        target[] <- rtmb_obs_simulate(name, state(rows))
        next
      }
      if (!observing) {
        mu <- state(rows)$mu
        if (name == "binomial" && !poisson_link) {
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
        devresid[rows, m] <- sqrt(-2 * rtmb_ifelse_positive(
          ifelse(is.na(y), 0, y), s$log_p, s$log_one_minus_p))
      }
      i <- group$observed
      if (!length(i)) next
      s <- state(i)
      y <- fit$y[i, m]
      log_density <- rtmb_obs_log_density(name, y, s, poisson_link)
      jnll_obs[i] <- jnll_obs[i] - fit$weights[i] * log_density
      if (!poisson_link) {
        # As in C++, some residuals can be NaN; keep them without R's
        # warnings.
        devresid[i, m] <- suppressWarnings(
          rtmb_obs_deviance(name, y, s, log_density))
      }
    }
  }
  c(out, list(y_i = y_i, jnll_obs = jnll_obs, devresid = devresid))
}
