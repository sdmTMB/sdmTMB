# Observation families for the RTMB backend. Each entry keeps one family's
# log density, simulator, and deviance residual together so they share one
# parameterization. Which families exist, and their auxiliary parameters, is
# declared separately in `.family_registry` (family-spec.R).
#
# Every function takes `mu`, the component mean, and `s`, the state built by
# rtmb_obs_state(), for any other parameter (`phi`, `ln_phi`, `size`, ...):
#
#   logpdf(y, mu, s)                 log density of each observation
#   simulate(mu, s)                  one draw per row, with plain numbers
#   deviance(y, mu, s, log_density)  deviance residuals; omit to report zero
#
# Optional fields:
#
#   logit_mu    TRUE if `mu` is logit(p) rather than the mean (binomial)
#   mean        function(mu, s): the response expectation, if not `mu`
#   log_nzprob  function(mu, phi): log P(Y > 0), for zero-truncated families
#   mixture     TRUE for two-component mixtures; `s` then carries
#               `p_extreme` and `mix_ratio`
#
# Use `mu`, not `s$mu`, inside these functions: the mixture wrapper calls the
# base family with the larger component's mean.

rtmb_obs_families <- list(
  gaussian = list(
    logpdf = function(y, mu, s) RTMB::dnorm(y, mu, s$phi, log = TRUE),
    simulate = function(mu, s) stats::rnorm(length(mu), mu, s$phi),
    deviance = function(y, mu, s, log_density) y - mu
  ),

  # `y` successes out of `size` trials. `mu` is logit(p), as in the C++
  # template, for accuracy at extreme probabilities.
  binomial = list(
    logit_mu = TRUE,
    logpdf = function(y, mu, s) RTMB::dbinom_robust(y, s$size, mu, log = TRUE),
    mean = function(mu, s) RTMB::plogis(mu) * s$size,
    simulate = function(mu, s) {
      stats::rbinom(length(mu), s$size, stats::plogis(mu))
    },
    deviance = function(y, mu, s, log_density) {
      n <- s$size
      xlogx <- function(x) ifelse(x > 0, x * log(x / n), 0) # data only
      log_p <- -RTMB::logspace_add(0, -mu)
      log_one_minus_p <- -RTMB::logspace_add(0, mu)
      deviance <- 2 * (xlogx(y) - y * log_p + xlogx(n - y) -
        (n - y) * log_one_minus_p)
      sign(y - n * exp(log_p)) * sqrt(deviance)
    }
  ),

  # Beta-binomial on shapes p * phi and (1 - p) * phi, with p from `eta`.
  betabinomial = list(
    logpdf = function(y, mu, s) rtmb_dbetabinom(y, s),
    mean = function(mu, s) mu * s$size,
    simulate = function(mu, s) {
      shape <- rtmb_betabinom_shapes(s)
      stats::rbinom(length(mu), s$size,
        stats::rbeta(length(mu), shape$a, shape$b))
    }
  ),

  # Beta-binomial with counts censored to [y, upr]; `upr` Inf means `size`.
  censored_betabinomial = list(
    logpdf = function(y, mu, s) rtmb_dcensbetabinom(y, s),
    mean = function(mu, s) mu * s$size,
    simulate = function(mu, s) {
      shape <- rtmb_betabinom_shapes(s)
      stats::rbinom(length(mu), s$size,
        stats::rbeta(length(mu), shape$a, shape$b))
    }
  ),

  # Binomial with counts censored to [y, upr]; `upr` Inf means `size`.
  censored_binomial = list(
    logit_mu = TRUE,
    logpdf = function(y, mu, s) rtmb_dcensbinom(y, mu, s),
    mean = function(mu, s) RTMB::plogis(mu) * s$size,
    simulate = function(mu, s) {
      stats::rbinom(length(mu), s$size, stats::plogis(mu))
    }
  ),

  poisson = list(
    logpdf = function(y, mu, s) RTMB::dpois(y, mu, log = TRUE),
    simulate = function(mu, s) stats::rpois(length(mu), mu),
    deviance = function(y, mu, s, log_density) {
      sign(y - mu) * sqrt(2 * (y * log((1e-10 + y) / mu) - (y - mu)))
    }
  ),

  censored_poisson = list(
    logpdf = function(y, mu, s) rtmb_dcenspois(y, mu, s$upr),
    simulate = function(mu, s) stats::rpois(length(mu), mu),
    deviance = function(y, mu, s, log_density) {
      rtmb_censpois_devresid(y, mu, s$upr, log_density)
    }
  ),

  # Variance mu * (1 + phi).
  nbinom1 = list(
    logpdf = function(y, mu, s) {
      RTMB::dnbinom_robust(y, log(mu), log(mu) + s$ln_phi, log = TRUE)
    },
    simulate = function(mu, s) {
      stats::rnbinom(length(mu), size = mu / s$phi, mu = mu)
    },
    deviance = function(y, mu, s, log_density) {
      rtmb_nbinom1_deviance(y, mu, s, log_density)
    }
  ),

  # Variance mu + mu^2 / phi.
  nbinom2 = list(
    logpdf = function(y, mu, s) {
      RTMB::dnbinom_robust(y, log(mu), 2 * log(mu) - s$ln_phi, log = TRUE)
    },
    simulate = function(mu, s) {
      stats::rnbinom(length(mu), size = s$phi, mu = mu)
    },
    deviance = function(y, mu, s, log_density) {
      rtmb_nb_deviance(y, mu, log_density, s$ln_phi)
    }
  ),

  # Compound Poisson-gamma with power `tweedie_p` in (1, 2).
  tweedie = list(
    logpdf = function(y, mu, s) {
      RTMB::dtweedie(y, mu, s$phi, s$tweedie_p, log = TRUE)
    },
    simulate = function(mu, s) {
      # A Poisson number of gamma summands.
      n <- length(mu)
      p <- s$tweedie_p
      count <- stats::rpois(n, mu^(2 - p) / (s$phi * (2 - p)))
      stats::rgamma(n, shape = count * (2 - p) / (p - 1),
        scale = s$phi * (p - 1) * mu^(p - 1))
    },
    deviance = function(y, mu, s, log_density) {
      p <- s$tweedie_p
      deviance <- 2 * (y^(2 - p) / (1 - p) / (2 - p) -
        y * mu^(1 - p) / (1 - p) + mu^(2 - p) / (2 - p))
      sign(y - mu) * sqrt(deviance)
    }
  ),

  # Shape phi, mean mu.
  Gamma = list(
    logpdf = function(y, mu, s) {
      RTMB::dgamma(y, shape = s$phi, scale = mu / s$phi, log = TRUE)
    },
    simulate = function(mu, s) {
      stats::rgamma(length(mu), shape = s$phi, scale = mu / s$phi)
    },
    deviance = function(y, mu, s, log_density) {
      sign(y - mu) * sqrt(2 * ((y - mu) / mu - log(y / mu)))
    }
  ),

  # Bias-corrected so that the mean is mu; phi is the SD of log(y).
  lognormal = list(
    logpdf = function(y, mu, s) {
      RTMB::dnorm(log(y), log(mu) - s$phi^2 / 2, s$phi, log = TRUE) - log(y)
    },
    simulate = function(mu, s) {
      exp(stats::rnorm(length(mu), log(mu) - s$phi^2 / 2, s$phi))
    },
    deviance = function(y, mu, s, log_density) {
      log(y) - (log(mu) - s$phi^2 / 2)
    }
  ),

  gengamma = list(
    logpdf = function(y, mu, s) rtmb_dgengamma(y, mu, s$phi, s$Q),
    simulate = function(mu, s) {
      w <- log(stats::rgamma(length(mu), s$Q^-2, 1))
      exp(w / (s$Q / s$phi) + rtmb_gengamma_log_theta(mu, s$phi, s$Q))
    },
    deviance = function(y, mu, s, log_density) {
      rtmb_gengamma_devresid(y, mu, s$phi, s$Q)
    }
  ),

  # Location mu, scale phi, `df` degrees of freedom.
  student = list(
    logpdf = function(y, mu, s) {
      RTMB::dt((y - mu) / s$phi, s$df, log = TRUE) - log(s$phi)
    },
    simulate = function(mu, s) mu + s$phi * stats::rt(length(mu), s$df),
    deviance = function(y, mu, s, log_density) {
      residual <- (y - mu) / s$phi
      sign(y - mu) * sqrt((s$df + 1) * log(1 + residual^2 / s$df))
    }
  ),

  # Mean mu, precision phi.
  Beta = list(
    logpdf = function(y, mu, s) {
      RTMB::dbeta(y, mu * s$phi, (1 - mu) * s$phi, log = TRUE)
    },
    simulate = function(mu, s) {
      stats::rbeta(length(mu), mu * s$phi, (1 - mu) * s$phi)
    }
  ),

  ordbeta = list(
    logpdf = function(y, mu, s) rtmb_dordbeta(y, mu, s),
    mean = function(mu, s) {
      p0 <- RTMB::plogis(s$psi[[1L]] - s$eta)
      p1 <- RTMB::plogis(s$eta - s$psi[[2L]])
      p1 + (1 - p0 - p1) * mu
    },
    simulate = function(mu, s) {
      # One uniform selects the component: P(0) = p0, P(1) = p1.
      n <- length(mu)
      u <- stats::runif(n)
      p0 <- stats::plogis(s$psi[[1L]] - s$eta)
      p1 <- stats::plogis(s$eta - s$psi[[2L]])
      ifelse(u < p0, 0, ifelse(u > 1 - p1, 1,
        stats::rbeta(n, mu * s$phi, (1 - mu) * s$phi)))
    }
  )
)

# Zero-truncated version of a negative binomial `base`. `size(mu, phi)` is
# the base family's size parameter, for simulation.
rtmb_truncated_nb <- function(base, size, log_nzprob) {
  force(base)
  force(size)
  force(log_nzprob)
  list(
    log_nzprob = log_nzprob,
    mean = function(mu, s) exp(log(mu) - log_nzprob(mu, s$phi)),
    logpdf = function(y, mu, s) {
      out <- base$logpdf(y, mu, s) - log_nzprob(mu, s$phi)
      out[y < 0.001] <- -Inf # zero counts are impossible
      out
    },
    simulate = function(mu, s) {
      # Invert the upper tail: draw a survival probability uniformly on
      # (0, P(Y > 0)) in log space, which stays finite when P(Y = 0) rounds
      # to one.
      size <- size(mu, s$phi)
      log_nonzero <- stats::pnbinom(0, size = size, mu = mu,
        lower.tail = FALSE, log.p = TRUE)
      stats::qnbinom(log(stats::runif(length(mu))) + log_nonzero,
        size = size, mu = mu, lower.tail = FALSE, log.p = TRUE)
    }
  )
}

# Mixture of `base` with mean mu (probability 1 - p_extreme) and with mean
# mu * mix_ratio (probability p_extreme).
rtmb_mixture <- function(base) {
  force(base)
  list(
    mixture = TRUE,
    mean = function(mu, s) (1 - s$p_extreme) * mu + s$p_extreme * mu * s$mix_ratio,
    logpdf = function(y, mu, s) {
      RTMB::logspace_add(log(1 - s$p_extreme) + base$logpdf(y, mu, s),
        log(s$p_extreme) + base$logpdf(y, mu * s$mix_ratio, s))
    },
    simulate = function(mu, s) {
      ifelse(stats::rbinom(length(mu), 1, s$p_extreme) == 0,
        base$simulate(mu, s), base$simulate(mu * s$mix_ratio, s))
    }
  )
}

rtmb_obs_families$truncated_nbinom1 <- rtmb_truncated_nb(
  rtmb_obs_families$nbinom1,
  size = function(mu, phi) mu / phi,
  log_nzprob = function(mu, phi) {
    RTMB::logspace_sub(0, -mu / phi * RTMB::logspace_add(0, log(phi)))
  })
rtmb_obs_families$truncated_nbinom2 <- rtmb_truncated_nb(
  rtmb_obs_families$nbinom2,
  size = function(mu, phi) phi,
  log_nzprob = function(mu, phi) {
    RTMB::logspace_sub(0, -phi * RTMB::logspace_add(0, log(mu) - log(phi)))
  })
rtmb_obs_families$censored_nbinom1 <- list(
  logpdf = function(y, mu, s) rtmb_dcensnb(y, mu, log(mu) + s$ln_phi, s$upr),
  simulate = rtmb_obs_families$nbinom1$simulate
)
rtmb_obs_families$censored_nbinom2 <- list(
  logpdf = function(y, mu, s) {
    rtmb_dcensnb(y, mu, 2 * log(mu) - s$ln_phi, s$upr)
  },
  simulate = rtmb_obs_families$nbinom2$simulate
)
rtmb_obs_families$gamma_mix <- rtmb_mixture(rtmb_obs_families$Gamma)
rtmb_obs_families$lognormal_mix <- rtmb_mixture(rtmb_obs_families$lognormal)
rtmb_obs_families$nbinom2_mix <- rtmb_mixture(rtmb_obs_families$nbinom2)

rtmb_mixture_families <- names(Filter(function(spec) isTRUE(spec$mixture),
  rtmb_obs_families))

# Component 1 of a Poisson-link delta model. Not a separate family:
# rtmb_obs_state() supplies the encounter probability as `mu`, with its logs
# `log_p` and `log_one_minus_p`.
rtmb_poisson_link_binomial <- list(
  logpdf = function(y, mu, s) {
    rtmb_ifelse_positive(y, s$log_p, s$log_one_minus_p)
  },
  # Positive for an encounter, negative for a zero
  deviance = function(y, mu, s, log_density) {
    ifelse(y > 0, 1, -1) * sqrt(-2 * log_density)
  },
  simulate = function(mu, s) stats::rbinom(length(mu), s$size, mu)
)

rtmb_obs_family <- function(name, poisson_link = FALSE) {
  if (poisson_link) return(rtmb_poisson_link_binomial)
  spec <- rtmb_obs_families[[name]]
  if (is.null(spec)) {
    cli_abort("Family {.val {name}} is not implemented in the RTMB backend.")
  }
  spec
}

# Distribution helpers ------------------------------------------------------

# Deviance residuals of a negative binomial with log size `log_theta`.
rtmb_nb_deviance <- function(y, mu, log_density, log_theta) {
  log_y <- log(y + 1e-10)
  saturated <- RTMB::dnbinom_robust(y, log_y, 2 * log_y - log_theta,
    log = TRUE)
  sign(y - mu) * sqrt(2 * (saturated - log_density))
}

# NB1 deviance residuals with phi fixed, following the C++
# `devresid_nbinom1()`. With size r = mu / phi, the saturated r solves
# digamma(y + r) - digamma(r) = log(1 + phi). The deviance is computed from
# gamma-function differences in y, which stay accurate as phi -> 0 (the
# Poisson limit, where r is huge), and takes its sign from the direction of
# the saturated mean, which differs from y.
rtmb_nbinom1_deviance <- function(y, mu, s, log_density) {
  positive <- y > 0
  y1 <- ifelse(positive, y, 1) # keeps the Newton steps finite for y = 0
  c <- log1p(s$phi)
  r_fit <- mu / s$phi
  # Newton steps on 1 / (digamma(y + r) - digamma(r)), which is increasing
  # and concave in r, from the exact y = 1 solution approach the root from
  # below. 10 steps reach machine precision for phi in [1e-4, 1e6] and y up
  # to 1e6; 15 leaves a margin.
  r <- 1 / c
  for (k in 1:15) {
    D <- rtmb_digamma_diff(r, y1)
    r <- r + (1 / D - 1 / c) * D * D / rtmb_trigamma_diff(r, y1)
  }
  deviance <- 2 * (rtmb_lgamma_diff(r, y1) - rtmb_lgamma_diff(r_fit, y1) -
    c * (r - r_fit))
  direction <- sign(r - r_fit)
  # For y = 0 the saturated mean -> 0
  deviance[!positive] <- 2 * (c * r_fit)[!positive]
  direction[!positive] <- -1
  deviance <- (deviance + abs(deviance)) / 2 # remove rounding below zero
  direction * sqrt(deviance)
}

# Differences in gamma functions over y >= 0 for AD types: lgamma(x + y) -
# lgamma(x), digamma(x + y) - digamma(x), and trigamma(x + y) - trigamma(x).
# Each is differenced term by term, so nothing cancels when x is large. As
# in the C++ `lgamma_diff()` etc.
#
# The recurrence Gamma(x + 1) = x Gamma(x) (the loops over k) shifts the
# argument to z = x + 10, where the asymptotic series below have relative
# error < 1e-12. Their coefficients come from the Bernoulli numbers
# B2 = 1/6, B4 = -1/30, B6 = 1/42, B8 = -1/30 (Abramowitz & Stegun 6.1.40,
# 6.3.18, 6.4.12):
#   lgamma(z)   ~ (z - 1/2) log(z) - z + log(2 pi) / 2
#                 + 1/(12 z) - 1/(360 z^3) + 1/(1260 z^5) - 1/(1680 z^7)
#   digamma(z)  ~ log(z) - 1/(2 z)
#                 - 1/(12 z^2) + 1/(120 z^4) - 1/(252 z^6) + 1/(240 z^8)
#   trigamma(z) ~ 1/z + 1/(2 z^2)
#                 + 1/(6 z^3) - 1/(30 z^5) + 1/(42 z^7) - 1/(30 z^9)
# Signs flip in the code because `rtmb_inv_pow_diff(z, y, k)` is
# 1 / z^k - 1 / (z + y)^k.
rtmb_inv_pow_diff <- function(z, y, k) {
  w <- z + y
  out <- 0
  for (j in 0:(k - 1)) out <- out + z^-(k - 1 - j) * w^-j
  out * y / (z * w)
}

rtmb_lgamma_diff <- function(x, y) {
  d <- function(k) rtmb_inv_pow_diff(x + 10, y, k)
  z <- x + 10
  out <- (z + y - 0.5) * log1p(y / z) + y * log(z) - y - d(1) / 12 +
    d(3) / 360 - d(5) / 1260 + d(7) / 1680
  for (k in 0:9) out <- out - log1p(y / (x + k))
  out
}

rtmb_digamma_diff <- function(x, y) {
  d <- function(k) rtmb_inv_pow_diff(x + 10, y, k)
  out <- log1p(y / (x + 10)) + d(1) / 2 + d(2) / 12 - d(4) / 120 +
    d(6) / 252 - d(8) / 240
  for (k in 0:9) out <- out + rtmb_inv_pow_diff(x + k, y, 1)
  out
}

rtmb_trigamma_diff <- function(x, y) {
  d <- function(k) rtmb_inv_pow_diff(x + 10, y, k)
  out <- -d(1) - d(2) / 2 - d(3) / 6 + d(5) / 30 - d(7) / 42 + d(9) / 30
  for (k in 0:9) out <- out - rtmb_inv_pow_diff(x + k, y, 2)
  out
}

# Choose between two AD vectors by an observed (data) condition.
rtmb_ifelse_positive <- function(y, positive, other) {
  out <- other
  out[y > 0] <- positive[y > 0]
  out
}

# Right-censored Poisson (`upr = Inf`), interval-censored Poisson
# (`y <= count <= upr`), or an exact count (`upr == y`).
rtmb_dcenspois <- function(y, lambda, upr) {
  out <- lambda * 0
  exact <- is.finite(upr) & upr == y
  out[exact] <- RTMB::dpois(y[exact], lambda[exact], log = TRUE)
  if (any(!exact)) {
    U <- ifelse(is.finite(upr), upr, Inf)[!exact]
    out[!exact] <- rtmb_censpois_logprob(y[!exact], U)(lambda[!exact])
  }
  out
}

# Censored Poisson deviance residuals. Exact counts give the Poisson
# deviance. Otherwise P(L <= Y <= U) is maximized as lambda -> Inf (U = Inf)
# or lambda -> 0 (L = 0), with saturated log likelihood 0, or else where its
# derivative p(L - 1) - p(U) = 0, at lambda^(U - L + 1) = U! / (L - 1)!. The
# saturated values depend only on data.
rtmb_censpois_devresid <- function(y, lambda, upr, log_density) {
  exact <- is.finite(upr) & upr == y
  interval <- is.finite(upr) & !exact & y > 0
  log_lambda_sat <- ifelse(is.finite(upr), -Inf, Inf)
  log_lambda_sat[interval] <- (lgamma(upr[interval] + 1) -
    lgamma(y[interval])) / (upr[interval] - y[interval] + 1)
  log_sat <- numeric(length(y))
  log_sat[interval] <- rtmb_censpois_logprob_value(
    exp(log_lambda_sat[interval]), y[interval], upr[interval])
  out <- sign(log_lambda_sat - log(lambda)) * sqrt(2 * (log_sat - log_density))
  out[exact] <- rtmb_obs_families$poisson$deviance(y[exact], lambda[exact])
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

# Q * w, where w = (log(x) - location) / sigma is the standardized log
# response.
rtmb_gengamma_qw <- function(x, mean, sigma, Q) {
  location <- rtmb_gengamma_log_theta(mean, sigma, Q) + log(Q^-2) * sigma / Q
  Q * (log(x) - location) / sigma
}

# The mean exists only if 1 + sigma * Q > 0; otherwise this returns NaN
# (log(v) - log(v) is an AD-safe check).
rtmb_dgengamma <- function(x, mean, sigma, Q) {
  k <- Q^-2
  qw <- rtmb_gengamma_qw(x, mean, sigma, Q)
  v <- 1 + sigma * Q
  -log(sigma * x) + 0.5 * log(Q^2) * (1 - 2 * k) + k * (qw - exp(qw)) -
    lgamma(k) + log(v) - log(v)
}

# The density peaks in the mean where qw = 0, so twice the log-likelihood
# ratio against the saturated model is 2 * Q^-2 * (exp(qw) - 1 - qw). Scaled
# by the dispersion sigma^2, as for Gamma() and lognormal(), this equals the
# Gamma deviance when Q = sigma and the lognormal deviance as Q -> 0.
rtmb_gengamma_devresid <- function(x, mean, sigma, Q) {
  qw <- rtmb_gengamma_qw(x, mean, sigma, Q)
  sign(qw / Q) * sigma / sqrt(Q^2) * sqrt(2 * (exp(qw) - 1 - qw))
}

# Beta-binomial on shape parameters p * phi and (1 - p) * phi.
rtmb_betabinom_shapes <- function(s) {
  logit_p <- rtmb_logit_inverse_link(s$eta, s$link)
  list(a = RTMB::plogis(logit_p) * s$phi, b = RTMB::plogis(-logit_p) * s$phi)
}

rtmb_dbetabinom <- function(y, s) {
  shape <- rtmb_betabinom_shapes(s)
  rtmb_lbetabinom(y, shape$a, shape$b, s$size)
}

# Beta-binomial log PMF of count `k` with shapes `a`, `b` and `n` trials.
rtmb_lbetabinom <- function(k, a, b, n) {
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + lgamma(a + b) +
    lgamma(k + a) + lgamma(n - k + b) - lgamma(n + a + b) - lgamma(a) -
    lgamma(b)
}

# log of a count PMF summed over counts k0, ..., k0 + len - 1 (len >= 1) for
# rows `i` of the parameter vectors in `d`, with successive terms from the
# PMF ratio p(k + 1) / p(k): `d$lpmf(k, i)` is the log PMF and
# `d$lratio(k, i)` the log ratio. Term `j` of every run longer than `j` is
# added at once, so the loop runs over run lengths rather than rows.
rtmb_count_logsum <- function(i, k0, len, d) {
  "[<-" <- RTMB::ADoverload("[<-")
  term <- d$lpmf(k0, i)
  out <- term
  for (j in seq_len(max(len, 1) - 1)) {
    r <- which(len > j)
    term[r] <- term[r] + d$lratio(k0[r] + j - 1, i[r])
    out[r] <- RTMB::logspace_add(out[r], term[r])
  }
  out
}

# log P(y <= Y <= upr) for a count Y on 0, ..., n with log PMF and ratio
# `d` (see rtmb_count_logsum()), with `upr` NA meaning `n` and `upr == y` an
# exact count. Rows with `cens_direct == 0` sum whichever side has fewer
# terms: the interval, or its complement as log(1 - P(outside)). The side is
# chosen from the data, so the tape is fixed. The complement loses accuracy
# when the interval probability is small, so after fitting,
# check_censored_betabinomial() switches any row that disagrees with the direct sum
# over the interval to `cens_direct = 1`. The PMF is summed rather than using
# 1 - F(y - 1), which cancels badly.
rtmb_dcenscount <- function(y, n, upr, cens_direct, d) {
  "[<-" <- RTMB::ADoverload("[<-")
  upr <- ifelse(is.finite(upr), upr, n) # whole numbers (normalized in R)
  full <- y == 0 & upr >= n # the whole support, set to 0 below
  n_low <- y # counts 0, ..., y - 1
  n_high <- n - upr # counts upr + 1, ..., n
  comp <- !full & cens_direct == 0 & n_low + n_high < upr - y + 1
  out <- rtmb_count_logsum(seq_along(y), y,
    ifelse(full | comp, 1, upr - y + 1), d)
  if (any(comp)) {
    lo <- which(comp & n_low > 0)
    hi <- which(comp & n_high > 0)
    both <- n_low[hi] > 0
    outside <- out
    outside[lo] <- rtmb_count_logsum(lo, rep(0, length(lo)), n_low[lo], d)
    high <- rtmb_count_logsum(hi, upr[hi] + 1, n_high[hi], d)
    outside[hi[!both]] <- high[!both]
    outside[hi[both]] <- RTMB::logspace_add(outside[hi[both]], high[both])
    out[comp] <- RTMB::logspace_sub(0 * outside[comp], outside[comp])
  }
  out[full] <- 0
  out
}

rtmb_dcensbetabinom <- function(y, s) {
  shape <- rtmb_betabinom_shapes(s)
  a <- shape$a
  b <- shape$b
  n <- s$size
  rtmb_dcenscount(y, n, s$upr, s$cens_direct, list(
    lpmf = function(k, i) rtmb_lbetabinom(k, a[i], b[i], n[i]),
    lratio = function(k, i) {
      log(k + a[i]) - log(n[i] - k - 1 + b[i]) + log(n[i] - k) - log(k + 1)
    }
  ))
}

# Binomial with `logit_p` the logit of the success probability, always summed
# over the interval: binomial tails are thin enough that the complement can
# round to 1 away from the estimate (e.g., in the inner optimization), giving
# log(0) before any precision check. RTMB's pbeta() isn't an option either,
# since its higher derivatives (needed by the Laplace approximation) aren't
# finite.
rtmb_dcensbinom <- function(y, logit_p, s) {
  n <- s$size
  rtmb_dcenscount(y, n, s$upr, 1, list(
    lpmf = function(k, i) RTMB::dbinom_robust(k, n[i], logit_p[i], log = TRUE),
    lratio = function(k, i) log(n[i] - k) - log(k + 1) + logit_p[i]
  ))
}

# Negative binomial counts censored to [y, upr] with mean `mu` and
# log(Var - mu) `log_vmm`; `upr` Inf is right-censored and a non-integer
# bound means the largest count it allows. Bounded intervals are summed
# directly. Right tails use P(Y > 0) * P(Y >= y | Y > 0): conditioning keeps
# the complement accurate when positives are rare but have a heavy tail.
# For small conditional tails, blend into a direct upper-tail sum instead.
# RTMB's pbeta() isn't an option: its higher derivatives (needed by the
# Laplace approximation) aren't finite.
rtmb_dcensnb <- function(y, mu, log_vmm, upr) {
  "[<-" <- RTMB::ADoverload("[<-")
  log_mu <- log(mu)
  # default non-censored density:
  out <- RTMB::dnbinom_robust(y, log_mu, log_vmm, log = TRUE)
  upr <- floor(upr)
  cens <- which(upr != y)
  if (!length(cens)) return(out)
  log_mu <- log_mu[cens]
  log_p <- -RTMB::logspace_add(0, log_vmm[cens] - log_mu)
  r <- exp(2 * log_mu - log_vmm[cens]) # size
  log_p0 <- r * log_p
  # Starting PMFs and ratios avoid differences of lgamma() or nearly equal
  # logs. Those lose both values and derivatives as size -> Inf (Poisson).
  d <- list(
    lpmf = function(k, i) {
      lp <- log_p0[i] + k * (log_mu[i] + log_p[i]) - lgamma(k + 1)
      for (j in seq_len(max(k, 1) - 1)) {
        use <- which(k > j)
        lp[use] <- lp[use] + log1p(j / r[i[use]])
      }
      lp
    },
    lratio = function(k, i) {
      log_mu[i] + log_p[i] + log1p(k / r[i]) - log(k + 1)
    }
  )
  y <- y[cens]
  upr <- upr[cens]
  val <- 0 * log_p # right-censored zeros stay 0; use log_p to make AD vector type
  bounded <- which(is.finite(upr))
  if (length(bounded)) { # sum PMF y<=Y<=upr
    val[bounded] <- rtmb_count_logsum(bounded, y[bounded],
      upr[bounded] - y[bounded] + 1, d)
  }
  # now do full right censored ones first for accuracy:
  log_positive <- RTMB::logspace_sub(0 * log_p0, log_p0)
  one <- which(!is.finite(upr) & y == 1)
  val[one] <- log_positive[one]
  # now for > 1:
  tail <- which(!is.finite(upr) & y > 1)
  if (length(tail)) {
    k <- y[tail]
    below <- rtmb_count_logsum(tail, 0 * k + 1, k - 1, d) - log_positive[tail]
    # Keep the unused complement finite after rounding of the lower sum.
    # A steep blend makes this floor negligible when the direct sum is used.
    below <- -sqrt(below^2 + 1e-32) # keep below negative; a smooth version of -abs(below)
    comp <- RTMB::logspace_sub(0 * below, below)
    # Switch before complement cancellation matters. Extend the direct sum
    # well past y to cover its tail throughout the smooth transition at 1e-3.
    direct <- rtmb_count_logsum(tail, k, pmax(200, 4 * k), d)
    # 6 sets width of transition
    # here weight goes from about 0.0025 to 0.9975
    w <- RTMB::plogis(6 * (comp - log(1e-3)))
    val[tail] <- w * (log_positive[tail] + comp) + (1 - w) * direct
  }
  out[cens] <- val
  out
}

# Ordered beta logit-scale cutpoints from the `psi` parameter
# c(lower cutpoint, log(upper - lower)), which keeps them ordered.
ordbeta_cutpoints <- function(psi) c(psi[[1L]], psi[[1L]] + exp(psi[[2L]]))

# Ordered beta (Kubinec 2023) with logit-scale cutpoints psi[1] < psi[2].
rtmb_dordbeta <- function(y, mu, s) {
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
  phi <- rtmb_pick(s$phi, mid)
  out[mid] <- RTMB::logspace_sub(log_F0, log_F1) +
    RTMB::dbeta(y[mid], mu[mid] * phi, (1 - mu[mid]) * phi, log = TRUE)
  out
}
