rtmb_family_data <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  mu <- exp(0.5 + 0.4 * d$z)
  d$tw <- ifelse(runif(n) < 0.3, 0, rgamma(n, 2, 2 / mu))
  d$cont <- 1 + d$z + 0.5 * rt(n, 4)
  d$prop <- rbeta(n, 2, 3)
  d$ordprop <- d$prop
  d$ordprop[seq_len(n / 10)] <- 0
  d$ordprop[n / 10 + seq_len(n / 10)] <- 1
  d$count <- pmax(rnbinom(n, size = 2, mu = mu), 1)
  d$count0 <- rnbinom(n, size = 2, mu = mu)
  d$pos <- rgamma(n, 2, 2 / mu)
  d$mixpos <- ifelse(runif(n) < 0.1, 5, 1) * d$pos
  d$trials <- rep(c(5, 10, 3), length.out = n)
  d$succ <- rbinom(n, d$trials, rbeta(n, 2, 3))
  d$hurdle_count <- ifelse(runif(n) < 0.4, 0, d$count)
  d$hurdle_prop <- ifelse(runif(n) < 0.4, 0, d$prop)
  d$hurdle_pos <- ifelse(runif(n) < 0.4, 0, d$pos)
  d
}

test_that("RTMB basic likelihoods match TMB at fixed parameters", {
  set.seed(103)
  n <- 36L
  x <- seq(-0.8, 0.8, length.out = n)
  offset <- rep(c(-0.1, 0.2), length.out = n)
  weights <- rep(c(1, 0.5, 2), length.out = n)
  cases <- list(
    list(family = gaussian(), y = rnorm(n)),
    list(family = gaussian(link = "log"), y = exp(0.3 + x) + rnorm(n, sd = 0.2)),
    list(family = gaussian(link = "inverse"), y = rnorm(n, mean = 0.5, sd = 0.2),
      par = c(2, 0.2, -1)),
    list(family = poisson(), y = rpois(n, exp(0.3 + x))),
    list(family = poisson(link = "identity"), y = rpois(n, 2 + x),
      par = c(2, 0.2)),
    list(family = nbinom1(), y = rnbinom(n, size = exp(0.3 + x) / 0.8,
      mu = exp(0.3 + x))),
    list(family = nbinom2(), y = rnbinom(n, size = 2.5,
      mu = exp(0.3 + x)))
  )
  for (case in cases) {
    d <- data.frame(x = x, y = case$y)
    d$y[[3L]] <- NA_real_
    common <- list(formula = y ~ x, data = d, spatial = "off",
      family = case$family, weights = weights, offset = offset, do_fit = FALSE)
    tmb <- do.call(sdmTMB, c(common,
      list(control = sdmTMBcontrol(multiphase = FALSE, backend = "tmb"))))
    rtmb <- do.call(sdmTMB, c(common,
      list(control = sdmTMBcontrol(backend = "rtmb", multiphase = FALSE))))
    p <- tmb$tmb_obj$par
    p[] <- if (is.null(case$par)) p + seq_along(p) * 0.07 else case$par
    expect_equal(rtmb$tmb_obj$fn(p), tmb$tmb_obj$fn(p), tolerance = 1e-7,
      info = case$family$family)
    expect_equal(rtmb$tmb_obj$gr(p), tmb$tmb_obj$gr(p), tolerance = 1e-6,
      info = case$family$family)
    expected <- tmb$tmb_obj$report(p)
    actual <- rtmb$tmb_obj$report(p)
    for (name in c("eta_fixed_i", "eta_i", "jnll_obs", "devresid")) {
      expect_equal(actual[[name]], expected[[name]], tolerance = 1e-6,
        info = paste(case$family$family, name))
    }
    if (case$family$family != "poisson") {
      expect_equal(actual$phi, expected$phi, tolerance = 1e-6,
        info = case$family$family)
    }
  }
})

test_that("RTMB matches TMB for other families", {
  d <- rtmb_family_data(90L, 81)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  d$dist <- rep(c("tw", "bb", "gg"), 30L)
  d$ymix <- ifelse(d$dist == "tw", d$tw,
    ifelse(d$dist == "bb", d$succ, d$hurdle_pos))
  cases <- list(
    tweedie = list(tw ~ z, family = tweedie()),
    student = list(cont ~ z, family = student(df = 3)),
    student_df = list(cont ~ z, family = student()),
    Beta = list(prop ~ z, family = Beta()),
    ordbeta = list(ordprop ~ z, family = ordbeta()),
    truncated_nbinom2 = list(count ~ z, family = truncated_nbinom2()),
    truncated_nbinom1 = list(count ~ z, family = truncated_nbinom1()),
    gamma_mix = list(mixpos ~ z, family = gamma_mix()),
    lognormal_mix = list(mixpos ~ z, family = lognormal_mix()),
    nbinom2_mix = list(count0 ~ z, family = nbinom2_mix()),
    gengamma = list(pos ~ z, family = gengamma()),
    betabinomial = list(cbind(succ, trials - succ) ~ z,
      family = betabinomial()),
    binomial_trials = list(cbind(succ, trials - succ) ~ z,
      family = binomial()),
    delta_truncated_nbinom2 = list(hurdle_count ~ z,
      family = delta_truncated_nbinom2()),
    delta_truncated_nbinom1 = list(hurdle_count ~ z,
      family = delta_truncated_nbinom1()),
    delta_gengamma = list(hurdle_pos ~ z, family = delta_gengamma()),
    delta_beta = list(hurdle_prop ~ z, family = delta_beta()),
    delta_gamma_mix = list(hurdle_pos ~ z, family = delta_gamma_mix()),
    delta_lognormal_mix = list(hurdle_pos ~ z,
      family = delta_lognormal_mix(type = "poisson-link")),
    multi_family = list(ymix ~ z, distribution_column = "dist",
      family = list(tw = tweedie(), bb = betabinomial(),
        gg = delta_gengamma()),
      weights = ifelse(d$dist == "bb", d$trials, 1))
  )
  project <- c("tweedie", "truncated_nbinom2", "gamma_mix",
    "delta_truncated_nbinom1", "delta_gamma_mix", "delta_lognormal_mix",
    "multi_family")
  for (name in names(cases)) {
    fit <- suppressMessages(do.call(sdmTMB, c(cases[[name]],
      list(data = d, mesh = mesh, do_fit = FALSE))))
    p <- rtmb_test_parameters(fit$tmb_params, fit$tmb_map)
    expect_rtmb_matches_tmb(fit$tmb_data, p, fit$tmb_map, fit$tmb_random,
      info = name)
    if (name %in% project) {
      projection <- predict(fit, newdata = d[c(1, 2, 3, 90), ],
        return_tmb_data = TRUE)
      expect_rtmb_matches_tmb(projection, p, fit$tmb_map, fit$tmb_random,
        info = paste(name, "projection"))
    }
  }
})

# The C++ template once evaluated the censored Poisson CDF with `asDouble()`,
# holding those terms fixed at the taping point. Compare both backends with a
# direct calculation away from it.
test_that("censored Poisson likelihood is correct away from the taping point", {
  skip_if_not_installed("numDeriv")
  set.seed(81)
  n <- 30L
  d <- data.frame(z = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$z))
  upr <- ifelse(seq_len(n) %% 3 == 0, NA,
    ifelse(seq_len(n) %% 3 == 1, d$y, d$y + 2))
  fit <- sdmTMB(y ~ z, data = d, spatial = "off", do_fit = FALSE,
    family = censored_poisson(),
    control = sdmTMBcontrol(censored_upper = upr))
  nll <- function(b) {
    lambda <- exp(b[[1L]] + b[[2L]] * d$z)
    ll <- ifelse(is.na(upr),
      ifelse(d$y == 0, 0, ppois(d$y - 1, lambda, lower.tail = FALSE,
        log.p = TRUE)),
      ifelse(upr > d$y,
        log(ppois(upr, lambda) - ifelse(d$y > 0, ppois(d$y - 1, lambda), 0)),
        dpois(d$y, lambda, log = TRUE)))
    -sum(ll)
  }
  for (backend in c("tmb", "rtmb")) {
    obj <- make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
      fit$tmb_random, backend = backend)
    for (b in list(c(0.3, 0.2), c(0.6, -0.1))) {
      expect_equal(obj$fn(b), nll(b), tolerance = 1e-8, info = backend)
      expect_equal(as.vector(obj$gr(b)), numDeriv::grad(nll, b),
        tolerance = 1e-6, info = backend)
    }
  }
})

# Independent tests for log P(L <= Y <= U), Y ~ Poisson(lambda), and its
# derivatives in eta = log(lambda), using r(k) = p(k) / P.
censpois_oracle <- function(lambda, L, U) {
  lp <- if (is.infinite(U)) {
    if (L == 0) 0 else ppois(L - 1, lambda, lower.tail = FALSE, log.p = TRUE)
  } else {
    ld <- dpois(L:U, lambda, log = TRUE)
    max(ld) + log(sum(exp(ld - max(ld))))
  }
  r <- function(k) {
    if (k < 0 || is.infinite(k)) 0 else exp(dpois(k, lambda, log = TRUE) - lp)
  }
  g <- r(L - 1) - r(U)
  dg <- (r(L - 2) - r(L - 1)) - (r(U - 1) - r(U)) - g^2
  c(value = lp, gradient = lambda * g, hessian = lambda * g + lambda^2 * dg)
}

# Rows far in either tail once gave -Inf because the log probability was
# taken after forming the CDF on the ordinary scale.
test_that("censored Poisson log probabilities are stable in the tails", {
  cases <- rbind(
    c(1, 100, Inf), c(1, 100, 102), c(1000, 0, 2), c(3, 0, Inf), c(3, 0, 5),
    c(3, 2, 6), c(1e-3, 50, Inf), c(50, 3, 3), c(5, 10, 500),
    c(1000, 10, 200), c(300, 200, 400), c(1e-8, 2, Inf), c(400, 2, Inf))
  upr <- ifelse(is.infinite(cases[, 3]), NA, cases[, 3])
  expect_equal(rtmb_dcenspois(cases[, 2], cases[, 1], upr),
    apply(cases, 1, function(x) censpois_oracle(x[1], x[2], x[3])[[1]]),
    tolerance = 1e-12)
  # The same tape at several log means, in both backends.
  for (k in seq_len(nrow(cases))) {
    lambda <- cases[k, 1]
    L <- cases[k, 2]
    d <- data.frame(y = L, o = log(lambda))
    fit <- sdmTMB(y ~ 1, offset = "o", data = d, spatial = "off",
      do_fit = FALSE, family = censored_poisson(),
      control = sdmTMBcontrol(censored_upper = upr[k]))
    for (backend in c("tmb", "rtmb")) {
      obj <- make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
        backend = backend)
      for (b in c(0, 0.3, -0.4)) {
        expected <- censpois_oracle(lambda * exp(b), L, cases[k, 3])
        got <- -c(obj$fn(b), obj$gr(b), obj$he(b))
        expect_lt(max(abs(got - expected) / pmax(1, abs(expected))), 1e-8,
          label = paste(backend, k, b))
      }
    }
  }
})

# The Laplace gradient needs third derivatives of the censored terms.
test_that("censored Poisson random-effect models are consistent", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n <- 60L
  d <- data.frame(g = factor(rep(1:10, each = 6)), z = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$z + rnorm(10, 0, 0.5)[d$g]))
  type <- seq_len(n) %% 3
  upr <- ifelse(type == 0, NA, ifelse(type == 1, d$y, d$y + 3))
  d$y[1:2] <- c(60, 80) # far upper tail at the starting values
  upr[1:2] <- c(NA, 85)
  fit <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off", do_fit = FALSE,
    family = censored_poisson(), control = sdmTMBcontrol(censored_upper = upr))
  fits <- list()
  for (backend in c("tmb", "rtmb")) {
    obj <- make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
      fit$tmb_random, backend = backend)
    for (p in list(obj$par, obj$par + c(0.3, -0.2, -0.5))) {
      expect_true(is.finite(obj$fn(p)), label = backend)
      expect_equal(as.vector(obj$gr(p)), numDeriv::grad(obj$fn, p),
        tolerance = 1e-6, info = backend)
    }
    fits[[backend]] <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off",
      family = censored_poisson(),
      control = sdmTMBcontrol(censored_upper = upr, backend = backend))
    expect_true(fits[[backend]]$sd_report$pdHess, label = backend)
  }
  expect_equal(tidy(fits$rtmb), tidy(fits$tmb), tolerance = 1e-6)
  expect_equal(tidy(fits$rtmb, "ran_pars"), tidy(fits$tmb, "ran_pars"),
    tolerance = 1e-6)
})

# Observation draws from each backend against closed-form means and
# variances at fixed parameters, at two linear-predictor values. Tolerances
# are 5 Monte Carlo standard errors; variance errors use the sample fourth
# central moment.
test_that("simulated observations match their distributions in both backends", {
  n <- 20000L
  d <- data.frame(x = rep(c(0, 1), each = n / 2), trials = 8)
  b_j <- c(0.3, 0.5)
  phi <- 1.5
  # Zero-truncated NB moments from the untruncated mean and variance.
  truncated <- function(size, v) function(mu) {
    q <- exp(stats::pnbinom(0, size = size(mu), mu = mu, lower.tail = FALSE,
      log.p = TRUE))
    c(mu / q, (v(mu) + mu^2) / q - (mu / q)^2)
  }
  # Two-component mixtures with P(large) = 0.1 and mean ratio 3.
  mixture <- function(v) function(mu) {
    m <- c(mu, 3 * mu)
    w <- c(0.9, 0.1)
    mean <- sum(w * m)
    c(mean, sum(w * (v(m) + m^2)) - mean^2)
  }
  mix_par <- list(logit_p_extreme = stats::qlogis(0.1), log_ratio_mix = log(2))
  nb1_var <- function(mu) mu * (1 + phi)
  nb2_var <- function(mu) mu + mu^2 / phi
  gamma_var <- function(mu) mu^2 / phi
  lognormal_var <- function(mu) mu^2 * (exp(0.5^2) - 1)
  # `moments()` takes the inverse-link mean; `NA` skips the variance check.
  cases <- list(
    gaussian = list(gaussian(), function(mu) c(mu, phi^2)),
    poisson = list(poisson(), function(mu) c(mu, mu)),
    nbinom1 = list(nbinom1(), function(mu) c(mu, nb1_var(mu))),
    nbinom2 = list(nbinom2(), function(mu) c(mu, nb2_var(mu))),
    truncated_nbinom1 = list(truncated_nbinom1(),
      truncated(function(mu) mu / phi, nb1_var)),
    truncated_nbinom2 = list(truncated_nbinom2(),
      truncated(function(mu) phi, nb2_var)),
    tweedie = list(tweedie(), function(mu) c(mu, phi * mu^1.5),
      par = list(thetaf = 0)),
    Gamma = list(Gamma(link = "log"), function(mu) c(mu, gamma_var(mu))),
    lognormal = list(lognormal(), function(mu) c(mu, lognormal_var(mu)),
      par = list(ln_phi = log(0.5))),
    student = list(student(df = 10), function(mu) c(mu, phi^2 * 10 / 8)),
    Beta = list(Beta(), function(mu) c(mu, mu * (1 - mu) / (phi + 1)),
      y = 0.5),
    betabinomial = list(betabinomial(),
      function(mu) c(8 * mu, 8 * mu * (1 - mu) * (phi + 8) / (phi + 1))),
    gengamma = list(gengamma(), function(mu) c(mu, NA),
      par = list(gengamma_Q = 0.5)),
    gamma_mix = list(gamma_mix(), mixture(gamma_var), par = mix_par),
    lognormal_mix = list(lognormal_mix(), mixture(lognormal_var),
      par = c(mix_par, ln_phi = log(0.5))),
    nbinom2_mix = list(nbinom2_mix(), mixture(nb2_var), par = mix_par)
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    family <- case[[1L]]
    d$y <- if (is.null(case$y)) 1 else case$y
    formula <- if (name == "betabinomial") cbind(y, trials - y) ~ x else y ~ x
    for (backend in c("tmb", "rtmb")) {
      obj <- suppressMessages(sdmTMB(formula, data = d, family = family,
        spatial = "off", do_fit = FALSE,
        control = sdmTMBcontrol(backend = backend)))$tmb_obj
      values <- utils::modifyList(list(b_j = b_j, ln_phi = log(phi)),
        as.list(case$par))
      par <- obj$par
      for (p in names(values)) par[names(par) == p] <- values[[p]]
      set.seed(1)
      y <- obj$simulate(par)$y_i
      for (x in c(0, 1)) {
        info <- paste(name, backend, x)
        yx <- y[d$x == x]
        eta <- b_j[[1L]] + b_j[[2L]] * x
        m <- case[[2L]](stats::make.link(family$link)$linkinv(eta))
        expect_lt(abs(mean(yx) - m[[1L]]), 5 * sd(yx) / sqrt(length(yx)),
          label = info)
        if (!is.na(m[[2L]])) {
          m4 <- mean((yx - mean(yx))^4)
          expect_lt(abs(var(yx) - m[[2L]]),
            5 * sqrt((m4 - var(yx)^2) / length(yx)), label = info)
        }
      }
    }
  }
})

test_that("RTMB truncated NB draws are finite and match conditional moments", {
  n <- 20000L
  for (family in c("truncated_nbinom1", "truncated_nbinom2")) {
    for (mu in c(1e-20, 1e-10, 0.1, 1, 10, 1e4)) {
      for (phi in list(1, rep(c(0.5, 2), n / 2))) {
        set.seed(1)
        y <- rtmb_obs_simulate(family, list(mu = rep(mu, n), phi = phi))
        info <- paste(family, mu, length(phi))
        expect_true(all(is.finite(y) & y >= 1 & y == round(y)), info = info)
        # Conditional moments given Y > 0, averaged over phi values.
        nb1 <- family == "truncated_nbinom1"
        size <- if (nb1) mu / phi else phi
        q <- exp(stats::pnbinom(0, size = size, mu = mu,
          lower.tail = FALSE, log.p = TRUE))
        v <- if (nb1) mu * (1 + phi) else mu + mu^2 / phi
        m <- rep(mu / q, length.out = n)
        s2 <- rep((v + mu^2) / q - (mu / q)^2, length.out = n)
        se <- sqrt(max(mean(s2), 0) / n) # cancellation at tiny means
        expect_lt(abs(mean(y) - mean(m)), 5 * se + 1e-8, label = info)
      }
    }
  }
  # Tiny-mean NB2 concentrates at one; NB1 at phi = 1 tends to 1 / log(2).
  y <- rtmb_obs_simulate("truncated_nbinom2",
    list(mu = rep(1e-20, 1000), phi = 1))
  expect_true(all(y == 1))
  set.seed(2)
  y <- rtmb_obs_simulate("truncated_nbinom1",
    list(mu = rep(1e-20, n), phi = 1))
  expect_equal(mean(y), 1 / log(2), tolerance = 0.03)
})
