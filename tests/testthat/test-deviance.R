test_that("Deviance calculations are correct", {
  # Poisson
  set.seed(123)
  x <- rnorm(10)
  y <- rpois(10, exp(x))
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ 1 + x, family = poisson(), spatial = "off", data = d)
  mglm <- glm(y ~ 1 + x, family = poisson(), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Gamma
  set.seed(123)
  x <- rnorm(10)
  y <- rgamma(10, shape = 0.4, scale = exp(x) / 0.4)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = Gamma(link = "log"), spatial = "off", data = d)
  mglm <- glm(y ~ x, family = Gamma(link = "log"), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Binomial
  set.seed(123)
  x <- rnorm(10)
  y <- rbinom(10, size = 1, prob = plogis(x))
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = binomial(link = "logit"), spatial = "off", data = d)
  mglm <- glm(y ~ x, family = binomial(link = "logit"), data = d)
  r1 <- unname(residuals(mglm))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.001)

  # Binomial, size > 1
  set.seed(123)
  x <- rnorm(10)
  w <- sample(1:5, size = 10, replace = TRUE)
  y <- rbinom(10, size = w, prob = plogis(x))
  d <- data.frame(x = x, y = y, prop = y / w)
  m <- sdmTMB(prop ~ x, family = binomial(link = "logit"), weights = w, spatial = "off", data = d)
  expect_error(residuals(m, type = "deviance"), regexp = "size")
  expect_error(deviance(m), regexp = "size")

  # # NB2
  # set.seed(123)
  # x <- rnorm(20)
  # y <- rnbinom(20, size = 0.4, mu = exp(x))
  # d <- data.frame(x = x, y = y)
  # m <- sdmTMB(y ~ x, family = nbinom2(), spatial = "off", data = d)
  # mglm <- MASS::glm.nb(y ~ x, data = d)
  # r1 <- unname(residuals(mglm))
  # r2 <- residuals(m, type = "deviance")
  # expect_equal(r1, r2, tolerance = 0.01)
  # expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # NB1: glmmTMB holds the size mu / phi fixed in the saturated model rather
  # than phi, so see the fixed-phi likelihood-ratio test below

  # Tweedie
  set.seed(1)
  x <- rnorm(20)
  y <- mgcv::rTweedie(exp(x), 1.5, 1.2)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, family = tweedie(), spatial = "off", data = d)
  library(mgcv)
  mglm <- mgcv::gam(y ~ x, family = mgcv::tw(), data = d)
  r1 <- unname(residuals(mglm, type = "deviance"))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.1)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # Gaussian
  set.seed(1)
  x <- rnorm(20)
  y <- rnorm(20, x * 0.3, sd = 0.6)
  d <- data.frame(x = x, y = y)
  m <- sdmTMB(y ~ x, spatial = "off", data = d)
  mglm <- glm(y ~ x, data = d)
  r1 <- unname(residuals(mglm, type = "deviance"))
  r2 <- residuals(m, type = "deviance")
  expect_equal(r1, r2, tolerance = 0.01)
  expect_equal(deviance(mglm), deviance(m), tolerance = 0.01)

  # # lognormal
  # set.seed(1)
  # x <- rnorm(20)
  # y <- rlnorm(20, x * 0.3, sdlog = 0.2)
  # d <- data.frame(x = x, y = y)
  # m <- sdmTMB(y ~ x, family = lognormal(), spatial = "off", data = d)
  # mglm <- tinyVAST::tinyVAST(y ~ x, data = d, family = lognormal())
  # r1 <- unname(residuals(mglm, type = "deviance"))
  # r2 <- residuals(m, type = "deviance")
  # expect_equal(r1, r2, tolerance = 0.01)
  # expect_equal(sum(r1^2), deviance(m), tolerance = 0.01)

  # # delta gamma
  # m <- sdmTMB(density ~ depth_scaled, family = sdmTMB::delta_gamma(), spatial = "off", data = pcod)
  # dev1 <- deviance(m)
  # m2 <- tinyVAST::tinyVAST(density ~ depth_scaled,
  #   family = tinyVAST::delta_gamma(),
  #   data = as.data.frame(pcod), delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report(m2$obj$env$last.par.best)
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # delta lognormal
  # m <- sdmTMB(density ~ depth_scaled, family = sdmTMB::delta_lognormal(), spatial = "off", data = pcod)
  # dev1 <- deviance(m)
  # m2 <- tinyVAST::tinyVAST(density ~ depth_scaled,
  #   family = tinyVAST::delta_lognormal(),
  #   data = as.data.frame(pcod), delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report(m2$obj$env$last.par.best)
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # delta gamma poisson link
  # mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
  # m0 <- sdmTMB(
  #   density ~ depth_scaled, mesh = mesh,
  #   family = delta_gamma(type = "poisson-link"), spatial = "on", data = pcod,
  # )
  # dev1 <- deviance(m0)
  # m2 <- tinyVAST::tinyVAST(
  #   density ~ depth_scaled,
  #   family = tinyVAST::delta_gamma(type = "poisson-link"), data = as.data.frame(pcod),
  #   delta_options = list(formula = ~depth_scaled)
  # )
  # r2 <- m2$obj$report()
  # dev2 <- r2$deviance
  # expect_equal(dev1, dev2)

  # # deviance explained
  # mnull <- sdmTMB(
  #   density ~ 1, spatial = "off", mesh = mesh,
  #   family = delta_gamma(type = "poisson-link"), data = pcod,
  # )
  # (devexplained <- 1 - deviance(m0) / deviance(mnull))
  #
  # expect_equal(devexplained, m2$deviance_explained, tolerance = 0.0001)
})

test_that("Generalized gamma deviance residuals are correct", {
  set.seed(1)
  y <- rgamma(30, shape = 2, scale = 1)
  mean <- exp(rnorm(30, 0.5, 0.3))

  # Definition: sign * sqrt(sigma^2 * 2 * (saturated - fitted log density)),
  # with the saturated mean found numerically.
  numerical <- function(y, mean, sigma, Q) {
    mapply(function(y, mean) {
      sat <- optimize(function(m) rtmb_dgengamma(y, exp(m), sigma, Q),
        c(log(y) - 10, log(y) + 10), maximum = TRUE, tol = 1e-12)
      dev <- 2 * (sat$objective - rtmb_dgengamma(y, mean, sigma, Q))
      sign(sat$maximum - log(mean)) * sigma * sqrt(dev)
    }, y, mean)
  }
  for (Q in c(-0.8, 0.3, 1.5)) {
    expect_equal(rtmb_gengamma_devresid(y, mean, 0.6, Q),
      numerical(y, mean, 0.6, Q), tolerance = 1e-5, info = Q)
  }

  # Q = sigma is a gamma with shape sigma^-2: matches Gamma() and glm().
  sigma <- 0.7
  expect_equal(rtmb_gengamma_devresid(y, mean, sigma, sigma),
    sign(y - mean) * sqrt(2 * ((y - mean) / mean - log(y / mean))))
  # Q -> 0 is the lognormal: matches lognormal().
  expect_equal(rtmb_gengamma_devresid(y, mean, sigma, 1e-4),
    log(y) - (log(mean) - sigma^2 / 2), tolerance = 1e-3)

  # Fitted models: TMB and RTMB agree with the formula and with deviance().
  x <- rnorm(30)
  d <- data.frame(x = x, y = rgamma(30, shape = 2, scale = exp(0.5 * x) / 2))
  for (backend in c("tmb", "rtmb")) {
    m <- sdmTMB(y ~ x, family = gengamma(), spatial = "off", data = d,
      control = sdmTMBcontrol(backend = backend))
    p <- as.list(m$sd_report, "Estimate", report = TRUE)
    Q <- as.list(m$sd_report, "Estimate")$gengamma_Q
    r <- residuals(m, type = "deviance")
    expect_equal(r, rtmb_gengamma_devresid(d$y,
      exp(predict(m)$est), c(p$phi), Q), tolerance = 1e-6, info = backend)
    expect_equal(deviance(m), sum(r^2), info = backend)
  }
})

test_that("Censored Poisson deviance residuals are correct", {
  # Interval rows: the closed-form saturated lambda maximizes P(L <= Y <= U).
  L <- c(1, 3, 2, 10)
  U <- c(4, 3 + 70, 2, 30)
  for (j in seq_along(L)[L < U]) {
    sat <- optimize(function(l) rtmb_censpois_logprob_value(exp(l), L[j], U[j]),
      c(-5, 10), maximum = TRUE, tol = 1e-10)
    closed <- (lgamma(U[j] + 1) - lgamma(L[j])) / (U[j] - L[j] + 1)
    expect_equal(sat$maximum, closed, tolerance = 1e-4, info = j)
  }

  set.seed(1)
  x <- rnorm(60)
  d <- data.frame(x = x, y = rpois(60, exp(1 + 0.5 * x)))
  upr <- d$y
  upr[1:15] <- NA # right censored
  upr[16:30] <- d$y[16:30] + 3 # interval censored
  d$y[31:35] <- 0 # interval starting at zero
  upr[31:35] <- 2
  mpois <- glm(y ~ x, family = poisson(), data = d)

  for (backend in c("tmb", "rtmb")) {
    control <- sdmTMBcontrol(backend = backend)
    # Without censoring, it is the Poisson deviance.
    m <- sdmTMB(y ~ x, family = censored_poisson(), spatial = "off",
      data = d, censored_upper = d$y, control = control)
    expect_equal(residuals(m, type = "deviance"), unname(residuals(mpois)),
      tolerance = 1e-4, info = backend)
    expect_equal(deviance(m), deviance(mpois), tolerance = 1e-4,
      info = backend)

    # With censoring, matches 2 * (saturated - fitted) log likelihood, with
    # the saturated likelihood maximized numerically.
    m <- sdmTMB(y ~ x, family = censored_poisson(), spatial = "off",
      data = d, censored_upper = upr, control = control)
    lambda <- exp(predict(m)$est)
    ll <- function(lambda) rtmb_dcenspois(d$y, lambda, upr)
    sat <- vapply(seq_len(nrow(d)), function(i) {
      optimize(function(l) rtmb_dcenspois(d$y[i], exp(l), upr[i]),
        c(-20, 20), maximum = TRUE, tol = 1e-10)$objective
    }, numeric(1))
    r <- residuals(m, type = "deviance")
    expect_equal(r^2, 2 * (sat - ll(lambda)), tolerance = 1e-5, info = backend)
    expect_true(all(r[1:15] >= 0), info = backend)
    expect_true(all(r[31:35] <= 0), info = backend)
    expect_equal(deviance(m), sum(r^2), info = backend)
  }
})

test_that("compare_deviance() holds shape parameters fixed", {
  set.seed(1)
  d <- data.frame(x = rnorm(300))
  d$y <- rnbinom(300, size = 1.5, mu = exp(1 + 0.8 * d$x))
  fit <- sdmTMB(y ~ x, data = d, spatial = "off", family = nbinom2())
  fit0 <- sdmTMB(y ~ 1, data = d, spatial = "off", family = nbinom2())
  expect_message(out <- compare_deviance(fit, fit0), "phi")
  expect_equal(attr(out, "reduced")$parlist$ln_phi, fit$parlist$ln_phi)

  # Matches a GLM with theta fixed at the full model's estimate:
  skip_if_not_installed("MASS")
  m <- MASS::glm.nb(y ~ x, data = d)
  m0 <- glm(y ~ 1, family = MASS::negative.binomial(m$theta), data = d)
  expect_equal(out$deviance, deviance(m), tolerance = 1e-4)
  expect_equal(out$deviance_reduced, deviance(m0), tolerance = 1e-4)
  expect_equal(out$deviance_explained, 1 - deviance(m) / deviance(m0),
    tolerance = 1e-4)

  # No refit needed when the deviance depends only on the mean:
  fit_p <- sdmTMB(y ~ x, data = d, spatial = "off", family = poisson())
  fit_p0 <- sdmTMB(y ~ 1, data = d, spatial = "off", family = poisson())
  expect_no_message(out <- compare_deviance(fit_p, fit_p0))
  expect_equal(out$deviance_reduced, deviance(fit_p0))

  expect_error(compare_deviance(fit_p, fit0), "same family")
  d$y2 <- d$y + 1L
  fit_y2 <- sdmTMB(y2 ~ 1, data = d, spatial = "off", family = nbinom2())
  expect_error(compare_deviance(fit, fit_y2), "same response")
})

test_that("compare_deviance() holds Tweedie p fixed", {
  fit <- sdmTMB(density ~ depth_scaled + depth_scaled2,
    data = pcod_2011, spatial = "off", family = tweedie())
  fit0 <- sdmTMB(density ~ 1,
    data = pcod_2011, spatial = "off", family = tweedie())
  expect_message(out <- compare_deviance(fit, fit0), "Tweedie p")
  reduced <- attr(out, "reduced")
  expect_equal(reduced$parlist$thetaf, fit$parlist$thetaf)
  expect_false(isTRUE(all.equal(reduced$parlist$ln_phi, fit$parlist$ln_phi)))
  expect_equal(out$deviance_reduced, deviance(reduced))
})

test_that("Deviance includes weights, NB1 holds phi fixed, and signs are kept", {
  # Deviance residuals at parameter vector `p` of an unfitted model
  devresid_at <- function(fit, p) {
    obj <- fit$tmb_obj
    if (backend_sdmTMB(fit) == "tmb") return(obj$report(p)$devresid)
    rtmb_report_values(fit$tmb_data, obj$env$parList(p),
      deviance = TRUE)$devresid
  }
  for (backend in c("tmb", "rtmb")) {
    ctl <- sdmTMBcontrol(backend = backend)
    # Weighted Poisson
    d <- data.frame(y = c(1, 5), w = c(1, 9))
    fit <- sdmTMB(y ~ 1, data = d, weights = d$w, family = poisson(),
      spatial = "off", control = ctl)
    mu <- exp(fit$model$par[["b_j"]])
    expect_equal(deviance(fit),
      sum(2 * d$w * (d$y * log(d$y / mu) - (d$y - mu))), label = backend)

    # NB1 deviance residuals are RTMB only; the TMB template omits them
    if (backend == "tmb") {
      fit <- sdmTMB(y ~ 1, data = data.frame(y = c(0, 1, 5, 12)),
        family = nbinom1(), spatial = "off", control = ctl)
      expect_error(deviance(fit), regexp = "NB1 deviance")
    } else {
      # NB1: the saturated mean maximizes the likelihood with phi fixed
      d <- data.frame(y = c(0, 1, 5, 12, 40, 0, 3))
      fit <- sdmTMB(y ~ 1, data = d, family = nbinom1(), spatial = "off",
        control = ctl)
      mu <- exp(fit$model$par[["b_j"]])
      phi <- exp(fit$model$par[["ln_phi"]])
      ll <- function(m, y) dnbinom(y, size = m / phi, mu = m, log = TRUE)
      saturated <- vapply(d$y, function(y) {
        if (y == 0) return(0)
        optimize(ll, c(1e-8, 1e4), y = y, maximum = TRUE, tol = 1e-12)$objective
      }, numeric(1))
      expect_equal(deviance(fit), sum(2 * (saturated - ll(mu, d$y))),
        tolerance = 1e-6, label = backend)

      # y = mu still has a deviance, since the saturated mean differs from y
      fit <- sdmTMB(y ~ 1, data = data.frame(y = c(1, 1)), family = nbinom1(),
        spatial = "off", do_fit = FALSE, control = ctl)
      p <- fit$tmb_obj$par
      p[names(p) == "b_j"] <- 0
      p[names(p) == "ln_phi"] <- 0
      expect_equal(devresid_at(fit, p)[, 1]^2, rep(0.1193202, 2),
        tolerance = 1e-6, label = backend)

      # Near the Poisson limit (phi ~ 1e-8), the deviance is the Poisson one
      fit <- sdmTMB(y ~ 1, data = data.frame(y = c(1, 1, 1, 3)),
        family = nbinom1(), spatial = "off", control = ctl)
      fit_pois <- sdmTMB(y ~ 1, data = data.frame(y = c(1, 1, 1, 3)),
        family = poisson(), spatial = "off", control = ctl)
      expect_equal(deviance(fit), deviance(fit_pois), tolerance = 1e-4,
        label = backend)
    }

    # Poisson-link delta encounter residuals: negative for a zero
    fit <- sdmTMB(y ~ 1, data = data.frame(y = c(0, 2)),
      family = delta_gamma(type = "poisson-link"), spatial = "off",
      do_fit = FALSE, control = ctl)
    p <- fit$tmb_obj$par
    p[names(p) == "b_j"] <- 0 # log(1 - p) = -1
    expect_equal(devresid_at(fit, p)[, 1],
      c(-sqrt(2), sqrt(-2 * log(1 - exp(-1)))), label = backend)
  }
})
