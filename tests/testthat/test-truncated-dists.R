test_that("Delta-truncated NB2 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))

  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = delta_truncated_nbinom2()
  )
  fit$sd_report

  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    ziformula = ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom2()
  )

  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est2, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)

  p <- predict(fit)
  p2 <- predict(fit2, type = "zlink")
  expect_equal(-p2, p$est1, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "conditional")
  expect_equal(p2, p$est2, tolerance = 1e-4)
})

test_that("Delta-truncated NB1 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))

  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = delta_truncated_nbinom1()
  )
  fit$sd_report

  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    ziformula = ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom1()
  )

  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est2, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)

  p <- predict(fit)
  p2 <- predict(fit2, type = "zlink")
  expect_equal(-p2, p$est1, tolerance = 1e-4)

  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "conditional")
  expect_equal(p2, p$est2, tolerance = 1e-4)
})

test_that("Truncated NB2 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))
  d0 <- subset(d0, y > 0)
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = truncated_nbinom2()
  )
  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom2()
  )
  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est, tolerance = 1e-4)
  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)
})

test_that("Truncated NB1 works with various types of predictions", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))
  d0 <- subset(d0, y > 0)
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    spatial = FALSE,
    family = truncated_nbinom1()
  )
  fit2 <- glmmTMB::glmmTMB(
    y ~ 1,
    data = d0,
    family = glmmTMB::truncated_nbinom1()
  )
  p <- predict(fit)
  p2 <- predict(fit2)
  expect_equal(p2, p$est, tolerance = 1e-4)
  p <- predict(fit, type = "response")
  p2 <- predict(fit2, type = "response")
  expect_equal(p2, p$est, tolerance = 1e-4)
})

test_that("Truncated NB1/2 indexes are right", {
  skip_on_cran()
  set.seed(1)
  d0 <- data.frame(y = rnbinom(n = 1e3, mu = 3, size = 5))
  d0$year <- 1

  # D-TNB2
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    time = "year",
    spatiotemporal = "off",
    spatial = FALSE,
    family = delta_truncated_nbinom2()
  )
  x <- get_index(fit, newdata = d0, area = 1/nrow(d0))
  pp <- predict(fit, type = "response")
  expect_equal(x$est, pp$est[1], tolerance = 1e-3)

  # D-TNB1
  fit <- sdmTMB(
    y ~ 1,
    data = d0,
    time = "year",
    spatiotemporal = "off",
    spatial = FALSE,
    family = delta_truncated_nbinom1()
  )
  x <- get_index(fit, newdata = d0, area = 1/nrow(d0))
  pp <- predict(fit, type = "response")
  expect_equal(x$est, pp$est[1], tolerance = 1e-3)

  # TNB2
  dd <- subset(d0, y > 0)
  fit <- sdmTMB(
    y ~ 1,
    data = dd,
    time = "year",
    spatiotemporal = "off",
    spatial = FALSE,
    family = truncated_nbinom2()
  )
  x <- get_index(fit, newdata = d0, area = 1/nrow(d0))
  pp <- predict(fit, type = "response")
  expect_equal(x$est, pp$est[1], tolerance = 1e-3)

  # TNB1
  fit <- sdmTMB(
    y ~ 1,
    data = dd,
    time = "year",
    spatiotemporal = "off",
    spatial = FALSE,
    family = truncated_nbinom2()
  )
  x <- get_index(fit, newdata = d0, area = 1/nrow(d0))
  pp <- predict(fit, type = "response")
  expect_equal(x$est, pp$est[1], tolerance = 1e-3)
})

test_that("Non-delta truncated NB families reject zeros", {
  d0 <- data.frame(y = c(0, 1, 2, 3))
  expect_error(sdmTMB(y ~ 1, data = d0, spatial = "off",
    family = truncated_nbinom1()), regexp = "response > 0")
  expect_error(sdmTMB(y ~ 1, data = d0, spatial = "off",
    family = truncated_nbinom2()), regexp = "response > 0")
})

test_that("Truncated NB response predictions use the estimated phi", {
  set.seed(1)
  d <- data.frame(y = rnbinom(300, mu = 3, size = 1.5))
  d <- d[d$y > 0, , drop = FALSE]
  fit <- sdmTMB(y ~ 1, data = d, spatial = "off", family = truncated_nbinom2())
  fit0 <- update(fit, control = sdmTMBcontrol(start = list(ln_phi = 3)))
  mu <- exp(fit$model$par[["b_j"]])
  phi <- exp(fit$model$par[["ln_phi"]])
  p <- predict(fit, type = "response")$est[1]
  expect_equal(p, mu / (1 - (1 + mu / phi)^(-phi)), tolerance = 1e-6)
  expect_equal(predict(fit0, type = "response")$est[1], p, tolerance = 1e-4)
})

test_that("Truncated NB response residuals use the truncated mean", {
  set.seed(1)
  d <- data.frame(y = rnbinom(300, mu = 3, size = 1.5))
  d <- d[d$y > 0, , drop = FALSE]
  for (family in list(truncated_nbinom1(), truncated_nbinom2())) {
    fit <- sdmTMB(y ~ 1, data = d, spatial = "off", family = family)
    expect_equal(residuals(fit, type = "response"),
      d$y - predict(fit, type = "response")$est)
  }
})

test_that("Truncated NB fits sharing a family object keep their own phi", {
  set.seed(1)
  d1 <- data.frame(y = rnbinom(300, mu = 3, size = 1.5))
  d1 <- d1[d1$y > 0, , drop = FALSE]
  d2 <- data.frame(y = rnbinom(300, mu = 3, size = 20))
  d2 <- d2[d2$y > 0, , drop = FALSE]
  fam <- truncated_nbinom2()
  fit1 <- sdmTMB(y ~ 1, data = d1, spatial = "off", family = fam)
  p1 <- predict(fit1, type = "response")$est[1]
  fit2 <- sdmTMB(y ~ 1, data = d2, spatial = "off", family = fam)
  expect_equal(predict(fit1, type = "response")$est[1], p1)
  mu <- exp(fit1$model$par[["b_j"]])
  phi <- exp(fit1$model$par[["ln_phi"]])
  expect_equal(p1, mu / (1 - (1 + mu / phi)^(-phi)), tolerance = 1e-6)

  # linkinv() needs phi supplied now that fits no longer store it
  expect_error(fam$linkinv(0), regexp = "phi")
  expect_equal(fam$linkinv(log(mu), phi = phi), p1, tolerance = 1e-6)
})
