test_that("Multi-likelihood residuals are blocked", {
  dat <- data.frame(
    y = c(1, 0, 2, 0),
    dist = c("poisson", "binomial", "poisson", "binomial")
  )
  fam <- list(
    poisson = poisson(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_error(
    residuals(fit),
    regexp = "multi-likelihood"
  )
})

test_that("Multi-likelihood delta simulations return combined responses", {
  dat <- data.frame(
    y = c(0, 2, 0, 3, 1, 0, 5),
    dist = c(
      "delta_gamma", "delta_gamma", "delta_gamma", "delta_gamma",
      "poisson", "poisson", "poisson"
    )
  )
  fam <- list(
    delta_gamma = delta_gamma(),
    poisson = poisson()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  sims <- simulate(fit, nsim = 3, silent = TRUE)
  expect_true(is.matrix(sims))
  expect_equal(nrow(sims), nrow(dat))
  expect_equal(ncol(sims), 3L)
})
