test_that("tidy works for multi-family models", {
  dat <- data.frame(
    y = c(1, 2.5, 1, 3.2),
    dist = c("poisson", "gaussian", "poisson", "gaussian")
  )
  fam <- list(
    poisson = poisson(),
    gaussian = gaussian()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0)
  )
  td <- tidy(fit, effects = "fixed")
  expect_true("(Intercept)" %in% td$term)

  td_re <- tidy(fit, effects = "ran_pars")
  expect_true("phi_gaussian" %in% td_re$term)
})
