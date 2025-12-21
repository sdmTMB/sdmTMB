test_that("Stage 2 multi-likelihood supports extra-parameter families", {
  dat <- data.frame(
    y = c(0, 2, 4, 0.5, 1.1, 2.2, 0, 3.5),
    dist = c(
      "nbinom2", "nbinom2", "nbinom2",
      "Gamma", "Gamma", "Gamma",
      "tweedie", "tweedie"
    )
  )
  fam <- list(
    nbinom2 = nbinom2(),
    Gamma = Gamma(link = "log"),
    tweedie = tweedie(link = "log")
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
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$ln_phi_len, c(1L, 1L, 1L))
  expect_equal(fit$tmb_data$ln_phi_start, c(0L, 1L, 2L))
  expect_equal(fit$tmb_data$thetaf_len, c(0L, 0L, 1L))
  expect_equal(fit$tmb_data$thetaf_start, c(-1L, -1L, 0L))
})
