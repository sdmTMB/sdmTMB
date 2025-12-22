test_that("Stage 2 multi-family supports extra-parameter families", {
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

test_that("Stage 2 multi-family supports Beta and betabinomial", {
  dat <- data.frame(
    y = c(0.2, 0.5, 0.8, 0.1, 0.4, 0.7),
    dist = c(
      "Beta", "Beta", "Beta",
      "betabinomial", "betabinomial", "betabinomial"
    ),
    w = c(1, 1, 1, 10, 12, 8)
  )
  fam <- list(
    Beta = Beta(),
    betabinomial = betabinomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    weights = dat$w,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$ln_phi_len, c(1L, 1L))
  expect_equal(fit$tmb_data$ln_phi_start, c(0L, 1L))
  expect_equal(
    fit$tmb_data$size[dat$dist == "betabinomial"],
    dat$w[dat$dist == "betabinomial"]
  )
  expect_true(all(fit$tmb_data$weights_i[dat$dist == "betabinomial"] == 1))
  expect_equal(
    fit$tmb_data$y_i[dat$dist == "betabinomial", 1],
    dat$y[dat$dist == "betabinomial"] * dat$w[dat$dist == "betabinomial"]
  )
})
