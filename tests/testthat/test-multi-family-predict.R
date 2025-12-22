test_that("Multi-family prediction requires distribution_column and maps e_g", {
  dat <- data.frame(
    y = c(1, 2, 0, 3),
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
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_error(
    predict(fit, newdata = dat[, "y", drop = FALSE]),
    regexp = "distribution_column"
  )
  nd_bad <- dat
  nd_bad$dist[1] <- "bad"
  expect_error(
    predict(fit, newdata = nd_bad),
    regexp = "Unknown family names"
  )
  pred <- predict(fit, newdata = dat, return_tmb_object = TRUE)
  expect_equal(
    pred$pred_tmb_data$e_g,
    match(dat$dist, names(fam)) - 1L
  )
})

test_that("Multi-family predictions transform per-family links", {
  dat <- data.frame(
    y = c(1, 2, 0, 3),
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
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  pred <- predict(fit, newdata = dat, type = "response", return_tmb_object = TRUE)
  eta <- pred$report$proj_eta[,1]
  expected <- eta
  expected[dat$dist == "poisson"] <- exp(expected[dat$dist == "poisson"])
  expect_equal(pred$data$est, expected, tolerance = 1e-6)
})

test_that("Multi-family predictions handle delta families", {
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
  pred <- predict(fit, newdata = dat, type = "response")
  expect_true(all(c("est1", "est2") %in% names(pred)))
  expect_true(all(is.na(pred$est2[dat$dist == "poisson"])))
  expect_true(all(!is.na(pred$est2[dat$dist == "delta_gamma"])))
  idx <- dat$dist == "delta_gamma"
  expect_equal(pred$est[idx], pred$est1[idx] * pred$est2[idx], tolerance = 1e-6)
})
