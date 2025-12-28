test_that("Multi-family family list validation works", {
  fam <- list(
    gauss = gaussian(),
    pois = poisson(link = "log")
  )
  res <- sdmTMB:::.validate_multi_family_list(fam)
  expect_identical(res$family_names, c("gauss", "pois"))
  expect_equal(res$family_enum[[1]], sdmTMB:::.enum_family("gaussian"))
  expect_equal(res$link_enum[[2]], sdmTMB:::.enum_link("log"))
})

test_that("Multi-family family list validation rejects mix families", {
  fam <- list(mix = gamma_mix())
  expect_error(
    sdmTMB:::.validate_multi_family_list(fam),
    regexp = "_mix"
  )
})

test_that("Multi-family family list validation requires names", {
  fam <- list(gaussian(), poisson())
  expect_error(
    sdmTMB:::.validate_multi_family_list(fam),
    regexp = "named list"
  )
})

test_that("Multi-family distribution_column maps to e_i", {
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  res <- sdmTMB:::.validate_multi_family_list(fam, data = dat, distribution_column = "dist")
  expect_equal(res$e_i, c(1L, 2L, 3L, 1L, 2L, 3L))
})

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
  pred_resp <- predict(fit, newdata = dat, type = "response")
  pred_link <- predict(fit, newdata = dat, type = "link")
  expected <- pred_link$est
  expected[dat$dist == "poisson"] <- exp(expected[dat$dist == "poisson"])
  expect_equal(pred_resp$est, expected, tolerance = 1e-6)
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
  pred <- predict(fit, newdata = dat, type = "link")
  expect_true(all(c("est1", "est2") %in% names(pred)))
  expect_true(all(is.na(pred$est2[dat$dist == "poisson"])))
  expect_true(all(!is.na(pred$est2[dat$dist == "delta_gamma"])))
  idx <- dat$dist == "delta_gamma"
  expect_equal(
    pred$est[idx],
    log(plogis(pred$est1[idx]) * exp(pred$est2[idx])),
    tolerance = 1e-6
  )
})
