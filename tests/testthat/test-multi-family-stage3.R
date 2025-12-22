test_that("Stage 3 multi-family supports delta families", {
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
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$delta_family_e, c(1L, 0L))
  expect_equal(
    fit$tmb_data$family_e1,
    c(sdmTMB:::.enum_family("binomial"), sdmTMB:::.enum_family("poisson"))
  )
  expect_equal(
    fit$tmb_data$family_e2,
    c(sdmTMB:::.enum_family("Gamma"), -1L)
  )
  expect_equal(ncol(fit$tmb_data$y_i), 2L)
  expect_true(all(is.na(fit$tmb_data$y_i[dat$dist == "poisson", 2])))
})
