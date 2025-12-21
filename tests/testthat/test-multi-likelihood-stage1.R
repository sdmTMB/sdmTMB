test_that("Stage 1 multi-likelihood builds TMB data", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
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
    do_fit = FALSE
  )
  expect_identical(fit$tmb_data$multi_family, 1L)
  expect_equal(fit$tmb_data$e_i, c(0L, 1L, 2L, 0L, 1L, 2L))
  expect_equal(
    fit$tmb_data$family_e1,
    c(
      sdmTMB:::.enum_family("gaussian"),
      sdmTMB:::.enum_family("poisson"),
      sdmTMB:::.enum_family("binomial")
    )
  )
  expect_equal(
    fit$tmb_data$link_e1,
    c(
      sdmTMB:::.enum_link("identity"),
      sdmTMB:::.enum_link("log"),
      sdmTMB:::.enum_link("logit")
    )
  )
  expect_equal(fit$tmb_data$family_e2, rep(-1L, 3))
  expect_equal(fit$tmb_data$link_e2, rep(-1L, 3))
  expect_equal(fit$tmb_data$delta_family_e, rep(0L, 3))
})

test_that("Stage 1 multi-likelihood can fit a simple model", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist"
  )
  expect_s3_class(fit, "sdmTMB")
  expect_true(is.list(fit$family))
})
