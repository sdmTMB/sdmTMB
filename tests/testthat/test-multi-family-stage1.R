test_that("Stage 1 multi-family builds TMB data", {
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

test_that("Stage 1 multi-family can fit a simple model", {
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

test_that("Stage 1 multi-family binomial proportions use weights as size", {
  dat <- data.frame(
    y = c(0.2, 0.7, 1, 0, 0.5, 0.3),
    dist = c("binomial", "binomial", "gaussian", "gaussian", "binomial", "binomial"),
    w = c(10, 8, 2, 3, 12, 6)
  )
  fam <- list(
    gaussian = gaussian(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    weights = dat$w,
    do_fit = FALSE
  )
  prop_rows <- dat$dist == "binomial" & dat$y > 0 & dat$y < 1
  expect_equal(
    fit$tmb_data$size[prop_rows],
    dat$w[prop_rows]
  )
  expect_equal(
    fit$tmb_data$y_i[prop_rows, 1],
    dat$y[prop_rows] * dat$w[prop_rows]
  )
  expect_equal(
    fit$tmb_data$weights_i[prop_rows],
    rep(1, sum(prop_rows))
  )
  binary_rows <- dat$dist == "binomial" & dat$y %in% c(0, 1)
  expect_equal(
    fit$tmb_data$size[binary_rows],
    rep(1, sum(binary_rows))
  )
  expect_equal(
    fit$tmb_data$weights_i[binary_rows],
    dat$w[binary_rows]
  )
})
