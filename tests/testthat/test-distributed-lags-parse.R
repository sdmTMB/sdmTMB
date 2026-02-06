make_dl_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("distributed_lags parses valid terms", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    x_lag = rnorm(8),
    year = rep(1:4, each = 2),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_dl_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = ~ spatial(x_num) + temporal(x_lag) + spatiotemporal(x_num),
    do_fit = FALSE
  )

  expect_equal(
    fit$distributed_lags_parsed$terms$component,
    c("spatial", "temporal", "spatiotemporal")
  )
  expect_equal(
    fit$distributed_lags_parsed$terms$variable,
    c("x_num", "x_lag", "x_num")
  )
  expect_equal(
    fit$distributed_lags_parsed$terms$coef_name,
    c("dl_spatial_x_num", "dl_temporal_x_lag", "dl_spatiotemporal_x_num")
  )
})

test_that("distributed_lags covariates do not need to be in the main formula", {
  dat <- data.frame(
    y = rnorm(8),
    x_only_lag = rnorm(8),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_dl_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = ~ spatial(x_only_lag),
    do_fit = FALSE
  )

  expect_equal(fit$distributed_lags_parsed$covariates, "x_only_lag")
  expect_false("x_only_lag" %in% colnames(fit$tmb_data$X_ij[[1]]))
})

test_that("distributed_lags errors clearly for unsupported specs", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    x_chr = letters[1:8],
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_dl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ banana(x_num),
      do_fit = FALSE
    ),
    regexp = "Unsupported wrapper"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ spatial(scale(x_num)),
      do_fit = FALSE
    ),
    regexp = "bare variable"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ spatial(not_in_data),
      do_fit = FALSE
    ),
    regexp = "Missing distributed lag covariate"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ spatial(x_chr),
      do_fit = FALSE
    ),
    regexp = "must be numeric"
  )
})

test_that("distributed_lags temporal terms require time", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_dl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ temporal(x_num),
      do_fit = FALSE
    ),
    regexp = "require a `time` argument"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ spatiotemporal(x_num),
      do_fit = FALSE
    ),
    regexp = "require a `time` argument"
  )
})

test_that("distributed_lags rejects delta and multi-family models", {
  dat <- data.frame(
    y = c(0, 1, 0, 2, 0.5, 1.2),
    x_num = rnorm(6),
    dist = c("gaussian", "poisson", "gaussian", "poisson", "gaussian", "poisson"),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_dl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      family = delta_gamma(),
      distributed_lags = ~ spatial(x_num),
      do_fit = FALSE
    ),
    regexp = "unsupported for delta"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      family = list(gaussian = gaussian(), poisson = poisson()),
      distribution_column = "dist",
      distributed_lags = ~ spatial(x_num),
      do_fit = FALSE
    ),
    regexp = "unsupported for multi-family"
  )
})

test_that("distributed_lags requires an explicit mesh", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8)
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ spatial(x_num),
      do_fit = FALSE
    ),
    regexp = "mesh.*distributed_lags"
  )
})
