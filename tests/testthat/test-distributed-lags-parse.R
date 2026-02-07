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
    distributed_lags = ~ space(x_num) + time(x_lag) + spacetime(x_num),
    do_fit = FALSE
  )

  expect_equal(
    fit$distributed_lags_parsed$terms$component,
    c("space", "time", "spacetime")
  )
  expect_equal(
    fit$distributed_lags_parsed$terms$variable,
    c("x_num", "x_lag", "x_num")
  )
  expect_equal(
    fit$distributed_lags_parsed$terms$coef_name,
    c("dl_space_x_num", "dl_time_x_lag", "dl_spacetime_x_num")
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
    distributed_lags = ~ space(x_only_lag),
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
      distributed_lags = ~ space(scale(x_num)),
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
      distributed_lags = ~ space(not_in_data),
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
      distributed_lags = ~ space(x_chr),
      do_fit = FALSE
    ),
    regexp = "must be numeric"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      distributed_lags = ~ space(x_num) + space(x_num),
      do_fit = FALSE
    ),
    regexp = "Duplicate"
  )
})

test_that("distributed_lags time terms require time", {
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
      distributed_lags = ~ time(x_num),
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
      distributed_lags = ~ spacetime(x_num),
      do_fit = FALSE
    ),
    regexp = "require a `time` argument"
  )
})

test_that("distributed_lags supports single-family delta and rejects multi-family", {
  dat <- data.frame(
    y = c(0, 1, 0, 2, 0.5, 1.2),
    x_num = rnorm(6),
    dist = c("gaussian", "poisson", "gaussian", "poisson", "gaussian", "poisson"),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_dl_mesh(dat)

  fit_delta <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    distributed_lags = ~ space(x_num),
    do_fit = FALSE
  )
  expect_true(all(c("dl_space_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[1]])))
  expect_true(all(c("dl_space_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[2]])))

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      family = list(gaussian = gaussian(), poisson = poisson()),
      distribution_column = "dist",
      distributed_lags = ~ space(x_num),
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
      distributed_lags = ~ space(x_num),
      do_fit = FALSE
    ),
    regexp = "mesh.*distributed_lags"
  )
})
