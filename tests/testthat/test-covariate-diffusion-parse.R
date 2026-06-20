make_dl_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("covariate_diffusion parses valid terms", {
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
    covariate_diffusion = ~ space(x_num) + time(x_lag),
    do_fit = FALSE
  )

  expect_equal(
    fit$covariate_diffusion_parsed$terms$component,
    c("space", "time")
  )
  expect_equal(
    fit$covariate_diffusion_parsed$terms$variable,
    c("x_num", "x_lag")
  )
  expect_equal(
    fit$covariate_diffusion_parsed$terms$coef_name,
    c("cov_diff_space_x_num", "cov_diff_time_x_lag")
  )
})

test_that("covariate_diffusion covariates do not need to be in the main formula", {
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
    covariate_diffusion = ~ space(x_only_lag),
    do_fit = FALSE
  )

  expect_equal(fit$covariate_diffusion_parsed$covariates, "x_only_lag")
  expect_false("x_only_lag" %in% colnames(fit$tmb_data$X_ij[[1]]))
})

test_that("covariate_diffusion errors clearly for unsupported specs", {
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
      covariate_diffusion = ~ banana(x_num),
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
      covariate_diffusion = ~ space(scale(x_num)),
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
      covariate_diffusion = ~ space(not_in_data),
      do_fit = FALSE
    ),
    regexp = "Missing covariate diffusion covariate"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      covariate_diffusion = ~ space(x_chr),
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
      covariate_diffusion = ~ space(x_num) + space(x_num),
      do_fit = FALSE
    ),
    regexp = "Duplicate"
  )
})

test_that("covariate_diffusion time terms require time", {
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
      covariate_diffusion = ~ time(x_num),
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
      covariate_diffusion = ~ spacetime(x_num),
      do_fit = FALSE
    ),
    regexp = "Unsupported wrapper"
  )
})

test_that("covariate_diffusion supports single-family delta", {
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
    covariate_diffusion = ~ space(x_num),
    do_fit = FALSE
  )
  expect_true(all(c("cov_diff_space_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[1]])))
  expect_true(all(c("cov_diff_space_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[2]])))

})

test_that("covariate_diffusion requires an explicit mesh", {
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
      covariate_diffusion = ~ space(x_num),
      do_fit = FALSE
    ),
    regexp = "mesh.*covariate_diffusion"
  )
})

test_that("covariate_diffusion time terms detect irregular time spacing", {
  dat <- data.frame(
    y = rnorm(6),
    x_num = rnorm(6),
    year = rep(c(1, 3, 4), each = 2),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_dl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      covariate_diffusion = ~ time(x_num),
      extra_time = 2,
      do_fit = FALSE
    ),
    regexp = "extra_time.*not sufficient"
  )
})
