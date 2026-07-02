make_nl_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("nonlocal_formula parses valid terms", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    x_lag = rnorm(8),
    year = rep(1:4, each = 2),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_nl_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    nonlocal_formula = ~ diffusion(x_num) + time_lag(x_lag),
    do_fit = FALSE
  )

  expect_equal(
    fit$nonlocal_formula_parsed$terms$component,
    c("diffusion", "time_lag")
  )
  expect_equal(
    fit$nonlocal_formula_parsed$terms$variable,
    c("x_num", "x_lag")
  )
  expect_equal(
    fit$nonlocal_formula_parsed$terms$coef_name,
    c("nl_diffusion_x_num", "nl_time_lag_x_lag")
  )
})

test_that("nonlocal_formula covariates do not need to be in the main formula", {
  dat <- data.frame(
    y = rnorm(8),
    x_only_lag = rnorm(8),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_nl_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    nonlocal_formula = ~ diffusion(x_only_lag),
    do_fit = FALSE
  )

  expect_equal(fit$nonlocal_formula_parsed$covariates, "x_only_lag")
  expect_false("x_only_lag" %in% colnames(fit$tmb_data$X_ij[[1]]))
})

test_that("nonlocal_formula errors clearly for unsupported specs", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    x_chr = letters[1:8],
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_nl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      nonlocal_formula = ~ banana(x_num),
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
      nonlocal_formula = ~ diffusion(scale(x_num)),
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
      nonlocal_formula = ~ diffusion(not_in_data),
      do_fit = FALSE
    ),
    regexp = "Missing nonlocal covariate"
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      nonlocal_formula = ~ diffusion(x_chr),
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
      nonlocal_formula = ~ diffusion(x_num) + diffusion(x_num),
      do_fit = FALSE
    ),
    regexp = "Duplicate"
  )
})

test_that("nonlocal_formula time terms require time", {
  dat <- data.frame(
    y = rnorm(8),
    x_num = rnorm(8),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_nl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "off",
      nonlocal_formula = ~ time_lag(x_num),
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
      nonlocal_formula = ~ spacetime(x_num),
      do_fit = FALSE
    ),
    regexp = "Unsupported wrapper"
  )
})

test_that("nonlocal_formula supports single-family delta", {
  dat <- data.frame(
    y = c(0, 1, 0, 2, 0.5, 1.2),
    x_num = rnorm(6),
    dist = c("gaussian", "poisson", "gaussian", "poisson", "gaussian", "poisson"),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_nl_mesh(dat)

  fit_delta <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    nonlocal_formula = ~ diffusion(x_num),
    do_fit = FALSE
  )
  expect_true(all(c("nl_diffusion_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[1]])))
  expect_true(all(c("nl_diffusion_x_num") %in% colnames(fit_delta$tmb_data$X_ij[[2]])))

})

test_that("nonlocal_formula requires an explicit mesh", {
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
      nonlocal_formula = ~ diffusion(x_num),
      do_fit = FALSE
    ),
    regexp = "mesh.*nonlocal_formula"
  )
})

test_that("nonlocal_formula time terms detect irregular time spacing", {
  dat <- data.frame(
    y = rnorm(6),
    x_num = rnorm(6),
    year = rep(c(1, 3, 4), each = 2),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_nl_mesh(dat)

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      nonlocal_formula = ~ time_lag(x_num),
      extra_time = 2,
      do_fit = FALSE
    ),
    regexp = "extra_time.*not sufficient"
  )
})
