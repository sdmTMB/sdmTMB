make_dl_fit_data <- function() {
  set.seed(42)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 2) + y / max(y)))
  eta <- 0.3 + 0.4 * x1 - 0.2 * x2
  data.frame(
    y = eta + rnorm(length(eta), sd = 0.15),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

make_dl_fit_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("distributed lag fits run for each wrapper and combined terms", {
  skip_on_cran()
  dat <- make_dl_fit_data()
  mesh <- make_dl_fit_mesh(dat)

  lag_forms <- list(
    spatial = ~ space(x1),
    temporal = ~ time(x1),
    spatiotemporal = ~ spacetime(x2),
    combined = ~ space(x1) + time(x1) + spacetime(x2)
  )

  for (nm in names(lag_forms)) {
    fit <- sdmTMB(
      y ~ x1 + x2,
      data = dat,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      family = gaussian(),
      distributed_lags = lag_forms[[nm]],
      control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
    )
    expect_true(is.finite(fit$model$objective), info = nm)
  }
})

test_that("distributed lag model matches no-lag model when lag coefficients are fixed at 0", {
  skip_on_cran()
  dat <- make_dl_fit_data()
  mesh <- make_dl_fit_mesh(dat)

  fit_base <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  dl_formula <- ~ space(x1) + time(x2)
  proto <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = dl_formula,
    do_fit = FALSE
  )

  lag_cols <- proto$distributed_lags_data$term_coef_name
  lag_idx <- match(lag_cols, colnames(proto$tmb_data$X_ij[[1]]))
  b_map <- seq_along(proto$tmb_params$b_j)
  b_map[lag_idx] <- NA_integer_
  b_start <- proto$tmb_params$b_j
  b_start[lag_idx] <- 0

  fit_lag_fixed <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = dl_formula,
    control = sdmTMBcontrol(
      start = list(b_j = b_start),
      map = list(
        b_j = factor(b_map)
      ),
      newton_loops = 0,
      getsd = FALSE
    )
  )

  expect_equal(fit_base$model$objective, fit_lag_fixed$model$objective, tolerance = 1e-6)
  expect_equal(fit_base$tmb_obj$report()$eta_i, fit_lag_fixed$tmb_obj$report()$eta_i, tolerance = 1e-6)
  expect_equal(fit_base$tmb_obj$report()$mu_i, fit_lag_fixed$tmb_obj$report()$mu_i, tolerance = 1e-6)
})
