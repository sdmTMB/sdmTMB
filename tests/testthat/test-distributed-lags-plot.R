make_dl_plot_data <- function() {
  set.seed(77)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  data.frame(
    y = 0.2 + 0.4 * x1 + rnorm(length(x1), sd = 0.1),
    x1 = x1,
    year = year,
    X = x,
    Y = y
  )
}

make_dl_plot_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("plot_distributed_lag_diffusion returns mesh-based diffusion fields", {
  dat <- make_dl_plot_data()
  mesh <- make_dl_plot_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_distributed_lag_diffusion(
    fit,
    covariate = "x1",
    component = "space",
    time_value = 1L,
    n_steps = 2,
    plot = FALSE
  )

  expect_true(is.matrix(out$impulse_vertex_time))
  expect_true(is.matrix(out$transformed_vertex_time))
  expect_equal(nrow(out$impulse_vertex_time), fit$spde$mesh$n)
  expect_equal(nrow(out$transformed_vertex_time), fit$spde$mesh$n)
  expect_equal(ncol(out$triangle_values), 2L)
  expect_equal(sum(out$impulse_vertex_time), 1)
  expect_equal(out$component, "space")
  expect_equal(out$covariate, "x1")
})

test_that("plot_distributed_lag_diffusion builds a ggplot object when ggplot2 is installed", {
  skip_if_not_installed("ggplot2")

  dat <- make_dl_plot_data()
  mesh <- make_dl_plot_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ time(x1),
    do_fit = FALSE
  )

  out <- plot_distributed_lag_diffusion(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    n_steps = 2,
    plot = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(is.data.frame(out$triangle_df))
})

test_that("plot_distributed_lag_diffusion does not require time_value for space-only distributed lags", {
  skip_if_not_installed("ggplot2")

  dat <- make_dl_plot_data()
  mesh <- make_dl_plot_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_distributed_lag_diffusion(
    fit,
    covariate = "x1",
    component = "space",
    n_steps = 2,
    plot = FALSE
  )

  expect_equal(out$impulse_time_index, 1L)
})

test_that("plot_distributed_lag_diffusion requires time_value for temporal distributed lags", {
  skip_if_not_installed("ggplot2")

  dat <- make_dl_plot_data()
  mesh <- make_dl_plot_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ time(x1),
    do_fit = FALSE
  )

  expect_error(
    plot_distributed_lag_diffusion(
      fit,
      covariate = "x1",
      component = "time",
      n_steps = 2,
      plot = FALSE
    ),
    regexp = "time_value.*required"
  )
})
