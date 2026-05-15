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

test_that("plot_covariate_diffusion returns mesh-based diffusion fields", {
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
    covariate_diffusion = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion(
    fit,
    covariate = "x1",
    component = "space",
    time_value = 1L,
    n_steps = 2
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

test_that("plot_covariate_diffusion builds a ggplot object when ggplot2 is installed", {
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
    covariate_diffusion = ~ time(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    n_steps = 2
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(is.data.frame(out$triangle_df))
})

test_that("plot_covariate_diffusion plots raw colour values", {
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
    covariate_diffusion = ~ time(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    n_steps = 2,
    common_scale = TRUE
  )

  expect_equal(out$triangle_df$value_plot, out$triangle_df$value, tolerance = 1e-12)
})

test_that("plot_covariate_diffusion signed_sqrt handles signed values", {
  x <- c(-4, -1, 0, 9)
  expect_equal(.dl_plot_transform_values(x, "signed_sqrt"), c(-2, -1, 0, 3))
  expect_error(.dl_plot_transform_values(x, "sqrt"), regexp = "non-negative")
})

test_that("plot_covariate_diffusion does not require time_value for space-only covariate diffusions", {
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
    covariate_diffusion = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion(
    fit,
    covariate = "x1",
    component = "space",
    n_steps = 2
  )

  expect_equal(out$impulse_time_index, 1L)
})

test_that("plot_covariate_diffusion defaults to the first time index for temporal covariate diffusions", {
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
    covariate_diffusion = ~ time(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion(
    fit,
    covariate = "x1",
    component = "time",
    n_steps = 2
  )

  expect_equal(out$impulse_time_index, 1L)
})

test_that("plot_covariate_diffusion_grid returns regular-grid diffusion fields", {
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
    covariate_diffusion = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "space",
    n_steps = 2,
    grid_resolution = 20,
    plot = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(is.data.frame(out$grid_df))
  expect_equal(ncol(out$projected_grid_values), 2L)
  expect_equal(length(levels(out$grid_df$panel)), 2L)
  expect_equal(out$panel_totals[[1L]], 1, tolerance = 1e-8)
  expect_true(all(is.finite(out$panel_msd)))
})

test_that("plot_covariate_diffusion_grid respects longer-axis grid_resolution", {
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
    covariate_diffusion = ~ space(x1),
    do_fit = FALSE
  )

  out <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "space",
    n_steps = 2,
    grid_resolution = 30,
    plot = FALSE
  )
  loc <- as.matrix(fit$spde$mesh$loc[, 1:2])
  spans <- apply(loc, 2, function(z) diff(range(z)))
  expected_short <- max(1L, as.integer(round(30 * min(spans) / max(spans))))

  expect_equal(max(out$grid_nx, out$grid_ny), 30L)
  expect_equal(min(out$grid_nx, out$grid_ny), expected_short)
  expect_equal(nrow(out$grid), out$grid_nx * out$grid_ny)
})

test_that("plot_covariate_diffusion_grid supports panel and common scales", {
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
    covariate_diffusion = ~ time(x1),
    do_fit = FALSE
  )

  out_panel <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    grid_resolution = 15,
    scale = "panel",
    plot = FALSE
  )
  out_common <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    grid_resolution = 15,
    scale = "common",
    plot = FALSE
  )

  expect_equal(out_panel$scale, "panel")
  expect_equal(out_common$scale, "common")
  expect_equal(out_common$grid_df$value_plot, out_common$grid_df$value, tolerance = 1e-12)
  expect_true(all(out_panel$grid_df$value_plot >= 0 & out_panel$grid_df$value_plot <= 1))
})

test_that("plot_covariate_diffusion_grid can plot combined space-time response", {
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
    covariate_diffusion = ~ space(x1) + time(x1) + spacetime(x1),
    control = sdmTMBcontrol(
      start = list(kappaT_dl_raw = 0.25, kappaST_dl_raw = 0)
    ),
    do_fit = FALSE
  )

  out_combined <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "combined",
    time_value = 1L,
    n_steps = 2,
    grid_resolution = 15,
    plot = FALSE
  )
  out_spacetime <- plot_covariate_diffusion_grid(
    fit,
    covariate = "x1",
    component = "spacetime",
    time_value = 1L,
    n_steps = 2,
    grid_resolution = 15,
    plot = FALSE
  )

  expect_equal(out_combined$component, "combined")
  expect_equal(ncol(out_combined$projected_grid_values), 3L)
  expect_false(isTRUE(all.equal(
    out_combined$projected_grid_values[, -1L],
    out_spacetime$projected_grid_values[, -1L],
    tolerance = 1e-8
  )))
})

test_that("plot_covariate_diffusion_grid errors cleanly for invalid inputs", {
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
    covariate_diffusion = ~ time(x1),
    do_fit = FALSE
  )

  expect_error(
    plot_covariate_diffusion_grid(
      fit,
      covariate = "x1",
      component = "time",
      time_value = 1L,
      grid_resolution = 1,
      plot = FALSE
    ),
    regexp = "grid_resolution"
  )
  expect_error(
    plot_covariate_diffusion_grid(
      fit,
      covariate = "x1",
      time_value = 1L,
      plot = FALSE
    ),
    regexp = "component.*required"
  )
  expect_error(
    plot_covariate_diffusion_grid(
      fit,
      covariate = "x1",
      component = "time",
      plot = FALSE
    ),
    regexp = "time_value.*required"
  )
})
