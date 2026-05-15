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

test_that("plot_diffusion_kernel returns mesh-based diffusion fields", {
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

  out <- plot_diffusion_kernel(
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

test_that("plot_diffusion_kernel builds a ggplot object when ggplot2 is installed", {
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

  out <- plot_diffusion_kernel(
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

test_that("plot_diffusion_kernel plots raw colour values", {
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

  out <- plot_diffusion_kernel(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 1L,
    n_steps = 2,
    common_scale = TRUE,
    plot = FALSE
  )

  expect_equal(out$triangle_df$value_plot, out$triangle_df$value, tolerance = 1e-12)
})

test_that("plot_diffusion_kernel signed_sqrt handles signed values", {
  x <- c(-4, -1, 0, 9)
  expect_equal(.dl_plot_transform_values(x, "signed_sqrt"), c(-2, -1, 0, 3))
  expect_error(.dl_plot_transform_values(x, "sqrt"), regexp = "non-negative")
})

test_that("plot_diffusion_kernel does not require time_value for space-only covariate diffusions", {
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

  out <- plot_diffusion_kernel(
    fit,
    covariate = "x1",
    component = "space",
    n_steps = 2,
    plot = FALSE
  )

  expect_equal(out$impulse_time_index, 1L)
})

test_that("plot_diffusion_kernel defaults to the first time index for temporal covariate diffusions", {
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

  out <- plot_diffusion_kernel(
    fit,
    covariate = "x1",
    component = "time",
    n_steps = 2,
    plot = FALSE
  )

  expect_equal(out$impulse_time_index, 1L)
})

test_that("plot_diffused_covariate returns original and diffused mesh fields", {
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

  out <- plot_diffused_covariate(
    fit,
    covariate = "x1",
    component = "space",
    plot = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_true(is.matrix(out$original_vertex_time))
  expect_true(is.matrix(out$transformed_vertex_time))
  expect_equal(dim(out$original_vertex_time), dim(out$transformed_vertex_time))
  expect_equal(nrow(out$original_vertex_time), fit$spde$mesh$n)
  expect_equal(ncol(out$triangle_values), 2L)
  expect_equal(levels(out$triangle_df$panel), paste0(c("original", "diffused"), " (t=", out$time_value, ")"))
  expect_equal(out$component, "space")
  expect_equal(out$covariate, "x1")
})

test_that("plot_diffused_covariate selects requested time slices", {
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

  out <- plot_diffused_covariate(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 3L,
    plot = FALSE
  )

  expect_equal(out$time_index, 3L)
  expect_equal(out$time_value, 3L)
  expect_equal(dim(out$transformed_vertex_time), dim(out$original_vertex_time))
  expect_equal(levels(out$triangle_df$panel), c("original (t=3)", "diffused (t=3)"))
})

test_that("plot_diffused_covariate plots lagged contributions from one time slice", {
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

  out <- plot_diffused_covariate(
    fit,
    covariate = "x1",
    component = "time",
    time_value = 2L,
    n_steps = 3L,
    plot = FALSE
  )

  expect_equal(out$transformed_time_index, 2:4)
  expect_equal(out$transformed_time_values, 2:4)
  expect_equal(out$source_vertex_time[, 2L], out$original_vertex_time[, 2L])
  expect_equal(out$source_vertex_time[, -2L], matrix(0, nrow = nrow(out$source_vertex_time), ncol = 4L))
  expect_equal(ncol(out$triangle_values), 4L)
  expect_equal(
    levels(out$triangle_df$panel),
    c(
      "original (t=2)",
      "diffused (t=2)",
      "lag+1 (t=3)",
      "lag+2 (t=4)"
    )
  )
})

test_that("plot_diffused_covariate can plot combined fitted transforms", {
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
    covariate_diffusion = ~ space(x1) + time(x1),
    control = sdmTMBcontrol(
      start = list(kappaT_dl_raw = 0.25)
    ),
    do_fit = FALSE
  )

  out_combined <- plot_diffused_covariate(
    fit,
    covariate = "x1",
    component = "combined",
    time_value = 2L,
    n_steps = 2L,
    plot = FALSE
  )
  out_space <- plot_diffused_covariate(
    fit,
    covariate = "x1",
    component = "space",
    time_value = 2L,
    plot = FALSE
  )

  expect_equal(out_combined$component, "combined")
  expect_equal(out_combined$transformed_time_index, 2:3)
  expect_false(isTRUE(all.equal(
    out_combined$transformed_vertex_time[, 3L],
    out_space$transformed_vertex_time[, 3L],
    tolerance = 1e-8
  )))
})

test_that("plot_diffused_covariate errors cleanly for invalid inputs", {
  skip_if_not_installed("ggplot2")

  dat <- make_dl_plot_data()
  dat$x2 <- rev(dat$x1)
  mesh <- make_dl_plot_mesh(dat)

  fit_multi <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ space(x1) + time(x2),
    do_fit = FALSE
  )

  expect_error(
    plot_diffused_covariate(fit_multi, component = "space", plot = FALSE),
    regexp = "Multiple covariate-diffusion covariates"
  )
  expect_error(
    plot_diffused_covariate(fit_multi, covariate = "x1", component = "time", plot = FALSE),
    regexp = "No term `time\\(x1\\)`"
  )
  expect_error(
    plot_diffused_covariate(fit_multi, covariate = "x1", component = "space", time_value = 99L, plot = FALSE),
    regexp = "Could not match `time_value`"
  )
})
