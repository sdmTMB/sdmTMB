make_nl_plot_data <- function() {
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

make_nl_plot_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

make_nl_plot_grid <- function(mesh, years = NULL, include_x2 = FALSE) {
  grid <- make_nl_covariate_grid(mesh, years, "x1")
  if (isTRUE(include_x2)) {
    grid$x2 <- rev(grid$x1)
  }
  grid
}

test_that("plot_nonlocal_kernel returns a ggplot on mesh vertices by default", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_kernel(
    fit,
    component = "diffusion",
    time_value = 1L,
    n_steps = 2
  )

  expect_s3_class(out, "ggplot")
  expect_s3_class(out$layers[[1L]]$geom, "GeomSegment")
  expect_s3_class(out$layers[[2L]]$geom, "GeomPoint")
  original <- out$data[out$data$panel == levels(out$data$panel)[[1L]], ]
  expect_equal(nrow(original), fit$spde$mesh$n)
  expect_equal(
    unname(as.matrix(original[, c("x", "y")])),
    unname(.nl_plot_extract_mesh(fit$spde$mesh)$loc)
  )
  expect_equal(sum(original$value == 1), 1L)
  expect_equal(sum(original$value == 0), fit$spde$mesh$n - 1L)
})

test_that("plot_nonlocal_kernel projects to newdata locations", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  newdata <- unique(dat[dat$year %in% 1:2, c("X", "Y")])
  newdata <- newdata[seq_len(4L), , drop = FALSE]
  out <- plot_nonlocal_kernel(
    fit,
    component = "time_lag",
    newdata = newdata,
    covariate = "x1",
    time_value = 1L,
    n_steps = 2
  )

  original <- out$data[out$data$panel == levels(out$data$panel)[[1L]], ]
  expect_s3_class(out$layers[[1L]]$geom, "GeomPoint")
  expect_equal(nrow(original), nrow(newdata))
  expect_equal(
    unname(as.matrix(original[, c("x", "y")])),
    unname(as.matrix(newdata))
  )
})

test_that("plot_nonlocal_kernel can use raster output with newdata", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  grid <- expand.grid(X = seq(1, 6, length.out = 4), Y = seq(1, 2, length.out = 3))
  out <- plot_nonlocal_kernel(
    fit,
    component = "time_lag",
    newdata = grid,
    type = "raster",
    covariate = "x1",
    time_value = 1L,
    n_steps = 2
  )

  expect_s3_class(out$layers[[1L]]$geom, "GeomRaster")
  expect_equal(nrow(out$data[out$data$panel == levels(out$data$panel)[[1L]], ]), nrow(grid))
  expect_error(
    plot_nonlocal_kernel(fit, component = "time_lag", type = "raster", covariate = "x1"),
    regexp = '`type = "raster"` requires `newdata`'
  )
})

test_that("plot_nonlocal_kernel plots raw colour values with common_scale", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_kernel(
    fit,
    component = "time_lag",
    covariate = "x1",
    time_value = 1L,
    n_steps = 2,
    common_scale = TRUE
  )

  expect_equal(out$data$value_plot, out$data$value, tolerance = 1e-12)
})

test_that("plot_nonlocal_kernel does not require time_value for space-only covariate diffusions", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_kernel(
    fit,
    component = "diffusion",
    covariate = "x1",
    n_steps = 2
  )

  expect_equal(levels(out$data$panel), c("original (t=0)", "diffused (t=0)"))
})

test_that("plot_nonlocal_kernel defaults to the first time index for temporal covariate diffusions", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_kernel(
    fit,
    component = "time_lag",
    covariate = "x1",
    n_steps = 2
  )

  expect_equal(levels(out$data$panel), c("original (t=1)", "diffused (t=1)", "lag+1 (t=2)"))
})

test_that("plot_nonlocal_covariate returns a ggplot on mesh vertices by default", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_covariate(
    fit,
    component = "diffusion",
    covariate = "x1"
  )

  expect_s3_class(out, "ggplot")
  expect_s3_class(out$layers[[1L]]$geom, "GeomSegment")
  expect_s3_class(out$layers[[2L]]$geom, "GeomPoint")
  expect_equal(levels(out$data$panel), c("original (t=0)", "diffused (t=0)"))
  expect_equal(nrow(out$data[out$data$panel == "original (t=0)", ]), fit$spde$mesh$n)
})

test_that("plot_nonlocal_covariate projects to newdata locations and raster output", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  newdata <- unique(dat[dat$year == 1L, c("X", "Y")])
  point_plot <- plot_nonlocal_covariate(
    fit,
    component = "diffusion",
    newdata = newdata,
    covariate = "x1"
  )
  raster_plot <- plot_nonlocal_covariate(
    fit,
    component = "diffusion",
    newdata = newdata,
    type = "raster",
    covariate = "x1"
  )

  expect_s3_class(point_plot$layers[[1L]]$geom, "GeomPoint")
  expect_s3_class(raster_plot$layers[[1L]]$geom, "GeomRaster")
  expect_equal(nrow(point_plot$data[point_plot$data$panel == levels(point_plot$data$panel)[[1L]], ]), nrow(newdata))
  expect_error(
    plot_nonlocal_covariate(fit, component = "diffusion", type = "raster", covariate = "x1"),
    regexp = '`type = "raster"` requires `newdata`'
  )
})

test_that("plot_nonlocal_covariate selects requested time slices", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_covariate(
    fit,
    component = "time_lag",
    covariate = "x1",
    time_value = 3L
  )

  expect_equal(levels(out$data$panel), c("original (t=3)", "diffused (t=3)"))
})

test_that("plot_nonlocal_covariate plots lagged contributions from one time slice", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  out <- plot_nonlocal_covariate(
    fit,
    component = "time_lag",
    covariate = "x1",
    time_value = 2L,
    n_steps = 3L
  )

  expect_equal(
    levels(out$data$panel),
    c(
      "original (t=2)",
      "diffused (t=2)",
      "lag+1 (t=3)",
      "lag+2 (t=4)"
    )
  )
})

test_that("plot_nonlocal_covariate can plot combined fitted transforms", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(
      start = list(kappaT_nl_raw = 0.25)
    ),
    do_fit = FALSE
  )

  out_combined <- plot_nonlocal_covariate(
    fit,
    component = "combined",
    covariate = "x1",
    time_value = 2L,
    n_steps = 2L
  )
  out_space <- plot_nonlocal_covariate(
    fit,
    component = "diffusion",
    covariate = "x1",
    time_value = 2L,
    n_steps = 2L
  )

  expect_equal(levels(out_combined$data$panel), c("original (t=2)", "diffused (t=2)", "lag+1 (t=3)"))
  combined_lag <- out_combined$data$value[out_combined$data$panel == "lag+1 (t=3)"]
  space_lag <- out_space$data$value[out_space$data$panel == "lag+1 (t=3)"]
  expect_false(isTRUE(all.equal(
    combined_lag,
    space_lag,
    tolerance = 1e-8
  )))
})

test_that("plot_nonlocal_covariate errors cleanly for invalid inputs", {
  skip_if_not_installed("ggplot2")

  dat <- make_nl_plot_data()
  dat$x2 <- rev(dat$x1)
  mesh <- make_nl_plot_mesh(dat)
  grid <- make_nl_plot_grid(mesh, sort(unique(dat$year)), include_x2 = TRUE)

  fit_multi <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  expect_error(
    plot_nonlocal_covariate(fit_multi, component = "diffusion"),
    regexp = "Multiple covariate-diffusion covariates"
  )
  expect_error(
    plot_nonlocal_covariate(fit_multi, covariate = "x1", component = "time_lag"),
    regexp = "No term `time_lag\\(x1\\)`"
  )
  expect_error(
    plot_nonlocal_covariate(fit_multi, covariate = "x1", component = "diffusion", time_value = 99L),
    regexp = "Could not match `time_value`"
  )
})
