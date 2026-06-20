make_dl_grid_data <- function() {
  set.seed(303)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  eta <- 0.2 + 0.5 * x1
  data.frame(
    y = eta + rnorm(length(eta), sd = 0.1),
    x1 = x1,
    year = year,
    X = x,
    Y = y
  )
}

make_dl_grid_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

make_dl_dense_grid <- function(years) {
  set.seed(404)
  n_grid <- 40L
  X <- runif(n_grid, 1, 6)
  Y <- runif(n_grid, 1, 2)
  grid <- expand.grid(grid_id = seq_len(n_grid), year = years)
  grid$X <- X[grid$grid_id]
  grid$Y <- Y[grid$grid_id]
  grid$x1 <- as.numeric(scale(sin(grid$year / 2) + grid$X / max(grid$X)))
  grid$grid_id <- NULL
  grid
}

test_that("covariate_diffusion_data can supply missing temporal slices", {
  skip_on_cran()
  skip_on_ci()

  dat <- make_dl_grid_data()
  mesh <- make_dl_grid_mesh(dat)
  dat_irregular <- dat[dat$year != 3, , drop = FALSE]
  grid <- make_dl_dense_grid(sort(unique(dat$year)))

  expect_error(
    sdmTMB(
      y ~ x1,
      data = dat_irregular,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      family = gaussian(),
      covariate_diffusion = ~ time(x1),
      extra_time = 3,
      control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
    ),
    regexp = "irregular"
  )

  fit_grid <- sdmTMB(
    y ~ x1,
    data = dat_irregular,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ time(x1),
    covariate_diffusion_data = grid,
    extra_time = 3,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  expect_true(isTRUE(fit_grid$covariate_diffusion_grid_supplied))
  expect_equal(fit_grid$covariate_diffusion_data$n_t, length(unique(dat$year)))
})

test_that("covariate_diffusion_data validates the new grid requirements", {
  skip_on_cran()
  skip_on_ci()

  dat <- make_dl_grid_data()
  mesh <- make_dl_grid_mesh(dat)
  grid <- make_dl_dense_grid(sort(unique(dat$year)))

  fit_args <- function(grid_data) {
    sdmTMB(
      y ~ x1,
      data = dat,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      family = gaussian(),
      covariate_diffusion = ~ space(x1) + time(x1),
      covariate_diffusion_data = grid_data,
      control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
    )
  }

  expect_error(
    fit_args(subset(grid, year != min(year))),
    regexp = "does not cover all"
  )

  grid_bad_xy <- grid
  grid_bad_xy$X[1] <- Inf
  expect_error(fit_args(grid_bad_xy), regexp = "coordinates must be finite")

  expect_error(
    sdmTMB(
      y ~ x1,
      data = dat,
      mesh = mesh,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      family = gaussian(),
      covariate_diffusion_data = grid,
      control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
    ),
    regexp = "covariate_diffusion_data.*was supplied"
  )
})

test_that("predict() reuses or overrides the stored covariate diffusion grid", {
  skip_on_cran()
  skip_on_ci()

  dat <- make_dl_grid_data()
  mesh <- make_dl_grid_mesh(dat)
  grid <- make_dl_dense_grid(sort(unique(dat$year)))

  fit_grid <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ space(x1) + time(x1),
    covariate_diffusion_data = grid,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  nd_full <- dat
  nd_xy_time_only <- dat[, c("X", "Y", "year"), drop = FALSE]
  nd_subset <- subset(nd_xy_time_only, year == min(year))

  p_full <- predict(fit_grid, newdata = nd_full)
  p_reused <- predict(fit_grid, newdata = nd_xy_time_only)
  p_subset <- predict(fit_grid, newdata = nd_subset)

  expect_equal(p_full$est, p_reused$est, tolerance = 1e-6)
  expect_equal(nrow(p_subset), nrow(nd_subset))

  grid_counterfactual <- grid
  grid_counterfactual$x1 <- grid_counterfactual$x1 + 5
  p_override <- predict(
    fit_grid,
    newdata = nd_xy_time_only,
    covariate_diffusion_data = grid_counterfactual
  )

  expect_false(isTRUE(all.equal(
    p_reused$diffusion_cov_space_x1,
    p_override$diffusion_cov_space_x1
  )))
})
