make_nl_predict_data <- function() {
  set.seed(101)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 3) + y / max(y)))
  eta <- 0.2 + 0.5 * x1 - 0.3 * x2
  data.frame(
    y = eta + rnorm(length(eta), sd = 0.1),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

make_nl_predict_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

make_nl_predict_grid <- function(mesh, years = NULL) {
  make_nl_covariate_grid(mesh, years, c("x1", "x2"))
}

make_nl_predict_delta_data <- function() {
  set.seed(202)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 3) + y / max(y)))
  positive <- exp(0.2 + 0.4 * x1)
  present <- as.integer(x2 > 0)
  data.frame(
    y = ifelse(present == 1L, positive, 0),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

test_that("covariate diffusion predict works for default and newdata pathways", {
  skip_on_cran()

  dat <- make_nl_predict_data()
  mesh <- make_nl_predict_mesh(dat)
  grid <- make_nl_predict_grid(mesh, sort(unique(dat$year)))

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0)
  )

  expect_silent(p_fit <- predict(fit))
  expect_silent(p_new <- predict(fit, newdata = dat))
  expect_equal(nrow(p_fit), nrow(dat))
  expect_equal(p_fit$est, p_new$est, tolerance = 1e-6)
  expect_equal(p_fit$est_non_rf, p_new$est_non_rf, tolerance = 1e-6)
  expect_true(all(c("nl_diffusion_x1", "nl_time_lag_x2") %in% names(p_fit)))
  expect_equal(p_fit$nl_diffusion_x1, p_new$nl_diffusion_x1, tolerance = 1e-6)
  expect_equal(p_fit$nl_time_lag_x2, p_new$nl_time_lag_x2, tolerance = 1e-6)

  p_se <- predict(fit, newdata = dat, re_form = NA, se_fit = TRUE)
  expect_true("est_se" %in% names(p_se))
  expect_true(all(is.finite(p_se$est_se)))

  expect_silent(p_pop <- predict(fit, newdata = dat, re_form = NA))
  expect_equal(nrow(p_pop), nrow(dat))

  fit_sim <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(
      start = list(log_kappaS_nl = 0),
      map = list(log_kappaS_nl = factor(NA)),
      newton_loops = 0
    )
  )
  sims <- predict(fit_sim, newdata = dat, nsim = 2)
  expect_equal(dim(sims), c(nrow(dat), 2L))
})

test_that("covariate diffusion newdata covariate changes prediction direction when lag beta is fixed positive", {
  skip_on_cran()

  dat <- make_nl_predict_data()
  mesh <- make_nl_predict_mesh(dat)
  grid <- make_nl_predict_grid(mesh, sort(unique(dat$year)))

  proto <- sdmTMB(
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

  lag_col <- proto$nonlocal_parsed$term_coef_name
  lag_idx <- match(lag_col, colnames(proto$tmb_data$X_ij[[1]]))
  b_map <- seq_along(proto$tmb_params$b_j)
  b_map[lag_idx] <- NA_integer_
  b_start <- proto$tmb_params$b_j
  b_start[lag_idx] <- 1

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
    control = sdmTMBcontrol(
      start = list(b_j = b_start),
      map = list(b_j = factor(b_map)),
      newton_loops = 0,
      getsd = FALSE
    )
  )

  grid_low <- grid
  grid_high <- grid
  grid_high$x1 <- grid_high$x1 + 0.5

  p_low <- predict(fit, newdata = dat, nonlocal_newdata = grid_low)
  p_high <- predict(fit, newdata = dat, nonlocal_newdata = grid_high)

  expect_true("nl_time_lag_x1" %in% names(p_low))
  expect_gt(mean(p_high$nl_time_lag_x1 - p_low$nl_time_lag_x1), 0)
  expect_gt(mean(p_high$est - p_low$est), 0)
})

test_that("delta covariate diffusion in component 2 changes combined response predictions", {
  skip_on_cran()

  dat <- make_nl_predict_delta_data()
  mesh <- make_nl_predict_mesh(dat)
  grid <- make_nl_predict_grid(mesh, sort(unique(dat$year)))

  proto <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  lag_col <- proto$nonlocal_parsed$term_coef_name
  lag_idx1 <- match(lag_col, colnames(proto$tmb_data$X_ij[[1]]))
  lag_idx2 <- match(lag_col, colnames(proto$tmb_data$X_ij[[2]]))

  b_map1 <- seq_along(proto$tmb_params$b_j)
  b_map1[lag_idx1] <- NA_integer_
  b_start1 <- proto$tmb_params$b_j
  b_start1[lag_idx1] <- 0

  b_map2 <- seq_along(proto$tmb_params$b_j2)
  b_map2[lag_idx2] <- NA_integer_
  b_start2 <- proto$tmb_params$b_j2
  b_start2[lag_idx2] <- 1

  fit <- suppressWarnings(sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    nonlocal_formula = ~ time_lag(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(
      start = list(b_j = b_start1, b_j2 = b_start2),
      map = list(b_j = factor(b_map1), b_j2 = factor(b_map2)),
      newton_loops = 0,
      getsd = FALSE
    )
  ))

  grid_low <- grid
  grid_high <- grid
  grid_high$x1 <- grid_high$x1 + 0.5

  p_low <- predict(fit, newdata = dat, nonlocal_newdata = grid_low, type = "response")
  p_high <- predict(fit, newdata = dat, nonlocal_newdata = grid_high, type = "response")

  expect_gt(mean(p_high$est - p_low$est), 0)
})

test_that("time-indexed covariate diffusion full-time guard is scoped to no-grid prediction", {
  object <- list(
    nonlocal_formula_parsed = list(time_indexed = TRUE),
    nonlocal_grid_supplied = FALSE
  )
  expect_true(.nonlocal_prediction_requires_full_time(object))

  object$nonlocal_grid_supplied <- TRUE
  expect_false(.nonlocal_prediction_requires_full_time(object))
  expect_false(.nonlocal_prediction_requires_full_time(object, nonlocal_newdata = data.frame()))

  object$nonlocal_grid_supplied <- FALSE
  object$nonlocal_formula_parsed$time_indexed <- FALSE
  expect_false(.nonlocal_prediction_requires_full_time(object))
})

test_that("covariate diffusion prediction errors for zero-support newdata", {
  skip_on_cran()

  dat <- make_nl_predict_data()
  mesh <- make_nl_predict_mesh(dat)
  grid <- make_nl_predict_grid(mesh, sort(unique(dat$year)))
  fit_grid <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  grid_sparse <- grid[grid$X == grid$X[1], , drop = FALSE]
  expect_error(
    predict(fit_grid, newdata = dat, nonlocal_newdata = grid_sparse),
    regexp = "zero mesh-vertex support"
  )
})

test_that("space-only covariate diffusion predict works without modeled time", {
  skip_on_cran()

  dat <- make_nl_predict_data()
  dat$year <- NULL
  mesh <- make_nl_predict_mesh(dat)
  grid <- make_nl_predict_grid(mesh)

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  p <- predict(fit)
  expect_equal(nrow(p), nrow(dat))
  expect_true("nl_diffusion_x1" %in% names(p))
})
