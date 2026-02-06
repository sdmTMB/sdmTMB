make_dl_plumbing_data <- function() {
  data.frame(
    y = rnorm(12),
    x1 = rnorm(12),
    x2 = rnorm(12),
    year = rep(1:4, each = 3),
    X = rep(1:4, each = 3),
    Y = rep(c(0, 1, 2), 4)
  )
}

make_dl_plumbing_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("distributed lag tmb_data includes safe defaults when feature is off", {
  dat <- make_dl_plumbing_data()
  mesh <- make_dl_plumbing_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    do_fit = FALSE
  )

  expect_equal(fit$tmb_data$distributed_lag_n_terms, 0L)
  expect_equal(fit$tmb_data$distributed_lag_n_covariates, 0L)
  expect_equal(dim(fit$tmb_data$distributed_lag_covariate_vertex_time), c(1L, 1L, 1L))
  expect_length(fit$tmb_data$distributed_lag_term_component, 0L)
  expect_length(fit$tmb_data$distributed_lag_term_covariate, 0L)

  expect_true(all(c("log_kappaS_dl", "log_kappaT_dl", "kappaST_dl_unscaled") %in% names(fit$tmb_params)))
  expect_length(fit$tmb_params$log_kappaS_dl, 0L)
  expect_length(fit$tmb_params$log_kappaT_dl, 0L)
  expect_length(fit$tmb_params$kappaST_dl_unscaled, 0L)
  expect_null(fit$tmb_map[["log_kappaS_dl", exact = TRUE]])
  expect_null(fit$tmb_map[["log_kappaT_dl", exact = TRUE]])
  expect_null(fit$tmb_map[["kappaST_dl_unscaled", exact = TRUE]])
})

test_that("distributed lag coefficient slots are appended and lag parameters are length-aware", {
  dat <- make_dl_plumbing_data()
  mesh <- make_dl_plumbing_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = ~ spatial(x1) + temporal(x2),
    do_fit = FALSE
  )

  x_mat <- fit$tmb_data$X_ij[[1]]
  expect_true(all(c("dl_spatial_x1", "dl_temporal_x2") %in% colnames(x_mat)))
  expect_equal(
    unname(colSums(abs(x_mat[, c("dl_spatial_x1", "dl_temporal_x2"), drop = FALSE]))),
    c(0, 0)
  )

  expect_equal(length(fit$tmb_params$b_j), ncol(x_mat))
  expect_equal(fit$tmb_data$distributed_lag_n_terms, 2L)
  expect_equal(fit$tmb_data$distributed_lag_n_covariates, 2L)
  expect_equal(fit$tmb_data$distributed_lag_term_component, c(0L, 1L))
  expect_equal(fit$tmb_data$distributed_lag_term_covariate, c(0L, 1L))
  expect_equal(
    dim(fit$tmb_data$distributed_lag_covariate_vertex_time),
    c(ncol(fit$tmb_data$A_st), fit$tmb_data$n_t, 2L)
  )

  expect_null(fit$tmb_map[["b_j", exact = TRUE]])
  expect_length(fit$tmb_params$log_kappaS_dl, 1L)
  expect_length(fit$tmb_params$log_kappaT_dl, 1L)
  expect_length(fit$tmb_params$kappaST_dl_unscaled, 0L)
})

test_that("distributed lag parameter lengths follow used components", {
  dat <- make_dl_plumbing_data()
  mesh <- make_dl_plumbing_mesh(dat)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = ~ spatiotemporal(x1),
    do_fit = FALSE
  )

  expect_length(fit$tmb_params$log_kappaS_dl, 1L)
  expect_length(fit$tmb_params$log_kappaT_dl, 0L)
  expect_length(fit$tmb_params$kappaST_dl_unscaled, 1L)
})

test_that("no-lag fit remains numerically identical with explicit distributed_lags = NULL", {
  skip_on_cran()
  set.seed(1)
  dat <- make_dl_plumbing_data()
  mesh <- make_dl_plumbing_mesh(dat)

  fit_base <- sdmTMB(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian()
  )

  fit_null <- sdmTMB(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = NULL
  )

  expect_equal(fit_base$model$objective, fit_null$model$objective, tolerance = 1e-8)
  expect_equal(fit_base$tmb_obj$report()$eta_i, fit_null$tmb_obj$report()$eta_i, tolerance = 1e-8)
  expect_equal(fit_base$tmb_obj$report()$mu_i, fit_null$tmb_obj$report()$mu_i, tolerance = 1e-8)
})

test_that("predict tmb_data keeps distributed lag columns aligned with b_j", {
  skip_on_cran()
  set.seed(1)
  dat <- make_dl_plumbing_data()
  mesh <- make_dl_plumbing_mesh(dat)

  fit <- sdmTMB(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ spatial(x1) + temporal(x2),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  td <- predict(fit, newdata = dat, return_tmb_data = TRUE)
  lag_cols <- fit$distributed_lags_data$term_coef_name

  expect_equal(colnames(td$proj_X_ij[[1]]), colnames(fit$tmb_data$X_ij[[1]]))
  expect_equal(ncol(td$proj_X_ij[[1]]), length(fit$tmb_params$b_j))
  expect_equal(
    unname(colSums(abs(td$proj_X_ij[[1]][, lag_cols, drop = FALSE]))),
    rep(0, length(lag_cols))
  )
})
