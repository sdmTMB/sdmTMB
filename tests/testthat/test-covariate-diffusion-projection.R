test_that(".build_vertex_time_covariates matches hand-computed normalization", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )

  dat <- data.frame(
    x1 = c(2, 4, 6, 8, 10, NA),
    x2 = c(1, 3, 5, 7, 9, 11)
  )

  out <- .build_vertex_time_covariates(
    covariate_data = dat,
    covariates = c("x1", "x2"),
    A_st = A_st,
    year_i = c(0L, 0L, 0L, 1L, 1L, 1L),
    A_spatial_index = c(0L, 1L, 2L, 0L, 1L, 2L),
    n_t = 2L
  )

  expect_equal(dim(out$covariate_vertex_time), c(2L, 2L, 2L))

  expect_equal(
    out$covariate_vertex_time[, , 1],
    matrix(c(
      8 / 3, 26 / 3,
      16 / 3, 10
    ), nrow = 2, byrow = TRUE),
    tolerance = 1e-8
  )

  expect_equal(
    out$covariate_vertex_time[, , 2],
    matrix(c(
      5 / 3, 23 / 3,
      13 / 3, 31 / 3
    ), nrow = 2, byrow = TRUE),
    tolerance = 1e-8
  )
})

test_that(".build_vertex_time_covariates errors for zero-support vertices", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )

  expect_error(
    .build_vertex_time_covariates(
      covariate_data = data.frame(x1 = c(2, 4)),
      covariates = "x1",
      A_st = A_st,
      year_i = c(0L, 1L),
      A_spatial_index = c(0L, 0L),
      n_t = 2L
    ),
    regexp = "zero mesh-vertex support"
  )

  expect_error(
    .build_vertex_time_covariates(
      covariate_data = data.frame(x1 = c(2, 4, 6, NA, NA, NA)),
      covariates = "x1",
      A_st = A_st,
      year_i = c(0L, 0L, 0L, 1L, 1L, 1L),
      A_spatial_index = c(0L, 1L, 2L, 0L, 1L, 2L),
      n_t = 2L
    ),
    regexp = "only `NA` values"
  )

  sparse_grid_A <- Matrix::Matrix(rbind(c(0.5, 0.5, 0)), sparse = TRUE)
  expect_error(
    .build_vertex_time_covariates(
      covariate_data = data.frame(x1 = 1),
      covariates = "x1",
      A_st = sparse_grid_A,
      year_i = 0L,
      A_spatial_index = 0L,
      n_t = 1L
    ),
    regexp = "1 of 3 mesh vertices"
  )
})

test_that(".build_vertex_time_covariates is isolated by time slice", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )

  dat <- data.frame(x1 = c(2, 4, 6, 8))
  year_i <- c(0L, 0L, 1L, 1L)
  spatial_i <- c(0L, 1L, 2L, 1L)

  out1 <- .build_vertex_time_covariates(
    covariate_data = dat,
    covariates = "x1",
    A_st = A_st,
    year_i = year_i,
    A_spatial_index = spatial_i,
    n_t = 2L
  )

  dat2 <- dat
  dat2$x1[year_i == 1L] <- dat2$x1[year_i == 1L] + 100

  out2 <- .build_vertex_time_covariates(
    covariate_data = dat2,
    covariates = "x1",
    A_st = A_st,
    year_i = year_i,
    A_spatial_index = spatial_i,
    n_t = 2L
  )

  expect_equal(
    out1$covariate_vertex_time[, 1, 1],
    out2$covariate_vertex_time[, 1, 1],
    tolerance = 1e-10
  )
  expect_false(
    isTRUE(all.equal(
      out1$covariate_vertex_time[, 2, 1],
      out2$covariate_vertex_time[, 2, 1]
    ))
  )
})

test_that(".build_vertex_time_covariates rejects non-integer indices", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )
  dat <- data.frame(x1 = c(2, 4, 6, 8))

  expect_error(
    .build_vertex_time_covariates(
      covariate_data = dat,
      covariates = "x1",
      A_st = A_st,
      year_i = c(0, 0, 1.5, 1),
      A_spatial_index = c(0L, 1L, 2L, 1L),
      n_t = 2L
    ),
    regexp = "whole-number indices"
  )

  expect_error(
    .build_vertex_time_covariates(
      covariate_data = dat,
      covariates = "x1",
      A_st = A_st,
      year_i = c(0L, 0L, 1L, 1L),
      A_spatial_index = c(0, 1, 2.2, 1),
      n_t = 2L
    ),
    regexp = "whole-number indices"
  )
})

test_that(".build_nonlocal_tmb_data returns term-covariate mapping", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )

  dat <- data.frame(
    x1 = c(2, 4, 6, 8),
    x2 = c(1, 3, 5, 7),
    year = c(1, 1, 2, 2)
  )

  parsed <- suppressWarnings(.parse_nonlocal_formula(
    ~ diffusion(x1) + time_lag(x2)
  ))
  parsed <- suppressWarnings(.validate_nonlocal_terms(
    parsed,
    data = dat,
    time = "year",
    multi_family = FALSE
  ))

  out <- .build_nonlocal_tmb_data(
    nonlocal_formula = parsed,
    data = dat,
    A_st = A_st,
    A_spatial_index = c(0L, 1L, 2L, 1L),
    year_i = c(0L, 0L, 1L, 1L),
    n_t = 2L
  )

  expect_equal(out$term_covariate_index, c(1L, 2L))
  expect_equal(out$term_component_id, c(1L, 2L))
  expect_equal(dim(out$covariate_vertex_time), c(2L, 2L, 2L))
})

test_that(".build_nonlocal_tmb_data represents a joint operator once", {
  A_st <- Matrix::Diagonal(2L)
  dat <- data.frame(x1 = 1:4)
  parsed <- .parse_nonlocal_formula(~ diffusion(x1) + time_lag(x1))

  out <- .build_nonlocal_tmb_data(
    nonlocal_formula = parsed,
    data = dat,
    A_st = A_st,
    A_spatial_index = c(0L, 1L, 0L, 1L),
    year_i = c(0L, 0L, 1L, 1L),
    n_t = 2L
  )

  expect_equal(out$n_terms, 1L)
  expect_equal(out$term_component, "combined")
  expect_equal(out$term_component_id, 3L)
  expect_equal(out$term_coef_name, "nl_diffusion_time_lag_x1")
  expect_equal(out$covariate_has_spatial, 1L)
  expect_equal(out$covariate_has_temporal, 1L)
})

test_that("joint R solver uses the stationary space-time recursion", {
  M0 <- Matrix::Diagonal(2L, c(1, 2))
  M1 <- Matrix::Diagonal(2L, c(2, 4))
  x <- matrix(c(1, 2, 3, 4), nrow = 2L)
  kappaS <- 2
  kappaT <- 0.5

  actual <- .solve_nonlocal_vertex_time(
    "combined", x, M0, M1, kappaS, kappaT,
    has_space = TRUE, has_time = TRUE
  )
  system <- (1 + kappaT) * M0 + kappaS^(-2) * M1
  expected <- matrix(0, nrow = 2L, ncol = 2L)
  expected[, 1L] <- as.numeric(Matrix::solve(system, M0 %*% x[, 1L]))
  expected[, 2L] <- as.numeric(Matrix::solve(
    system, M0 %*% x[, 2L] + kappaT * M0 %*% expected[, 1L]
  ))

  expect_equal(actual, expected)

  # This is the stationary operator used by spacetime_lag.cpp with
  # options_z(0) = 1 and kappaST = 0:
  # P = kappaS^-2 (I_T x P_S) + kappaT ((L - I_T) x I_S),
  # where P_S = -M0^-1 M1. R matrices are column-major, so the vectorized
  # ordering is time slices of spatial vertices, as in the C++ template.
  P_s <- -as.matrix(Matrix::solve(M0, M1))
  L <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  L[row(L) == col(L) + 1L] <- 1
  P <- kappaS^(-2) * kronecker(diag(ncol(x)), P_s) +
    kappaT * kronecker(L - diag(ncol(x)), diag(nrow(x)))
  authoritative <- matrix(
    solve(diag(length(x)) - P, as.vector(x)),
    nrow = nrow(x)
  )
  expect_equal(actual, authoritative, tolerance = 1e-12)

  old_sum <- .solve_nonlocal_vertex_time(
    "diffusion", x, M0, M1, kappaS, kappaT,
    has_space = TRUE, has_time = FALSE
  ) + .solve_nonlocal_vertex_time(
    "time_lag", x, M0, M1, kappaS, kappaT,
    has_space = FALSE, has_time = TRUE
  )
  expect_gt(max(abs(actual - old_sum)), 1e-6)

  expect_error(
    .solve_nonlocal_vertex_time(
      "combined", x, M0, M1, kappaS, kappaT,
      has_space = TRUE, has_time = FALSE
    ),
    "requires both"
  )
})

test_that("fitted and prediction C++ paths use the joint operator", {
  skip_on_cran()
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    year = rep(1:4, each = 2),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
  grid <- make_nl_covariate_grid(mesh, sort(unique(dat$year)), "x")

  fit <- suppressWarnings(sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x) + time_lag(x),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  ))

  fitted_params <- fit$tmb_obj$env$parList(fit$model$par)
  kappaS <- exp(fitted_params$log_kappaS_nl[[1L]])
  kappaT <- fitted_params$kappaT_nl_raw[[1L]]
  expected_vertex_time <- .solve_nonlocal_vertex_time(
    component = "combined",
    vertex_time_input = fit$nonlocal_parsed$covariate_vertex_time[, , 1L],
    M0 = fit$tmb_data$spde$M0,
    M1 = fit$tmb_data$spde$M1,
    kappaS = kappaS,
    kappaT = kappaT,
    has_space = TRUE,
    has_time = TRUE
  )
  expected <- .project_nonlocal_vertex_time(
    transformed_vertex_time = expected_vertex_time,
    A_st = fit$tmb_data$A_st,
    A_spatial_index = fit$tmb_data$A_spatial_index,
    year_i = fit$tmb_data$year_i,
    n_t = fit$tmb_data$n_t
  )
  fit_report <- fit$tmb_obj$report()
  expect_equal(
    as.numeric(fit_report$covariate_diffusion_values[, 1L]),
    expected,
    tolerance = 1e-8
  )

  prediction_report <- predict(
    fit, newdata = dat, return_tmb_report = TRUE
  )
  expect_equal(
    as.numeric(prediction_report$proj_covariate_diffusion_values[, 1L]),
    expected,
    tolerance = 1e-6
  )
})

test_that("sdmTMB builds nonlocal_parsed in fit path", {
  dat <- data.frame(
    y = rnorm(8),
    x1 = rnorm(8),
    x2 = rnorm(8),
    year = rep(1:4, each = 2),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
  grid <- merge(
    data.frame(year = sort(unique(dat$year))),
    setNames(as.data.frame(mesh$mesh$loc[, 1:2, drop = FALSE]), c("X", "Y"))
  )
  grid$x1 <- as.numeric(scale(grid$year + grid$X))
  grid$x2 <- as.numeric(scale(grid$year + grid$Y))

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
    nonlocal_data = grid,
    do_fit = FALSE
  )

  expect_true(!is.null(fit$nonlocal_parsed))
  expect_equal(fit$nonlocal_parsed$n_terms, 2L)
  expect_equal(fit$nonlocal_parsed$n_covariates, 2L)
  expect_equal(
    dim(fit$nonlocal_parsed$covariate_vertex_time),
    c(ncol(fit$tmb_data$A_st), fit$tmb_data$n_t, 2L)
  )
})
