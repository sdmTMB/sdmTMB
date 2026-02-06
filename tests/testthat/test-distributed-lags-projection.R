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
    x1 = c(2, 4, 6, NA),
    x2 = c(1, 3, 5, 7)
  )

  out <- sdmTMB:::.build_vertex_time_covariates(
    covariate_data = dat,
    covariates = c("x1", "x2"),
    A_st = A_st,
    year_i = c(0L, 0L, 1L, 1L),
    A_spatial_index = c(0L, 1L, 2L, 1L),
    n_t = 2L
  )

  expect_equal(dim(out$covariate_vertex_time), c(2L, 2L, 2L))

  expect_equal(
    out$covariate_vertex_time[, , 1],
    matrix(c(
      8 / 3, 0,
      4, 6
    ), nrow = 2, byrow = TRUE),
    tolerance = 1e-8
  )

  expect_equal(
    out$covariate_vertex_time[, , 2],
    matrix(c(
      5 / 3, 7,
      3, 17 / 3
    ), nrow = 2, byrow = TRUE),
    tolerance = 1e-8
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

  out1 <- sdmTMB:::.build_vertex_time_covariates(
    covariate_data = dat,
    covariates = "x1",
    A_st = A_st,
    year_i = year_i,
    A_spatial_index = spatial_i,
    n_t = 2L
  )

  dat2 <- dat
  dat2$x1[year_i == 1L] <- dat2$x1[year_i == 1L] + 100

  out2 <- sdmTMB:::.build_vertex_time_covariates(
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
    sdmTMB:::.build_vertex_time_covariates(
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
    sdmTMB:::.build_vertex_time_covariates(
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

test_that(".build_distributed_lag_tmb_data returns term-covariate mapping", {
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

  parsed <- sdmTMB:::.parse_distributed_lags_formula(
    ~ spatial(x1) + temporal(x2) + spatiotemporal(x1)
  )
  parsed <- sdmTMB:::.validate_distributed_lag_terms(
    parsed,
    data = dat,
    time = "year",
    delta = FALSE,
    multi_family = FALSE
  )

  out <- sdmTMB:::.build_distributed_lag_tmb_data(
    distributed_lags = parsed,
    data = dat,
    A_st = A_st,
    A_spatial_index = c(0L, 1L, 2L, 1L),
    year_i = c(0L, 0L, 1L, 1L),
    n_t = 2L
  )

  expect_equal(out$term_covariate_index, c(1L, 2L, 1L))
  expect_equal(out$term_component_id, c(1L, 2L, 3L))
  expect_equal(dim(out$covariate_vertex_time), c(2L, 2L, 2L))
})

test_that("sdmTMB builds distributed_lags_data in fit path", {
  dat <- data.frame(
    y = rnorm(8),
    x1 = rnorm(8),
    x2 = rnorm(8),
    year = rep(1:4, each = 2),
    X = rep(1:4, each = 2),
    Y = rep(c(0, 1), 4)
  )

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)

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

  expect_true(!is.null(fit$distributed_lags_data))
  expect_equal(fit$distributed_lags_data$n_terms, 2L)
  expect_equal(fit$distributed_lags_data$n_covariates, 2L)
  expect_equal(
    dim(fit$distributed_lags_data$covariate_vertex_time),
    c(ncol(fit$tmb_data$A_st), fit$tmb_data$n_t, 2L)
  )
})
