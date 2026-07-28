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
