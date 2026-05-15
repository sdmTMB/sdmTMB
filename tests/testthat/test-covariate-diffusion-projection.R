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

  out <- .build_vertex_time_covariates(
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

test_that(".build_covariate_diffusion_tmb_data returns term-covariate mapping", {
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

  parsed <- suppressWarnings(.parse_covariate_diffusion_formula(
    ~ space(x1) + time(x2) + spacetime(x1)
  ))
  parsed <- suppressWarnings(.validate_covariate_diffusion_terms(
    parsed,
    data = dat,
    time = "year",
    delta = FALSE,
    multi_family = FALSE
  ))

  out <- .build_covariate_diffusion_tmb_data(
    covariate_diffusion = parsed,
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

test_that(".build_covariate_diffusion_tmb_data accepts explicit vertex-time overrides", {
  A_st <- Matrix::Matrix(
    rbind(
      c(1, 0),
      c(0.5, 0.5),
      c(0, 1)
    ),
    sparse = TRUE
  )
  dat <- data.frame(x1 = c(2, 4, 6, 8), x2 = c(1, 3, 5, 7), year = c(1, 1, 2, 2))

  parsed_one <- .parse_covariate_diffusion_formula(~ space(x1))
  parsed_one <- .validate_covariate_diffusion_terms(parsed_one, dat, time = "year", delta = FALSE, multi_family = FALSE)
  vertex_mat <- matrix(1:4, nrow = 2L)
  out_mat <- .build_covariate_diffusion_tmb_data(
    covariate_diffusion = parsed_one,
    data = dat,
    A_st = A_st,
    A_spatial_index = c(0L, 1L, 2L, 1L),
    year_i = c(0L, 0L, 1L, 1L),
    n_t = 2L,
    covariate_vertex_time = vertex_mat
  )
  expect_equal(out_mat$covariate_vertex_time[, , 1], vertex_mat)

  parsed_two <- .parse_covariate_diffusion_formula(~ space(x1) + time(x2))
  parsed_two <- .validate_covariate_diffusion_terms(parsed_two, dat, time = "year", delta = FALSE, multi_family = FALSE)
  vertex_arr <- array(1:8, dim = c(2L, 2L, 2L))
  out_arr <- .build_covariate_diffusion_tmb_data(
    covariate_diffusion = parsed_two,
    data = dat,
    A_st = A_st,
    A_spatial_index = c(0L, 1L, 2L, 1L),
    year_i = c(0L, 0L, 1L, 1L),
    n_t = 2L,
    covariate_vertex_time = vertex_arr
  )
  expect_equal(unname(out_arr$covariate_vertex_time), vertex_arr)
})

test_that("sdmTMB builds covariate_diffusion_data in fit path", {
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
    covariate_diffusion = ~ space(x1) + time(x2),
    do_fit = FALSE
  )

  expect_true(!is.null(fit$covariate_diffusion_data))
  expect_equal(fit$covariate_diffusion_data$n_terms, 2L)
  expect_equal(fit$covariate_diffusion_data$n_covariates, 2L)
  expect_equal(
    dim(fit$covariate_diffusion_data$covariate_vertex_time),
    c(ncol(fit$tmb_data$A_st), fit$tmb_data$n_t, 2L)
  )
})

test_that("temporal covariate-diffusion vertex override can use extra_time without padded rows", {
  dat <- data.frame(
    y = rnorm(6),
    x1 = rnorm(6),
    year = c(1, 1, 2, 2, 4, 4),
    X = rep(1:3, each = 2),
    Y = rep(c(0, 1), 3)
  )
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
  vertex_time <- matrix(1, nrow = ncol(mesh$A_st), ncol = 4L)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    extra_time = 3,
    spatial = "off",
    spatiotemporal = "off",
    covariate_diffusion = ~ time(x1),
    experimental = list(covariate_diffusion_covariate_vertex = vertex_time),
    do_fit = FALSE
  )

  expect_equal(dim(fit$covariate_diffusion_data$covariate_vertex_time), c(ncol(mesh$A_st), 4L, 1L))
})
