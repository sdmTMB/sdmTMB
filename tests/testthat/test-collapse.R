test_that("collapsing spatial and spatiotemporal fields works", {
  skip_on_cran()
  skip_on_ci()

  set.seed(123)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:5, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0,
    phi = 0.2,
    sigma_O = 0,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  # create some fake 0s
  sim_dat$observed[sample(1:500, size = 50, replace = FALSE)] <- 0

  # test spatial collapse with gaussian family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off",
    control = sdmTMBcontrol(collapse_spatial_variance = FALSE)
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off",
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial / spatiotemporal collapse with spatiotemporal on
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off"
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial collapse with delta family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off", family = delta_gamma()
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", family = delta_gamma(),
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))

  # test spatial / spatiotemporal collapse with delta family
  fit_nospatial <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatiotemporal = "off", spatial = "off", family = delta_gamma()
  )
  fit <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    family = delta_gamma(),
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
  )
  expect_equal(tidy(fit_nospatial), tidy(fit))
  expect_equal(tidy(fit_nospatial, "ran_pars"), tidy(fit, "ran_pars"))
  expect_equal(tidy(fit_nospatial, "ran_vals"), tidy(fit, "ran_vals"))
})

test_that("custom collapse threshold works", {
  skip_on_cran()

  set.seed(456)

  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000), year = sample(1:5, size = 1000, replace = TRUE)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = 2,
    sigma_E = 0.2,
    phi = 0.2,
    sigma_O = 0,
    seed = 43,
    B = c(0.2, -0.4)
  )

  # With 0.001, should not collapse
  fit_default <- sdmTMB(observed ~ a1,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off",
    control = sdmTMBcontrol(
      collapse_spatial_variance = TRUE,
      collapse_spatial_variance_threshold = 0.001
    )
  )

  # With higher threshold (0.3), should collapse
  # Original call arguments can be objects in the caller's environment.
  fit_collapse <- local({
    form <- observed ~ a1
    family_obj <- gaussian()
    data_obj <- sim_dat
    mesh_obj <- mesh
    ctrl <- sdmTMBcontrol(
      collapse_spatial_variance = TRUE,
      collapse_spatial_variance_threshold = 0.3
    )
    sdmTMB(form,
      data = data_obj, mesh = mesh_obj, time = "year",
      family = family_obj, spatial = "off", control = ctrl
    )
  })

  expect_true(all(fit_default$spatiotemporal == "iid"))
  expect_true(all(fit_collapse$spatiotemporal == "off"))
})

test_that("collapse validation works", {
  ctrl <- sdmTMBcontrol()
  expect_equal(ctrl$collapse_spatial_variance_threshold, 0.01)
  expect_false("collapse_threshold" %in% names(ctrl))
  expect_false(ctrl$collapse_spatiotemporal_ar1)
  expect_equal(ctrl$collapse_ar1_threshold, 0.01)

  # Test that collapse_spatial_variance_threshold must be positive
  expect_error(
    sdmTMBcontrol(collapse_spatial_variance_threshold = -0.01),
    "collapse_spatial_variance_threshold not greater than 0"
  )

  expect_error(
    sdmTMBcontrol(collapse_spatial_variance_threshold = 0),
    "collapse_spatial_variance_threshold not greater than 0"
  )

  # Valid thresholds should work
  ctrl <- sdmTMBcontrol(collapse_spatial_variance_threshold = 0.001)
  expect_equal(ctrl$collapse_spatial_variance_threshold, 0.001)

  ctrl <- sdmTMBcontrol(collapse_spatial_variance_threshold = 0.1)
  expect_equal(ctrl$collapse_spatial_variance_threshold, 0.1)

  expect_warning(
    ctrl <- sdmTMBcontrol(collapse_threshold = 0.1),
    "deprecated"
  )
  expect_equal(ctrl$collapse_spatial_variance_threshold, 0.1)

  expect_warning(
    ctrl <- sdmTMBcontrol(
      collapse_threshold = 0.1,
      collapse_spatial_variance_threshold = 0.2
    ),
    "deprecated"
  )
  expect_equal(ctrl$collapse_spatial_variance_threshold, 0.2)

  ctrl <- sdmTMBcontrol(
    collapse_spatiotemporal_ar1 = TRUE,
    collapse_ar1_threshold = 0.02
  )
  expect_true(ctrl$collapse_spatiotemporal_ar1)
  expect_equal(ctrl$collapse_ar1_threshold, 0.02)

  for (x in list(NA, c(TRUE, FALSE), 1)) {
    expect_error(sdmTMBcontrol(collapse_spatiotemporal_ar1 = x))
  }
  for (x in list(NA_real_, Inf, numeric(0), c(0.01, 0.02), 0, -0.01, 0.5, 1)) {
    expect_error(sdmTMBcontrol(collapse_ar1_threshold = x))
  }
})

test_that("collapse checks ignore disabled delta fields", {
  args <- list(
    report_vals = list(
      sigma_O = matrix(c(0, 0.5), nrow = 1L),
      sigma_E = matrix(c(rep(0.001, 3L), rep(0.5, 3L)), nrow = 3L)
    ),
    spatial = c("off", "on"),
    spatiotemporal = c("off", "ar1"),
    n_m = 2L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = TRUE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = FALSE,
    collapse_ar1_threshold = 0.01,
    delta = TRUE,
    silent = TRUE
  )

  result <- do.call(check_and_collapse_spatial_fields, args)

  expect_false(result$do_refit)
  expect_equal(result$spatial_arg, list("off", "on"))
  expect_equal(result$spatiotemporal_arg, list("off", "ar1"))
})

test_that("collapse checks can disable one delta field", {
  args <- list(
    report_vals = list(
      sigma_O = matrix(c(0.001, 0.5), nrow = 1L),
      sigma_E = matrix(c(rep(0.001, 3L), rep(0.5, 3L)), nrow = 3L)
    ),
    spatial = c("on", "on"),
    spatiotemporal = c("ar1", "ar1"),
    n_m = 2L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = TRUE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = FALSE,
    collapse_ar1_threshold = 0.01,
    delta = TRUE,
    silent = TRUE
  )

  result <- do.call(check_and_collapse_spatial_fields, args)

  expect_true(result$do_refit)
  expect_equal(result$spatial_arg, list("off", "on"))
  expect_equal(result$spatiotemporal_arg, list("off", "ar1"))

  args$spatial <- unlist(result$spatial_arg)
  args$spatiotemporal <- unlist(result$spatiotemporal_arg)
  args$report_vals$sigma_O[1L] <- 0
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_false(result$do_refit)
})

test_that("spatiotemporal AR1 collapse respects rho boundaries", {
  args <- list(
    report_vals = list(
      sigma_O = 1,
      sigma_E = matrix(rep(1, 3L), nrow = 3L),
      rho = 0
    ),
    spatial = "off",
    spatiotemporal = "ar1",
    n_m = 1L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = FALSE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = TRUE,
    collapse_ar1_threshold = 0.01,
    delta = FALSE,
    silent = TRUE
  )

  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_true(result$do_refit)
  expect_identical(result$spatiotemporal_arg, "iid")

  args$report_vals$rho <- -0.01
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_identical(result$spatiotemporal_arg, "iid")

  args$report_vals$rho <- 0.99
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_identical(result$spatiotemporal_arg, "rw")

  for (rho in c(-0.0101, 0.0101, 0.9899, -0.999)) {
    args$report_vals$rho <- rho
    result <- do.call(check_and_collapse_spatial_fields, args)
    expect_false(result$do_refit)
    expect_identical(result$spatiotemporal_arg, "ar1")
  }

  args$collapse_spatiotemporal_ar1 <- FALSE
  args$report_vals$rho <- 0
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_false(result$do_refit)
  expect_identical(result$spatiotemporal_arg, "ar1")
})

test_that("spatiotemporal AR1 collapse handles delta components", {
  args <- list(
    report_vals = list(
      sigma_O = c(1, 1),
      sigma_E = matrix(rep(1, 6L), nrow = 3L),
      rho = c(0, 0.5)
    ),
    spatial = c("off", "off"),
    spatiotemporal = c("ar1", "ar1"),
    n_m = 2L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = FALSE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = TRUE,
    collapse_ar1_threshold = 0.01,
    delta = TRUE,
    silent = TRUE
  )

  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_equal(result$spatiotemporal_arg, list("iid", "ar1"))

  args$report_vals$rho <- c(0.999, -0.001)
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_equal(result$spatiotemporal_arg, list("rw", "iid"))

  args$spatiotemporal <- c("off", "ar1")
  args$report_vals$rho <- c(0, 0.5)
  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_false(result$do_refit)
  expect_equal(result$spatiotemporal_arg, list("off", "ar1"))
})

test_that("variance collapse takes precedence over AR1 collapse", {
  args <- list(
    report_vals = list(
      sigma_O = 1,
      sigma_E = matrix(rep(0.001, 3L), nrow = 3L),
      rho = 0
    ),
    spatial = "off",
    spatiotemporal = "ar1",
    n_m = 1L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = TRUE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = TRUE,
    collapse_ar1_threshold = 0.01,
    delta = FALSE,
    silent = TRUE
  )

  result <- do.call(check_and_collapse_spatial_fields, args)
  expect_true(result$do_refit)
  expect_identical(result$spatiotemporal_arg, "off")
})

test_that("spatiotemporal AR1 collapse ignores unavailable rho estimates", {
  args <- list(
    report_vals = list(
      sigma_O = c(1, 1),
      sigma_E = matrix(rep(1, 6L), nrow = 3L)
    ),
    spatial = c("off", "off"),
    spatiotemporal = c("ar1", "ar1"),
    n_m = 2L,
    n_t = 3L,
    omit_spatial_intercept = FALSE,
    collapse_spatial_variance = FALSE,
    collapse_spatial_variance_threshold = 0.01,
    collapse_spatiotemporal_ar1 = TRUE,
    collapse_ar1_threshold = 0.01,
    delta = TRUE,
    silent = TRUE
  )

  for (rho in list(NULL, numeric(0), NA_real_, c(NaN, Inf))) {
    args$report_vals$rho <- rho
    result <- do.call(check_and_collapse_spatial_fields, args)
    expect_false(result$do_refit)
    expect_equal(result$spatiotemporal_arg, list("ar1", "ar1"))
  }
})

test_that("spatiotemporal AR1 collapse refits and predicts", {
  skip_on_cran()
  skip_on_ci()

  fit_at_rho <- function(rho) {
    control <- sdmTMBcontrol(
      multiphase = FALSE,
      newton_loops = 0,
      getsd = FALSE,
      collapse_spatiotemporal_ar1 = TRUE,
      collapse_ar1_threshold = 0.01,
      start = list(ar1_phi = stats::qlogis((rho + 1) / 2)),
      map = list(ar1_phi = factor(NA))
    )
    sdmTMB(
      density ~ 1,
      data = pcod_2011,
      time = "year",
      mesh = pcod_mesh_2011,
      spatial = "off",
      spatiotemporal = "ar1",
      family = gaussian(),
      control = control,
      silent = TRUE
    )
  }

  fit_iid <- fit_at_rho(0)
  expect_identical(fit_iid$spatiotemporal, "iid")
  expect_equal(nrow(predict(fit_iid)), nrow(pcod_2011))

  fit_rw <- fit_at_rho(0.99)
  expect_identical(fit_rw$spatiotemporal, "rw")
  expect_equal(nrow(predict(fit_rw)), nrow(pcod_2011))
})
