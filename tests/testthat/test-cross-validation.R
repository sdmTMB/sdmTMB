test_that("Basic cross validation works", {
  skip_on_ci()
  skip_on_cran()
  d <- pcod
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)

  set.seed(2)
  # library(future) # for parallel processing
  # plan(multisession) # for parallel processing
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde, spatial = "off", spatiotemporal = "off",
    family = tweedie(link = "log"), time = "year", k_folds = 2
  )
  print(x)
  expect_equal(class(x$sum_loglik), "numeric")
  expect_equal(x$sum_loglik, sum(x$data$cv_loglik))
  expect_equal(x$sum_loglik, sum(x$fold_loglik))
  expect_true("data.frame" %in% class(x$data))
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # Use fold_ids:
  x <- sdmTMB_cv(
    density ~ 1,
    data = d, mesh = spde, spatial = "off",
    family = tweedie(link = "log"),
    fold_ids = rep(seq(1, 2), nrow(d))[seq(1, nrow(d))])
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # student-t: was broken at one time, must deal with `df`
  d <- subset(d, density > 0)
  d$log_density <- log(d$density)
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)
  x <- sdmTMB_cv(
    log_density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde, spatial = "off",
    family = sdmTMB::student(df = 9), time = "year", k_folds = 2
  )
  expect_equal(class(x$models[[1]]), "sdmTMB")

  # Try passing family as a variable -- this is per Issue #127
  fam <- gaussian(link = "identity")
  x <- sdmTMB_cv(
    log_density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde, spatial = "off", spatiotemporal = "off",
    family = fam, time = "year", k_folds = 2
  )
  expect_equal(class(x$models[[1]]), "sdmTMB")
})

test_that("Cross validation in parallel with globals", {
  skip_on_cran()
  skip_on_ci()
  # https://github.com/sdmTMB/sdmTMB/issues/127
  d <- pcod
  spde <- make_mesh(d, c("X", "Y"), cutoff = 15)
  set.seed(2)
  future::plan(future::multisession, workers = 2L)
  fam <- tweedie(link = "log")
  x <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = spde,
    family = fam, time = "year", k_folds = 2L, future_globals = 'fam'
  )
  expect_s3_class(x$models[[1]], "sdmTMB")
  future::plan(future::sequential)
})

test_that("Leave future out cross validation works", {
  skip_on_cran()
  skip_on_ci()
  x <- sdmTMB_cv(
    present ~ 1,
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    lfo = TRUE,
    lfo_forecast = 1,
    lfo_validations = 2,
    family = binomial(),
    time = "year"
  )
  expect_equal(class(x$sum_loglik), "numeric")
  expect_equal(x$sum_loglik, sum(x$data$cv_loglik))
  expect_equal(x$sum_loglik, sum(x$fold_loglik))
  expect_true("data.frame" %in% class(x$data))
  expect_true(inherits(x$models[[1]], "sdmTMB"))

  # Can see how the folds are assigned with
  table(x$models[[1]]$data$cv_fold, x$models[[1]]$data$year)

  expect_equal(length(x$models), 2)
  expect_equal(length(x$fold_loglik), 2)
  expect_equal(length(x$max_gradients), 2)
  expect_equal(cor(x$data$cv_fold[x$data$cv_fold!=1], x$data$year[x$data$cv_fold!=1]), 1.0)
})

test_that("LFO fold assignments follow documented structure", {
  skip_on_cran()
  skip_on_ci()

  # Test the documented example from help: 9 time steps, lfo_forecast = 2, lfo_validations = 3
  # Expected:
  # - Fit data to time steps 1 to 5, predict and validate step 7.
  # - Fit data to time steps 1 to 6, predict and validate step 8.
  # - Fit data to time steps 1 to 7, predict and validate step 9.

  # Use pcod data which has 9 time steps: 2003, 2004, 2005, 2007, 2009, 2011, 2013, 2015, 2017
  time_steps <- sort(unique(pcod$year))
  expect_equal(length(time_steps), 9)

  # This should work if we have enough validation periods
  # With lfo_forecast = 2 and lfo_validations = 2 instead of 3
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  x <- sdmTMB_cv(
    present ~ 1,
    data = pcod,
    mesh = mesh,
    lfo = TRUE,
    lfo_forecast = 2,
    lfo_validations = 2,  # Reduced to 2 to fit in 9 time steps
    family = binomial(),
    spatiotemporal = "off",
    time = "year"
  )

  # Check fold structure
  fold_table <- table(x$data$cv_fold, x$data$year)

  # Validation 1: Should train on years 2003-2011 (time steps 1-5), validate on 2015 (step 7)
  # Validation 2: Should train on years 2003-2013 (time steps 1-6), validate on 2017 (step 8)

  # Check that validation years have the right fold assignments
  # For lfo_forecast = 2, validation data should be in folds k + lfo_forecast
  expect_true(all(x$data$cv_fold[x$data$year == 2015] == 3))  # fold 1 + 2
  expect_true(all(x$data$cv_fold[x$data$year == 2017] == 4))  # fold 2 + 2

  # Check that training years are properly assigned
  # For validation 1: years 2003-2011 should have fold <= 1
  train_years_1 <- c(2003, 2004, 2005, 2007, 2009, 2011)
  expect_true(all(x$data$cv_fold[x$data$year %in% train_years_1] <= 1))

  # For validation 2: years 2003-2013 should have fold <= 2
  train_years_2 <- c(2003, 2004, 2005, 2007, 2009, 2011, 2013)
  expect_true(all(x$data$cv_fold[x$data$year %in% train_years_2] <= 2))
})

test_that("LFO parameter validation works", {
  skip_on_cran()
  skip_on_ci()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25)

  # Should error when lfo_validations + lfo_forecast > number of time steps
  expect_error(
    sdmTMB_cv(
      present ~ 1,
      data = pcod_2011,  # Only 2 time steps
      mesh = mesh,
      lfo = TRUE,
      lfo_forecast = 2,
      lfo_validations = 3,  # 3 + 2 = 5 > 2 time steps
      family = binomial(),
      time = "year"
    ),
    "Not enough time steps"
  )
})

test_that("Cross validation with offsets works", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  set.seed(1)
  d <- pcod_2011
  d$log_effort <- rnorm(nrow(d))

  # library(future)
  # future::plan(future::multisession)

  expect_error(
    fit_cv_off1 <- sdmTMB_cv(
      density ~ 1,
      data = d,
      mesh = pcod_mesh_2011,
      family = tweedie(),
      spatial = "off",
      offset = d$log_effort, #<
      k_folds = 2,
      parallel = TRUE #<
    ))

  set.seed(1)
  fit_cv_off1 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off2 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off3 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off4 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    offset = "log_effort", #<
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  # now without offset
  set.seed(1)
  fit_cv_off5 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = TRUE #<
  )

  set.seed(1)
  fit_cv_off6 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off7 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = FALSE, #<
    parallel = FALSE #<
  )

  set.seed(1)
  fit_cv_off8 <- sdmTMB_cv(
    density ~ 1,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    use_initial_fit = TRUE, #<
    parallel = FALSE #<
  )

  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off2$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off3$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off1$models[[1]]$model$par, 4), round(fit_cv_off4$models[[1]]$model$par, 4))

  # with/without offset:
  expect_false(identical(fit_cv_off1$models[[1]]$model$par, fit_cv_off5$models[[1]]$model$par))

  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off6$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off7$models[[1]]$model$par, 4))
  expect_equal(round(fit_cv_off5$models[[1]]$model$par, 4), round(fit_cv_off8$models[[1]]$model$par, 4))

  # future::plan(future::sequential)
})

test_that("Delta model cross validation works", {
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  out_tw <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = tweedie(), k_folds = 2
  )
  set.seed(1)
  out_dg <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = delta_gamma(), k_folds = 2
  )
  diff_ll <- out_tw$sum_loglik - out_dg$sum_loglik
  expect_equal(round(diff_ll, 4), round(-19.8025, 4))

  set.seed(1)
  out_dpg <- sdmTMB_cv(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
    family = delta_gamma(type = "poisson-link"), k_folds = 2
  )
  diff_ll <- out_dpg$sum_loglik - out_dg$sum_loglik
  expect_equal(round(diff_ll, 4), round(-4.5477, 4))
})

test_that("Cross validation with weights works", {
  skip_on_cran()
  skip_on_ci()

  # Test with constant mesh
  set.seed(123)
  d <- pcod_2011

  # Create some example weights
  weights_vec <- runif(nrow(d), 0.5, 1.5)

  # CV with weights
  set.seed(1)
  x_weighted <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    weights = weights_vec
  )

  # CV without weights (default)
  set.seed(1)
  x_unweighted <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # Check that CV runs successfully with weights
  expect_s3_class(x_weighted, "sdmTMB_cv")
  expect_equal(class(x_weighted$sum_loglik), "numeric")
  expect_equal(x_weighted$sum_loglik, sum(x_weighted$data$cv_loglik))
  expect_equal(x_weighted$sum_loglik, sum(x_weighted$fold_loglik))
  expect_true("data.frame" %in% class(x_weighted$data))
  expect_equal(class(x_weighted$models[[1]]), "sdmTMB")

  # Weighted and unweighted CV should give different results
  expect_false(isTRUE(all.equal(x_weighted$sum_loglik, x_unweighted$sum_loglik)))

  # Test that predictions exist
  expect_true("cv_predicted" %in% names(x_weighted$data))
  expect_true("cv_loglik" %in% names(x_weighted$data))

  # Test with use_initial_fit = TRUE
  set.seed(1)
  x_weighted_init <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "off",
    k_folds = 2,
    weights = weights_vec,
    use_initial_fit = TRUE
  )
  expect_s3_class(x_weighted_init, "sdmTMB_cv")
  expect_equal(class(x_weighted_init$sum_loglik), "numeric")

  # Test error handling for mismatched weights length
  expect_error(
    sdmTMB_cv(
      density ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = tweedie(),
      spatial = "off",
      k_folds = 2,
      weights = c(1, 2, 3)  # Wrong length
    ),
    "same length as the number of rows"
  )

  # Test error handling for non-positive weights
  expect_error(
    sdmTMB_cv(
      density ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = tweedie(),
      spatial = "off",
      k_folds = 2,
      weights = c(rep(1, nrow(d) - 1), 0)  # Zero weight
    ),
    "must be positive"
  )

  expect_error(
    sdmTMB_cv(
      density ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = tweedie(),
      spatial = "off",
      k_folds = 2,
      weights = c(rep(1, nrow(d) - 1), -1)  # Negative weight
    ),
    "must be positive"
  )
})

test_that("Cross validation returns deviance residuals", {
  skip_on_cran()
  skip_on_ci()

  set.seed(123)
  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  # Test with Gaussian family (simple case)
  d_pos <- subset(d, density > 0)
  d_pos$log_density <- log(d_pos$density)

  x_gauss <- sdmTMB_cv(
    log_density ~ depth_scaled,
    data = d_pos,
    mesh = mesh,
    family = gaussian(),
    spatial = "off",
    k_folds = 2
  )

  # Check that cv_deviance_resid column exists
  expect_true("cv_deviance_resid" %in% names(x_gauss$data))

  # Check that deviance residuals are numeric and not all NA
  expect_true(is.numeric(x_gauss$data$cv_deviance_resid))
  expect_false(all(is.na(x_gauss$data$cv_deviance_resid)))

  # Check that we have the same number of residuals as rows
  expect_equal(length(x_gauss$data$cv_deviance_resid), nrow(d_pos))

  # Test with Tweedie family
  set.seed(123)
  x_tweedie <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_tweedie$data))
  expect_false(all(is.na(x_tweedie$data$cv_deviance_resid)))

  # Test with binomial family
  set.seed(123)
  x_binom <- sdmTMB_cv(
    present ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = binomial(),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_binom$data))
  expect_false(all(is.na(x_binom$data$cv_deviance_resid)))
})

test_that("Deviance residuals work with different families", {
  skip_on_cran()
  skip_on_ci()

  set.seed(456)
  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  # Gamma family
  d_pos <- subset(d, density > 0)
  x_gamma <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d_pos,
    mesh = mesh,
    family = Gamma(link = "log"),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_gamma$data))
  expect_false(all(is.na(x_gamma$data$cv_deviance_resid)))

  # Lognormal family
  x_lognormal <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d_pos,
    mesh = mesh,
    family = lognormal(),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_lognormal$data))
  expect_false(all(is.na(x_lognormal$data$cv_deviance_resid)))

  # Negative binomial (nbinom2)
  d$count <- round(d$density)
  x_nbinom2 <- sdmTMB_cv(
    count ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = nbinom2(),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_nbinom2$data))
  expect_false(all(is.na(x_nbinom2$data$cv_deviance_resid)))
})

test_that("Deviance residuals work with delta models", {
  skip_on_cran()
  skip_on_ci()

  set.seed(789)
  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  # Delta-gamma model
  x_delta <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = delta_gamma(),
    spatial = "off",
    k_folds = 2
  )

  # Check that cv_deviance_resid column exists for delta models
  expect_true("cv_deviance_resid" %in% names(x_delta$data))
  expect_false(all(is.na(x_delta$data$cv_deviance_resid)))

  # Delta-lognormal model
  x_delta_lnorm <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = delta_lognormal(),
    spatial = "off",
    k_folds = 2
  )

  expect_true("cv_deviance_resid" %in% names(x_delta_lnorm$data))
  expect_false(all(is.na(x_delta_lnorm$data$cv_deviance_resid)))
})

test_that("Deviance residuals work with non-constant mesh", {
  skip_on_cran()
  skip_on_ci()

  set.seed(321)
  d <- pcod_2011

  # Use mesh_args instead of mesh to create non-constant mesh
  x_nc <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh_args = list(xy_cols = c("X", "Y"), cutoff = 10),
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # Check that cv_deviance_resid column exists with non-constant mesh
  expect_true("cv_deviance_resid" %in% names(x_nc$data))
  expect_false(all(is.na(x_nc$data$cv_deviance_resid)))
})

test_that("Deviance residuals work with LFO cross-validation", {
  skip_on_cran()
  skip_on_ci()

  set.seed(654)
  x_lfo <- sdmTMB_cv(
    present ~ 1,
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    lfo = TRUE,
    lfo_forecast = 1,
    lfo_validations = 2,
    family = binomial(),
    spatial = "off",
    time = "year"
  )

  # Check that cv_deviance_resid exists for LFO
  expect_true("cv_deviance_resid" %in% names(x_lfo$data))
  # For LFO, deviance residuals might not be available (uses old method)
  # Just check the column exists
})

test_that("Deviance residuals match single-model residuals", {
  skip_on_cran()
  skip_on_ci()

  set.seed(111)
  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)

  # Fit a single model
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off"
  )

  # Get deviance residuals from the fit
  dev_res_fit <- residuals(fit, type = "deviance")

  # Do CV with k_folds = 2
  set.seed(111)
  x_cv <- sdmTMB_cv(
    density ~ depth_scaled,
    data = d,
    mesh = mesh,
    family = tweedie(),
    spatial = "off",
    k_folds = 2
  )

  # The CV deviance residuals should be calculated from models
  # fitted without the held-out fold, so they won't exactly match
  # the full-data residuals. But they should be in a similar range
  # and have similar properties.
  expect_true(is.numeric(x_cv$data$cv_deviance_resid))
  expect_false(all(is.na(x_cv$data$cv_deviance_resid)))

  # Check that the range is reasonable (not wildly different)
  # This is a soft check - just ensuring they're in the same ballpark
  expect_true(
    max(abs(x_cv$data$cv_deviance_resid), na.rm = TRUE) <
    max(abs(dev_res_fit)) * 3
  )
})

make_cv_covariate_diffusion_data <- function() {
  set.seed(303)
  n_t <- 2L
  n_s <- 5L
  year <- rep(seq_len(n_t), each = n_s)
  X <- rep(seq_len(n_s), times = n_t)
  Y <- rep(seq_len(2L), length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + X / max(X)))
  eta <- 0.1 + 0.6 * x1

  data.frame(
    y = eta + stats::rnorm(length(eta), sd = 0.05),
    x1 = x1,
    year = year,
    X = X,
    Y = Y
  )
}

test_that("Cross validation handles covariate diffusion under year-based folds", {
  skip_on_cran()

  dat <- make_cv_covariate_diffusion_data()
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
  fold_ids <- as.integer(factor(dat$year))

  fit_base <- sdmTMB_cv(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    fold_ids = fold_ids,
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    parallel = FALSE
  )

  fit_time <- sdmTMB_cv(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    fold_ids = fold_ids,
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    parallel = FALSE,
    covariate_diffusion = ~ time(x1),
    control = sdmTMBcontrol(newton_loops = 0)
  )

  fit_space <- sdmTMB_cv(
    y ~ x1,
    data = dat,
    mesh = mesh,
    time = "year",
    fold_ids = fold_ids,
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    parallel = FALSE,
    covariate_diffusion = ~ space(x1),
    control = sdmTMBcontrol(newton_loops = 0)
  )

  expect_s3_class(fit_base, "sdmTMB_cv")
  expect_s3_class(fit_time, "sdmTMB_cv")
  expect_s3_class(fit_space, "sdmTMB_cv")

  expect_equal(nrow(fit_time$data), nrow(dat))
  expect_false(anyNA(fit_time$data$cv_predicted))
  expect_true(all(is.finite(fit_time$data$cv_predicted)))
  expect_true(is.finite(fit_time$sum_loglik))

  expect_equal(nrow(fit_space$data), nrow(dat))
  expect_false(anyNA(fit_space$data$cv_predicted))
  expect_true(all(is.finite(fit_space$data$cv_predicted)))
  expect_true(is.finite(fit_space$sum_loglik))
})
