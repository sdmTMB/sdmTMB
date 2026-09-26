test_that("tidy() works with basic spatiotemporal model", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)
})

test_that("tidy() works with spatial varying coefficients", {
  skip_on_cran()

  d <- pcod_2011
  d$year_scaled <- as.numeric(scale(d$year))
  fit <- sdmTMB(
    density ~ year_scaled,
    data = d, mesh = pcod_mesh_2011,
    family = tweedie(link = "log"), time = "year",
    spatiotemporal = "off", spatial_varying = ~ 0 + year_scaled
  )

  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)
})

test_that("tidy() works with delta models", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)

  x <- tidy(fit, "ran_pars", conf.int = TRUE, model = 2)
  expect_true(sum(is.na(x$std.error)) == 0L)
})

test_that("tidy() works with separate range parameters", {
  skip_on_cran()
  skip_on_ci()

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod_2011, mesh = pcod_mesh_2011, share_range = FALSE,
    family = tweedie(link = "log"), time = "year",
  )
  x <- tidy(fit, "ran_pars", conf.int = TRUE)
  expect_true(sum(is.na(x$std.error)) == 0L)
})

test_that("tidy() works with time-varying coefficients", {
  skip_on_cran()
  skip_on_ci()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ 0 + as.factor(year),
    time_varying = ~ 0 + depth_scaled + depth_scaled2,
    data = pcod_2011, time = "year", mesh = mesh,
    family = tweedie()
  )
  pars <- tidy(fit, "ran_vals")
  expect_equal(pars$estimate, c(
    -1.75, -1.62, -1.50, -2.21,
    -1.92, -0.92, -1.59, -2.20
  ), tolerance = 0.01)
  expect_equal(pars$term, c(
    "depth_scaled:2011", "depth_scaled:2013", "depth_scaled:2015", "depth_scaled:2017",
    "depth_scaled2:2011", "depth_scaled2:2013", "depth_scaled2:2015", "depth_scaled2:2017"
  ))
})

test_that("tidy() works with smooth terms", {
  skip_on_cran()
  skip_on_ci()

  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(link = "log")
  )
  x <- tidy(fit)
  expect_equal(x$term, c("(Intercept)", "sdepth"))
})

test_that("tidy() works with sdmTMB_cv objects", {
  skip_on_cran()
  skip_on_ci()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  m_cv <- sdmTMB_cv(
    density ~ 0 + depth_scaled + depth_scaled2,
    data = pcod, mesh = mesh,
    family = tweedie(link = "log"), k_folds = 2
  )
  model1 <- tidy(m_cv$models[[1]])
  allmodels <- tidy(m_cv)
  expect_equal(allmodels[1:2,1:5], model1[,1:5])
  expect_true("cv_split" %in% names(allmodels))
  expect_equal(allmodels$cv_split, sort(rep(1:2,2)))
})

test_that("tidy() correctly handles random effects standard deviations", {
  skip_on_cran()
  skip_on_ci()

  .pcod <- pcod
  .pcod$binned_lon <- round((pcod$lon - min(pcod$lon)) * 5)
  fit <- sdmTMB(data = .pcod,
                formula = density ~ (depth_scaled+1 + binned_lon|year) + (1 | binned_lon),
                family = tweedie(link = "log"), spatial = "off")
  t1 <- tidy(fit, "ran_pars")

  # Test parameter structure
  expect_identical(t1$term, c("phi", "tweedie_p", "sd__(Intercept)","sd__depth_scaled","sd__binned_lon","sd__(Intercept)"))
  expect_identical(t1$group_name[3:6], c("year", "year", "year", "binned_lon"))

  # Test that all SD estimates are positive
  sd_rows <- grepl("^sd__", t1$term)
  expect_true(all(t1$estimate[sd_rows] > 0), "All SD estimates should be positive")
  expect_true(all(t1$conf.low[sd_rows] > 0), "All SD confidence intervals should be positive")
  expect_true(all(t1$conf.high[sd_rows] > 0), "All SD confidence intervals should be positive")

  # Test that confidence intervals are properly ordered
  expect_true(all(t1$conf.low[sd_rows] < t1$estimate[sd_rows]), "CI lower bounds should be less than estimates")
  expect_true(all(t1$estimate[sd_rows] < t1$conf.high[sd_rows]), "Estimates should be less than CI upper bounds")

  # Test specific estimate values and CIs (regression tests)
  expect_equal(t1$estimate[t1$term == "sd__(Intercept)" & t1$group_name == "year"], 0.925, tolerance = 1e-3)
  expect_equal(t1$conf.low[t1$term == "sd__(Intercept)" & t1$group_name == "year"], 0.374, tolerance = 1e-2)
  expect_equal(t1$conf.high[t1$term == "sd__(Intercept)" & t1$group_name == "year"], 2.29, tolerance = 1e-2)

  expect_equal(t1$estimate[t1$term == "sd__depth_scaled"], 0.565, tolerance = 1e-3)
  expect_equal(t1$conf.low[t1$term == "sd__depth_scaled"], 0.333, tolerance = 1e-2)
  expect_equal(t1$conf.high[t1$term == "sd__depth_scaled"], 0.959, tolerance = 1e-2)
})

test_that("tidy() correctly handles smoother standard deviations", {
  skip_on_cran()
  skip_on_ci()

  .pcod <- pcod
  .pcod$fyear <- as.factor(pcod$year)
  fit <- sdmTMB(data = .pcod,
                formula = density ~ s(depth_scaled, by = fyear),
                family = tweedie(link = "log"), spatial = "off")
  t2 <- tidy(fit, "ran_par")

  # Test parameter structure
  expect_identical(t2$term, c("phi", "tweedie_p", "sd__s(depth_scaled):fyear2003", "sd__s(depth_scaled):fyear2004",
    "sd__s(depth_scaled):fyear2005", "sd__s(depth_scaled):fyear2007",
    "sd__s(depth_scaled):fyear2009", "sd__s(depth_scaled):fyear2011",
    "sd__s(depth_scaled):fyear2013", "sd__s(depth_scaled):fyear2015",
    "sd__s(depth_scaled):fyear2017"))

  # Test that all smooth SD estimates are positive
  smooth_sd_rows <- grepl("^sd__s\\(", t2$term)
  expect_true(all(t2$estimate[smooth_sd_rows] > 0), "All smooth SD estimates should be positive")
  expect_true(all(t2$conf.low[smooth_sd_rows] > 0), "All smooth SD confidence intervals should be positive")
  expect_true(all(t2$conf.high[smooth_sd_rows] > 0), "All smooth SD confidence intervals should be positive")

  # Test that confidence intervals are properly ordered for smooth SDs
  expect_true(all(t2$conf.low[smooth_sd_rows] < t2$estimate[smooth_sd_rows]), "CI lower bounds should be less than estimates")
  expect_true(all(t2$estimate[smooth_sd_rows] < t2$conf.high[smooth_sd_rows]), "Estimates should be less than CI upper bounds")

  # Test that std.error is NA for smooth SDs (since they're transformed from log-space)
  expect_true(all(is.na(t2$std.error[smooth_sd_rows])), "Smooth SD std.error should be NA")

  # Test specific smooth SD estimate values and CIs (regression tests)
  expect_equal(t2$estimate[t2$term == "sd__s(depth_scaled):fyear2003"], 10.01841, tolerance = 1e-3)
  expect_equal(t2$conf.low[t2$term == "sd__s(depth_scaled):fyear2003"], 3.414768, tolerance = 1e-2)
  expect_equal(t2$conf.high[t2$term == "sd__s(depth_scaled):fyear2003"], 29.39249, tolerance = 1e-2)

  expect_equal(t2$estimate[t2$term == "sd__s(depth_scaled):fyear2011"], 13.56432, tolerance = 1e-3)
  expect_equal(t2$conf.low[t2$term == "sd__s(depth_scaled):fyear2011"], 6.740783, tolerance = 1e-2)
  expect_equal(t2$conf.high[t2$term == "sd__s(depth_scaled):fyear2011"], 27.29518, tolerance = 1e-2)
})

test_that("tidy() works with delta model with random intercepts and AR1 time series", {
  skip_on_cran()

  d <- pcod
  d$fake <- rep(c("a", "b", "c"), 9999)[1:nrow(d)]
  # This fixture is a marginal fit (non-positive-definite Hessian); it's only
  # used to check tidy() output structure, not convergence quality.
  fit <- suppressWarnings(sdmTMB(
    density ~ breakpt(depth_scaled) + (1|fake),
    data = d,
    time = "year",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    time_varying = ~1,
    time_varying_type = "ar1",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(type = "poisson-link")
  ))
  b <- tidy(fit, effects = "ran_pars")
  b <- tidy(fit, effects = "ran_vals")
  expect_identical(b$term, c("(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)",
    "(Intercept)", "(Intercept)", "(Intercept):2003", "(Intercept):2004",
    "(Intercept):2005", "(Intercept):2006", "(Intercept):2007", "(Intercept):2008",
    "(Intercept):2009", "(Intercept):2010", "(Intercept):2011", "(Intercept):2012",
    "(Intercept):2013", "(Intercept):2014", "(Intercept):2015", "(Intercept):2016",
    "(Intercept):2017", "(Intercept):2003", "(Intercept):2004", "(Intercept):2005",
    "(Intercept):2006", "(Intercept):2007", "(Intercept):2008", "(Intercept):2009",
    "(Intercept):2010", "(Intercept):2011", "(Intercept):2012", "(Intercept):2013",
    "(Intercept):2014", "(Intercept):2015", "(Intercept):2016", "(Intercept):2017"
  ))
})

test_that("tidy() correctly handles anisotropic ranges", {
  skip_on_cran()
  skip_on_ci()

  # Test 1: Spatial-only anisotropic model (should have range_min and range_max)
  fit_sp_only <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "on",
    spatiotemporal = "off",
    anisotropy = TRUE
  )

  t1 <- tidy(fit_sp_only, "ran_pars", conf.int = TRUE)

  # Should not have "range" term
  expect_false("range" %in% t1$term)

  # Should have range_min and range_max
  expect_true("range_min" %in% t1$term)
  expect_true("range_max" %in% t1$term)

  # std.error and CIs should be NA for anisotropic ranges
  expect_true(is.na(t1$std.error[t1$term == "range_min"]))
  expect_true(is.na(t1$std.error[t1$term == "range_max"]))
  expect_true(is.na(t1$conf.low[t1$term == "range_min"]))
  expect_true(is.na(t1$conf.high[t1$term == "range_min"]))

  # Check that values match plot_anisotropy
  aniso_df <- plot_anisotropy(fit_sp_only, return_data = TRUE)
  aniso_dat_sp <- aniso_df[aniso_df$random_field == "spatial", ][1, ]
  expect_equal(t1$estimate[t1$term == "range_min"], aniso_dat_sp$b)
  expect_equal(t1$estimate[t1$term == "range_max"], aniso_dat_sp$a)

  # Test 2: Model with separate spatial and spatiotemporal ranges
  test_mesh <- make_mesh(data = dogfish, xy_cols = c("X", "Y"), cutoff = 10)
  fit_separate <- sdmTMB(
    data = dogfish,
    formula = catch_weight ~ 1,
    offset = log(dogfish$area_swept),
    mesh = test_mesh,
    family = tweedie(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 0)
  )

  t2 <- tidy(fit_separate, "ran_pars", conf.int = TRUE)

  # Should not have "range" term
  expect_false("range" %in% t2$term)

  # Should have separate spatial and spatiotemporal ranges
  expect_true("range_spatial_min" %in% t2$term)
  expect_true("range_spatial_max" %in% t2$term)
  expect_true("range_spatiotemporal_min" %in% t2$term)
  expect_true("range_spatiotemporal_max" %in% t2$term)

  # Check that values match plot_anisotropy
  aniso_df2 <- plot_anisotropy(fit_separate, return_data = TRUE)
  aniso_dat2_sp <- aniso_df2[aniso_df2$random_field == "spatial", ][1, ]
  aniso_dat2_st <- aniso_df2[aniso_df2$random_field == "spatiotemporal", ][1, ]
  expect_equal(t2$estimate[t2$term == "range_spatial_min"], aniso_dat2_sp$b)
  expect_equal(t2$estimate[t2$term == "range_spatial_max"], aniso_dat2_sp$a)
  expect_equal(t2$estimate[t2$term == "range_spatiotemporal_min"], aniso_dat2_st$b)
  expect_equal(t2$estimate[t2$term == "range_spatiotemporal_max"], aniso_dat2_st$a)

  # Test 3: Delta model with shared range and anisotropy
  fit_delta_shared <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = delta_gamma(),
    share_range = TRUE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1L)
  )

  # Test model 1
  t3_m1 <- tidy(fit_delta_shared, "ran_pars", conf.int = TRUE, model = 1)
  expect_false("range" %in% t3_m1$term)
  expect_true("range_min" %in% t3_m1$term)
  expect_true("range_max" %in% t3_m1$term)

  aniso_df3_m1 <- plot_anisotropy(fit_delta_shared, return_data = TRUE)
  aniso_dat3_m1_sp <- aniso_df3_m1[aniso_df3_m1$model_num == 1 & aniso_df3_m1$random_field == "spatial", ][1, ]
  expect_equal(t3_m1$estimate[t3_m1$term == "range_min"], aniso_dat3_m1_sp$b)
  expect_equal(t3_m1$estimate[t3_m1$term == "range_max"], aniso_dat3_m1_sp$a)

  # Test model 2
  t3_m2 <- tidy(fit_delta_shared, "ran_pars", conf.int = TRUE, model = 2)
  expect_false("range" %in% t3_m2$term)
  expect_true("range_min" %in% t3_m2$term)
  expect_true("range_max" %in% t3_m2$term)

  aniso_dat3_m2_sp <- aniso_df3_m1[aniso_df3_m1$model_num == 2 & aniso_df3_m1$random_field == "spatial", ][1, ]
  expect_equal(t3_m2$estimate[t3_m2$term == "range_min"], aniso_dat3_m2_sp$b)
  expect_equal(t3_m2$estimate[t3_m2$term == "range_max"], aniso_dat3_m2_sp$a)

  # Test 4: Non-anisotropic model should still have "range" term
  fit_iso <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    anisotropy = FALSE
  )

  t4 <- tidy(fit_iso, "ran_pars", conf.int = TRUE)
  expect_true("range" %in% t4$term)
  expect_false("range_min" %in% t4$term)
  expect_false("range_max" %in% t4$term)
})
