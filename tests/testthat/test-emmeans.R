test_that("emmeans works with parametric models", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Basic parametric model without spatial component  
  fit <- sdmTMB(
    present ~ as.factor(year),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans creation
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test pairwise comparisons
  pairs_result <- pairs(emm)  # pairs is a method, not namespaced
  expect_s4_class(pairs_result, "emmGrid")
  
  # Test response scale
  emm_resp <- emmeans::emmeans(fit, ~ year, type = "response")
  expect_s4_class(emm_resp, "emmGrid")
  
  # Check that results have correct dimensions
  emm_summary <- summary(emm)
  expect_equal(nrow(emm_summary), length(unique(pcod_2011$year)))
  expect_true(all(c("year", "emmean", "SE", "df", "lower.CL", "upper.CL") %in% names(emm_summary)))
})

test_that("emmeans works with continuous covariates", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  fit <- sdmTMB(
    present ~ as.factor(year) + depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans with continuous covariate
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test emtrends for slopes
  emtrends_result <- emmeans::emtrends(fit, ~ year, var = "depth_scaled")
  expect_s4_class(emtrends_result, "emmGrid")
  
  # Test interaction
  fit_int <- sdmTMB(
    present ~ as.factor(year) * depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  emm_int <- emmeans::emmeans(fit_int, ~ year)
  expect_s4_class(emm_int, "emmGrid")
})

test_that("emmeans works with smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Model with smoother - this was previously broken
  fit_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test that emmeans works with smoothers
  emm_smooth <- emmeans::emmeans(fit_smooth, ~ year)
  expect_s4_class(emm_smooth, "emmGrid")
  
  # Test pairwise comparisons with smoothers
  pairs_smooth <- pairs(emm_smooth)
  expect_s4_class(pairs_smooth, "emmGrid")
  
  # Test response scale with smoothers  
  emm_smooth_resp <- emmeans::emmeans(fit_smooth, ~ year, type = "response")
  expect_s4_class(emm_smooth_resp, "emmGrid")
  
  # Check that results are reasonable
  emm_summary <- summary(emm_smooth)
  expect_equal(nrow(emm_summary), length(unique(pcod_2011$year)))
  expect_true(all(is.finite(emm_summary$emmean)))
  expect_true(all(emm_summary$SE > 0))
})

test_that("emmeans works with multiple smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Model with multiple smoothers
  fit_multi_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled) + s(X),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test that emmeans still works with multiple smoothers
  emm_multi <- emmeans::emmeans(fit_multi_smooth, ~ year)
  expect_s4_class(emm_multi, "emmGrid")
  
  # Test pairwise comparisons
  pairs_multi <- pairs(emm_multi)
  expect_s4_class(pairs_multi, "emmGrid")
  
  # Results should be sensible
  emm_summary <- summary(emm_multi)
  expect_true(all(is.finite(emm_summary$emmean)))
})

test_that("emmeans results are consistent between parametric and smooth models", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Parametric model  
  fit_param <- sdmTMB(
    present ~ as.factor(year) + depth_scaled,
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Smooth model
  fit_smooth <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off", 
    family = binomial()
  )
  
  # Get emmeans for both
  emm_param <- emmeans::emmeans(fit_param, ~ year)
  emm_smooth <- emmeans::emmeans(fit_smooth, ~ year)
  
  # Results should be reasonably similar (smooth vs linear depth effect)
  param_est <- summary(emm_param)$emmean
  smooth_est <- summary(emm_smooth)$emmean
  
  # They shouldn't be identical but should be correlated
  expect_true(cor(param_est, smooth_est) > 0.8)
})

test_that("emmeans works with delta models", {
  skip_if_not_installed("emmeans")
  library(emmeans)

  # Delta model
  fit_delta <- sdmTMB(
    density ~ as.factor(year),
    data = pcod_2011,
    spatial = "off",
    family = delta_gamma()
  )

  # Test binomial component (model = 1)
  emm1 <- emmeans::emmeans(fit_delta, ~ year, model = 1)
  expect_s4_class(emm1, "emmGrid")

  # Test positive component (model = 2)
  emm2 <- emmeans::emmeans(fit_delta, ~ year, model = 2)
  expect_s4_class(emm2, "emmGrid")

  # Check that results have correct dimensions
  emm1_summary <- summary(emm1)
  emm2_summary <- summary(emm2)
  expect_equal(nrow(emm1_summary), length(unique(pcod_2011$year)))
  expect_equal(nrow(emm2_summary), length(unique(pcod_2011$year)))

  # Test pairwise comparisons work
  pairs1 <- pairs(emm1)
  pairs2 <- pairs(emm2)
  expect_s4_class(pairs1, "emmGrid")
  expect_s4_class(pairs2, "emmGrid")

  # Test response scale works
  emm1_resp <- emmeans::emmeans(fit_delta, ~ year, model = 1, type = "response")
  emm2_resp <- emmeans::emmeans(fit_delta, ~ year, model = 2, type = "response")
  expect_s4_class(emm1_resp, "emmGrid")
  expect_s4_class(emm2_resp, "emmGrid")
})

test_that("emmeans works with binomial family and smoothers", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  
  # Test binomial with smoother - this is the main use case that was broken
  fit <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011,
    spatial = "off",
    family = binomial()
  )
  
  # Test emmeans works
  emm <- emmeans::emmeans(fit, ~ year)
  expect_s4_class(emm, "emmGrid")
  
  # Test that we can extract meaningful results
  emm_summary <- summary(emm)
  expect_true(nrow(emm_summary) == 4)  # 4 years
  expect_true(all(c("year", "emmean") %in% names(emm_summary)))
})
test_that("emmeans with smoothers matches predict()", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  yrs <- sort(unique(pcod_2011$year))
  check <- function(fit, specs, at, nd, ...) {
    e <- as.data.frame(emmeans::emmeans(fit, specs, at = at, ...))
    p <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA, ...)
    expect_equal(e$emmean, p$est, tolerance = 1e-6)
    expect_equal(e$SE, p$est_se, tolerance = 1e-6)
  }
  fit <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled),
    data = pcod_2011, spatial = "off", family = binomial()
  )
  check(fit, ~ year, list(depth_scaled = 0.5),
    data.frame(year = yrs, depth_scaled = 0.5))
  # the smoothed covariate is in the reference grid:
  check(fit, ~ depth_scaled | year, list(depth_scaled = c(-1, 1), year = 2013),
    data.frame(year = 2013, depth_scaled = c(-1, 1)))

  fit_by <- sdmTMB(
    present ~ as.factor(year) + s(depth_scaled, by = as.factor(year), k = 5),
    data = pcod_2011, spatial = "off", family = binomial()
  )
  check(fit_by, ~ year, list(depth_scaled = 0.5),
    data.frame(year = yrs, depth_scaled = 0.5))

  fit_delta <- sdmTMB(
    density ~ as.factor(year) + s(depth_scaled, k = 5),
    data = pcod_2011, spatial = "off", family = delta_gamma()
  )
  for (m in 1:2) {
    check(fit_delta, ~ year, list(depth_scaled = 0.5),
      data.frame(year = yrs, depth_scaled = 0.5), model = m)
  }
})

test_that("emmeans handles poly(), offsets, and delta smoothers in one component", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  d <- pcod_2011
  d$fyear <- as.factor(d$year)
  # poly() keeps its fitted coefficients and spans the same space as I():
  fit_poly <- sdmTMB(present ~ fyear + poly(depth_scaled, 2),
    data = d, spatial = "off", family = binomial())
  fit_I <- sdmTMB(present ~ fyear + depth_scaled + I(depth_scaled^2),
    data = d, spatial = "off", family = binomial())
  e_poly <- as.data.frame(emmeans::emmeans(fit_poly, ~ fyear, at = list(depth_scaled = 0.5)))
  e_I <- as.data.frame(emmeans::emmeans(fit_I, ~ fyear, at = list(depth_scaled = 0.5)))
  expect_equal(e_poly$emmean, e_I$emmean, tolerance = 1e-5)
  expect_equal(e_poly$SE, e_I$SE, tolerance = 1e-4)

  # EMMs are for an offset of 0 unless emmeans' `offset` is given:
  d$off <- log(seq(0.5, 2, length.out = nrow(d)))
  fit_off <- sdmTMB(present ~ fyear, offset = "off",
    data = d, spatial = "off", family = poisson())
  e0 <- as.data.frame(emmeans::emmeans(fit_off, ~ fyear))
  e2 <- as.data.frame(emmeans::emmeans(fit_off, ~ fyear, offset = log(2)))
  nd <- data.frame(fyear = levels(d$fyear), year = 2011)
  expect_equal(e0$emmean, predict(fit_off, newdata = nd)$est, tolerance = 1e-6)
  expect_equal(e2$emmean, e0$emmean + log(2), tolerance = 1e-6)

  # a smoother in only one delta component:
  fit_delta <- sdmTMB(
    list(density ~ fyear + s(depth_scaled, k = 5), density ~ fyear),
    data = d, spatial = "off", family = delta_gamma()
  )
  nd <- data.frame(fyear = levels(d$fyear), year = 2011, depth_scaled = 0.5)
  for (m in 1:2) {
    e <- as.data.frame(emmeans::emmeans(fit_delta, ~ fyear,
      at = list(depth_scaled = 0.5), model = m))
    p <- predict(fit_delta, newdata = nd, se_fit = TRUE, model = m)
    expect_equal(e$emmean, p$est, tolerance = 1e-6)
    expect_equal(e$SE, p$est_se, tolerance = 1e-6)
  }
})

test_that("emmeans errors on terms it can't represent", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  fit <- sdmTMB(present ~ as.factor(year) + breakpt(depth_scaled),
    data = pcod_2011, spatial = "off", family = binomial(), do_fit = FALSE)
  expect_error(emmeans::emmeans(fit, ~ year), "breakpt")
  fit <- sdmTMB(present ~ 0 + as.factor(year), time_varying = ~ 0 + depth_scaled,
    time = "year", data = pcod_2011, spatial = "off", family = binomial(),
    do_fit = FALSE)
  expect_error(emmeans::emmeans(fit, ~ year), "time_varying")
})
