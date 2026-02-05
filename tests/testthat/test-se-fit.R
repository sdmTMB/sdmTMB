# Related to issue #496
test_that("Delta model predictions are consistent with and without se_fit (re_form = NA)", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011,
    spatial = "off",
    family = delta_gamma()
  )

  # Test with re_form = NA
  p_no_se <- predict(fit, re_form = NA, newdata = pcod_2011)
  p_with_se <- predict(fit, se_fit = TRUE, re_form = NA, newdata = pcod_2011)

  # Component estimates should be identical
  expect_equal(p_no_se$est1, p_with_se$est1)
  expect_equal(p_no_se$est2, p_with_se$est2)

  # Combined estimate should be identical (convert to numeric to handle array vs vector)
  expect_equal(as.numeric(p_no_se$est), as.numeric(p_with_se$est))

  # Combined estimate should match manual calculation: log(plogis(est1) * exp(est2))
  manual_est <- log(plogis(p_no_se$est1) * exp(p_no_se$est2))
  expect_equal(as.numeric(p_no_se$est), manual_est, tolerance = 1e-4)
  expect_equal(as.numeric(p_with_se$est), manual_est, tolerance = 1e-4)
})

test_that("Delta model component predictions consistent with and without se_fit (default re_form)", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011,
    spatial = "off",
    family = delta_gamma()
  )

  # Test with default re_form
  p_no_se <- predict(fit, newdata = pcod_2011)
  p_with_se <- predict(fit, se_fit = TRUE, newdata = pcod_2011)

  # Component estimates should be identical
  expect_equal(p_no_se$est1, p_with_se$est1)
  expect_equal(p_no_se$est2, p_with_se$est2)

  # When se_fit = TRUE with default re_form, est should match manual calculation
  manual_est <- log(plogis(p_with_se$est1) * exp(p_with_se$est2))
  expect_equal(as.numeric(p_with_se$est), manual_est, tolerance = 1e-4)
})

test_that("Delta-lognormal predictions are consistent with and without se_fit", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod_2011,
    spatial = "off",
    family = delta_lognormal()
  )

  # Test with re_form = NA
  p_no_se <- predict(fit, re_form = NA, newdata = pcod_2011)
  p_with_se <- predict(fit, se_fit = TRUE, re_form = NA, newdata = pcod_2011)

  # Component estimates should be identical
  expect_equal(p_no_se$est1, p_with_se$est1)
  expect_equal(p_no_se$est2, p_with_se$est2)

  # Combined estimate should be identical
  expect_equal(as.numeric(p_no_se$est), as.numeric(p_with_se$est))

  # Combined estimate should match manual calculation
  manual_est <- log(plogis(p_no_se$est1) * exp(p_no_se$est2))
  expect_equal(as.numeric(p_no_se$est), manual_est, tolerance = 1e-4)
  expect_equal(as.numeric(p_with_se$est), manual_est, tolerance = 1e-4)
})
