test_that("RSR coefficients are computed and differ from b_j", {
  skip_on_cran()
  skip_on_ci()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod,
    mesh = mesh,
    family = tweedie(),
    control = sdmTMBcontrol(get_rsr = TRUE)
  )

  s <- as.list(fit$sd_report, "Estimate")
  sr <- as.list(fit$sd_report, "Estimate", report = TRUE)

  expect_length(sr$b_j_prime, 2L)
  expect_false(identical(s$b_j, sr$b_j_prime))
  expect_equal(s$b_j, c(2.62917, -0.74587), tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(sr$b_j_prime, c(2.6291, -0.7283), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("tidy() works with effects = 'rsr'", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod,
    mesh = mesh,
    family = tweedie(),
    control = sdmTMBcontrol(get_rsr = TRUE)
  )

  rsr_tidy <- tidy(fit, effects = "rsr")
  expect_s3_class(rsr_tidy, "tbl_df")
  expect_equal(nrow(rsr_tidy), 2L)
  expect_equal(rsr_tidy$term, c("(Intercept)", "depth_scaled"))
  expect_true(all(c("term", "estimate", "std.error", "conf.low", "conf.high") %in% names(rsr_tidy)))

  regular_tidy <- tidy(fit, effects = "fixed")
  expect_false(identical(rsr_tidy$estimate, regular_tidy$estimate))

  rsr_tidy_no_ci <- tidy(fit, effects = "rsr", conf.int = FALSE)
  expect_false("conf.low" %in% names(rsr_tidy_no_ci))
  expect_false("conf.high" %in% names(rsr_tidy_no_ci))

  rsr_tidy_exp <- tidy(fit, effects = "rsr", exponentiate = TRUE)
  expect_false("std.error" %in% names(rsr_tidy_exp))
  expect_true(all(rsr_tidy_exp$estimate > 0))
})

test_that("tidy() RSR aborts if fit without get_rsr = TRUE", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod,
    mesh = mesh,
    family = tweedie()
  )
  expect_error(tidy(fit, effects = "rsr"), regexp = "get_rsr")
})

test_that("tidy() RSR works on delta models", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = pcod,
    mesh = mesh,
    family = delta_gamma(),
    control = sdmTMBcontrol(get_rsr = TRUE)
  )

  rsr_m1 <- tidy(fit, effects = "rsr", model = 1)
  expect_s3_class(rsr_m1, "tbl_df")
  expect_equal(nrow(rsr_m1), 2L)
  expect_equal(rsr_m1$term, c("(Intercept)", "depth_scaled"))

  rsr_m2 <- tidy(fit, effects = "rsr", model = 2)
  expect_s3_class(rsr_m2, "tbl_df")
  expect_equal(nrow(rsr_m2), 2L)
  expect_equal(rsr_m2$term, c("(Intercept)", "depth_scaled"))
})
