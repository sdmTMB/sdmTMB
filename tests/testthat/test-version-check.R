test_that("reloading a model fit with another sdmTMB version warns", {
  withr::local_options(rlib_warning_verbosity = "verbose")
  fit <- sdmTMB(density ~ 1, data = pcod_2011, spatial = "off")
  fit$version <- package_version("1.2.0")
  reloaded <- unserialize(serialize(fit, NULL))
  expect_warning(capture.output(print(reloaded)), "fit with sdmTMB 1.2.0")

  fit$version <- utils::packageVersion("sdmTMB")
  reloaded <- unserialize(serialize(fit, NULL))
  expect_no_warning(capture.output(print(reloaded)))
})

test_that("fits from before sdmTMB 1.2.0 error early", {
  fit <- sdmTMB(density ~ 1, data = pcod_2011, spatial = "off")
  old <- fit
  old$version <- package_version("1.1.0")
  expect_error(predict(old), "before 1.2.0")
  expect_error(fitted(old), "before 1.2.0")
  expect_error(residuals(old), "before 1.2.0")
  expect_error(tidy(old), "before 1.2.0")
  expect_error(get_index(old, newdata = pcod_2011), "before 1.2.0")
  expect_error(cAIC(old), "before 1.2.0")
  expect_error(run_extra_optimization(old), "before 1.2.0")
  expect_error(update(old), "before 1.2.0")
  old$version <- NULL
  expect_error(print(old), "Refit the model")
})
