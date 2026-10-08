test_that("sanity() runs", {
  skip_on_cran()

  fit <- sdmTMB(density ~ 1, time = "year",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie()
  )
  s <- sanity(fit)
  expect_true(all(unlist(s)))

  fit$gradients[1] <- 0.9
  expect_message(sanity(fit), regexp = "gradient")

  fit$gradients[1] <- 0.00001
  fit$bad_eig <- TRUE
  expect_message(sanity(fit), regexp = "model may not have converged")

  fit$bad_eig <- FALSE
  fit$model$convergence <- 1L
  expect_message(sanity(fit), regexp = "did not converge")
  fit$model$convergence <- 0L

  fit$pos_def_hessian <- FALSE
  expect_message(sanity(fit), regexp = "Non-positive-definite")

  fit <- sdmTMB(density ~ 1, time = "year",
    spatiotemporal = "off", spatial = "off",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  s <- sanity(fit)
  expect_true(all(unlist(s)))

  set.seed(1)
  pcod_2011$depth_fake <- pcod_2011$depth + runif(nrow(pcod_2011), -0.0001, 0.0001)
  suppressWarnings(
    fit <- sdmTMB(density ~ depth + depth_fake,
      data = pcod_2011, mesh = pcod_mesh_2011, spatial = "off",
      family = delta_gamma()
    )
  )
  expect_message(s <- sanity(fit), regexp = "may be large")
  # expect_message(s <- sanity(fit, se_ratio = 2), regexp = "2x")
  expect_false(s$se_magnitude_ok)
  expect_false(s$all_ok)

  expect_false(sanity(NA))
  expect_false(sanity(NULL))

  x <- try(stop())
  expect_false(sanity(x))
})

test_that("sanity() doesn't warn for non-spatial fits without a mesh", {
  skip_on_cran()
  fit <- sdmTMB(present ~ depth_scaled, data = pcod_2011, spatial = "off",
    family = binomial())
  expect_no_warning(s <- sanity(fit, silent = TRUE))
  expect_true(s$range_ok)
})
