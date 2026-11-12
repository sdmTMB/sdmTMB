edf_by_term <- function(fit) {
  as.numeric(tapply(fit$smooth_edf[, 1], fit$smoothers$zs_basis, sum))
}

test_that("smoother EDF matches mgcv", {
  skip_on_cran()
  set.seed(1)
  n <- 400
  d <- data.frame(x = runif(n), z = runif(n), f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- sin(2 * pi * d$x) + d$z^2 + as.numeric(d$f) * d$x + rnorm(n, 0, 0.5)
  d$count <- rpois(n, exp(1 + sin(2 * pi * d$x)))

  check <- function(formula, family = gaussian(), reml = FALSE, tolerance = 1e-3) {
    g <- mgcv::gam(formula, data = d, family = family,
      method = if (reml) "REML" else "ML")
    fit <- sdmTMB(formula, data = d, spatial = "off", family = family, reml = reml)
    expect_equal(edf_by_term(fit), as.numeric(summary(g)$edf), tolerance = tolerance)
  }
  check(y ~ s(x) + s(z))
  check(y ~ s(x) + s(z), reml = TRUE)
  check(y ~ f + s(x, by = f))
  check(y ~ t2(x, z))
  check(y ~ s(x, bs = "cc") + s(z, bs = "cs"))
  check(count ~ s(x) + s(z, k = 5), family = poisson())

  # finite-difference steps must scale with the design, not the coefficient
  d$z_big <- d$z * 1e4
  g <- mgcv::gam(count ~ s(x) + z_big, data = d, family = poisson(), method = "ML")
  # the optimizer warns about the poorly scaled covariate; the EDF is still right
  fit <- suppressWarnings(sdmTMB(count ~ s(x) + z_big, data = d,
    spatial = "off", family = poisson()))
  expect_equal(edf_by_term(fit), as.numeric(summary(g)$edf), tolerance = 1e-3)
})

test_that("smoother EDF is stored and shown in print() and tidy()", {
  skip_on_cran()
  fit <- sdmTMB(density ~ s(depth) + s(year, k = 5), data = pcod,
    spatial = "off", family = tweedie())
  expect_equal(dim(fit$smooth_edf), c(2L, 1L))
  expect_true(all(fit$smooth_edf > 1 & fit$smooth_edf < c(9, 4)))

  td <- tidy(fit, "ran_pars")
  sm <- grepl("^sd__s\\(", td$term)
  expect_equal(td$edf[sm], as.numeric(fit$smooth_edf[, 1]))
  expect_true(all(is.na(td$edf[!sm])))

  out <- capture.output(print(fit))
  expect_true(any(grepl("EDF", out)))

  fit2 <- sdmTMB(density ~ depth_scaled, data = pcod_2011, spatial = "off",
    family = tweedie())
  expect_null(fit2$smooth_edf)
  expect_false("edf" %in% names(tidy(fit2, "ran_pars")))
})

test_that("smoother EDF works for delta models", {
  skip_on_cran()
  fit <- sdmTMB(density ~ s(depth), data = pcod_2011, spatial = "off",
    family = delta_gamma())
  expect_equal(dim(fit$smooth_edf), c(1L, 2L))
  expect_equal(tidy(fit, "ran_pars", model = 2)$edf[2], fit$smooth_edf[1, 2])
})
