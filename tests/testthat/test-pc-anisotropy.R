test_that("pc_anisotropy() parses its arguments", {
  expect_equal(pc_anisotropy(3), c(3, 0.05), ignore_attr = TRUE)
  expect_error(pc_anisotropy(0.5))
  expect_error(pc_anisotropy(3, 1.5))
})

test_that("Anisotropy PC prior is a normalized density of ln_H_input", {
  aniso <- rtmb_aniso_prior_inputs(list(ratio_lt = 3, ratio_prob = 0.05),
    list(matern_s = pc_matern(2, 1)))
  g <- seq(-8, 8, length.out = 401)
  h <- diff(g)[[1L]]
  grid <- expand.grid(x1 = g, x2 = g)
  cosh_v <- rtmb_aniso_cosh(list(grid$x1, grid$x2))
  dens <- exp(rtmb_pc_aniso(cosh_v, aniso$lambda_v))
  expect_equal(sum(dens) * h^2, 1, tolerance = 1e-3)
  # ratio of major to minor axis ranges is exp(|v|)
  expect_equal(sum(dens[exp(acosh(cosh_v)) > 3]) * h^2, 0.05,
    tolerance = 0.01)

  # marginal P(range < range_gt) = range_prob
  set.seed(1)
  f0 <- pc_aniso_f(1)
  f_v <- f0 - log(stats::runif(1e5)) / aniso$lambda_v
  range <- aniso$range_s * f_v / stats::rexp(1e5)
  expect_equal(mean(range < 2), 0.05, tolerance = 0.05)
})

test_that("Anisotropy PC prior reduces to the isotropic range prior", {
  aniso <- rtmb_aniso_prior_inputs(list(ratio_lt = 1.01,
    ratio_prob = 0.05), list(matern_s = pc_matern(2, 1)))
  expect_equal(aniso$range_s * pc_aniso_f(1), -log(0.05) * 2,
    tolerance = 1e-3)
  expect_equal(lambert_w0_log(log(2 * exp(2))), 2)
  expect_equal(lambert_w0_log(log(1e-8)), 1e-8, tolerance = 1e-6)
})

test_that("Anisotropy PC prior fits and regularizes the anisotropy", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = tweedie(),
      anisotropy = TRUE, ...)
  }
  ratio <- function(m) {
    e <- eigen(m$tmb_obj$report()$H)$values
    sqrt(e[[1L]] / e[[2L]])
  }
  m0 <- fit()
  m1 <- fit(priors = sdmTMBpriors(matern_s = pc_matern(10, 5),
    aniso = pc_anisotropy(3)))
  expect_equal(m1$model$convergence, 0L)
  expect_lt(max(abs(m1$gradients)), 1e-4)
  expect_lt(ratio(m1), ratio(m0))
  expect_lt(ratio(m1), 3)
  expect_true(sanity(m1, silent = TRUE)$all_ok)

  mb <- fit(priors = sdmTMBpriors(matern_s = pc_matern(10, 5),
    aniso = pc_anisotropy(3)), bayesian = TRUE)
  expect_equal(mb$model$convergence, 0L)
})

test_that("Anisotropy PC prior counts a shared H once", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = delta_gamma(),
      anisotropy = TRUE, priors = sdmTMBpriors(aniso = pc_anisotropy(3)),
      do_fit = FALSE, ...)
  }
  expect_equal(fit()$tmb_data$priors_aniso$components, 1L)
  separate <- fit(control = sdmTMBcontrol(map = list(ln_H_input = factor(1:4))))
  expect_equal(separate$tmb_data$priors_aniso$components, 1:2)
})

test_that("Anisotropy PC prior checks its requirements", {
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = tweedie(),
      priors = sdmTMBpriors(aniso = pc_anisotropy(3)), do_fit = FALSE, ...)
  }
  expect_error(fit(), "anisotropy = TRUE")
  expect_error(fit(anisotropy = TRUE,
    control = sdmTMBcontrol(backend = "tmb")), "rtmb")
})
