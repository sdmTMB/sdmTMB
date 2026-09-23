# Optional RTMB checks for packages that sdmTMB doesn't declare as Suggests
# (tmbstan/Stan and INLAspacetime/INLAtools). This directory is excluded from
# the built package, so R CMD check doesn't run it. Run with
# `make test-optional`, or with the package loaded via devtools::load_all():
#   testthat::test_dir("tests/optional")

test_that("tmbstan samples an RTMB fit and MCMC prediction uses the draws", {
  skip_if_not_installed("tmbstan")
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 30)
  fit <- function(backend) {
    sdmTMB(present ~ depth_scaled, data = pcod_2011, mesh = mesh,
      family = binomial(),
      priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 15, sigma_lt = 5),
        b = normal(c(0, 0), c(5, 5))),
      bayesian = TRUE, control = sdmTMBcontrol(backend = backend))
  }
  cpp <- fit("tmb")
  rt <- fit("rtmb")
  expect_identical(names(rt$tmb_obj$env$last.par.best),
    names(cpp$tmb_obj$env$last.par.best))
  # A short smoke run; convergence warnings are expected.
  stan <- suppressWarnings(tmbstan::tmbstan(rt$tmb_obj, chains = 1,
    iter = 100, seed = 1, refresh = 0))
  draws <- t(as.matrix(stan))
  draws <- draws[rownames(draws) != "lp__", , drop = FALSE]
  expect_identical(nrow(draws), length(rt$tmb_obj$env$last.par.best))
  expect_true(all(is.finite(draws)))
  pred <- predict(rt, mcmc_samples = draws[, 1:3])
  for (i in 1:3) {
    eta <- rt$tmb_obj$report(draws[, i])$eta_i[, 1L]
    expect_equal(pred[, i], eta, ignore_attr = TRUE, tolerance = 1e-10)
  }
})

# The unit-SD barrier precision matches INLAspacetime's compiled model.
# Passing `libpath` uses INLAspacetime's own library, so INLA is not needed.
test_that("RTMB barrier precision matches INLAspacetime", {
  skip_if_not_installed("sdmTMBextra")
  skip_if_not_installed("sf")
  skip_if_not_installed("INLAspacetime")
  skip_if_not_installed("INLAtools")
  set.seed(72)
  d <- data.frame(x = runif(90, 0, 5), y = runif(90, 0, 5), z = rnorm(90))
  d$gaussian <- d$z + sin(d$x) + rnorm(90, sd = 0.3)
  mesh <- make_mesh(d, c("x", "y"), cutoff = 0.7)
  barrier <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_polygon(list(
    matrix(c(2.3, -0.2, 2.7, -0.2, 2.7, 5.2, 2.3, 5.2, 2.3, -0.2),
      ncol = 2, byrow = TRUE)))))
  mesh <- suppressMessages(suppressWarnings(
    sdmTMBextra::add_barrier_mesh(mesh, barrier, range_fraction = 0.2,
      plot = FALSE)))
  fit <- sdmTMB(gaussian ~ z, data = d, mesh = mesh, do_fit = FALSE)
  model <- INLAspacetime::barrierModel.define(mesh$mesh,
    mesh$barrier_triangles, prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5), range.fraction = 0.2, useINLAprecomp = FALSE,
    libpath = INLAtools::cgeneric_shlib_path(package = "INLAspacetime",
      useINLAprecomp = FALSE))
  inputs <- rtmb_precision_inputs(fit$tmb_data)
  for (range in c(1, 2.5)) {
    Q <- rtmb_precision(inputs, list(range = matrix(range, 2L, 1L)), 1L, 1L)
    # theta = (log range, log sigma)
    Q_inla <- INLAtools::cgeneric_Q(model, theta = c(log(range), 0))
    expect_equal(as.matrix(Q), as.matrix(Q_inla), tolerance = 1e-8,
      info = range)
  }
})
