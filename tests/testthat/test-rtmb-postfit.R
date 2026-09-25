rtmb_postfit_fits <- function(backend) {
  d <- pcod_2011
  d$year_f <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  control <- sdmTMBcontrol(backend = backend)
  list(
    tweedie = sdmTMB(density ~ depth_scaled + (1 | year_f), data = d,
      mesh = mesh, time = "year", family = tweedie(),
      spatiotemporal = "iid", control = control),
    delta = sdmTMB(density ~ depth_scaled + (1 | year_f), data = d,
      mesh = mesh, time = "year", family = delta_gamma(type = "poisson-link"),
      spatiotemporal = "iid", control = control)
  )
}

test_that("RTMB cAIC, EDF, and MVN residuals match TMB", {
  skip_on_cran()
  cpp <- rtmb_postfit_fits("tmb")
  rt <- rtmb_postfit_fits("rtmb")
  for (name in names(cpp)) {
    expect_equal(cAIC(rt[[name]]), cAIC(cpp[[name]]), tolerance = 1e-6,
      info = name)
    expect_equal(cAIC(rt[[name]], what = "EDF"),
      cAIC(cpp[[name]], what = "EDF"), tolerance = 1e-5, info = name)
    # `env$MC()` draws the same random effects from the same seed.
    mvn <- function(fit) {
      set.seed(1)
      suppressMessages(residuals(fit, type = "mle-mvn"))
    }
    expect_equal(mvn(rt[[name]]), mvn(cpp[[name]]), tolerance = 1e-5,
      info = name)
  }
})

test_that("RTMB sdmTMB_simulate() output matches TMB and its moments", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(X = runif(300), Y = runif(300), a1 = rnorm(300),
    year = rep(1:6, each = 50))
  mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 0.1)
  sim <- function(backend, seed) {
    sdmTMB_simulate(~ 1 + a1, data = d, time = "year", mesh = mesh,
      family = gaussian(), range = 0.5, sigma_E = 0.1, phi = 0.1,
      sigma_O = 0.2, seed = seed, B = c(0.2, -0.4),
      control = sdmTMBcontrol(backend = backend))
  }
  cpp <- sim("tmb", 1)
  rt <- sim("rtmb", 1)
  expect_identical(names(rt), names(cpp))
  expect_equal(rt$mu - rt$omega_s - rt$epsilon_st,
    cpp$mu - cpp$omega_s - cpp$epsilon_st, tolerance = 1e-10)
  draws <- vapply(seq_len(100), function(i) sim("rtmb", i)$observed,
    numeric(nrow(d)))
  expect_equal(mean(draws), mean(0.2 - 0.4 * d$a1), tolerance = 0.03)
  expect_equal(mean(apply(draws, 1, sd)), sqrt(0.2^2 + 0.1^2 + 0.1^2),
    tolerance = 0.1)
})

test_that("RTMB project() matches TMB forecast moments", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  grid <- qcs_grid[seq(1, nrow(qcs_grid), by = 20), ]
  nd <- replicate_df(grid, "year", c(2011, 2013, 2015, 2017, 2019, 2021))
  moments <- function(backend) {
    fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
      family = tweedie(), time = "year", spatiotemporal = "ar1",
      control = sdmTMBcontrol(backend = backend))
    set.seed(1)
    p <- project(fit, newdata = nd, nsim = 300, silent = TRUE)
    rbind(
      mean = tapply(rowMeans(p$est), nd$year, mean),
      sd = tapply(apply(p$est, 1, sd), nd$year, mean)
    )
  }
  cpp <- moments("tmb")
  rt <- moments("rtmb")
  expect_equal(rt["mean", ], cpp["mean", ], tolerance = 0.05)
  expect_equal(rt["sd", ], cpp["sd", ], tolerance = 0.1)
  # The forecast years are more uncertain than the fitted ones.
  expect_true(all(rt["sd", c("2019", "2021")] > rt["sd", "2017"]))
})

# The tmbstan check is in tests/optional/, since tmbstan is not in Suggests.

test_that("RTMB mle-mcmc residuals use supplied linear predictors", {
  set.seed(105)
  d <- data.frame(x = seq(-1, 1, length.out = 40L))
  d$y <- rnorm(nrow(d), 0.4 + 0.7 * d$x, 0.5)
  fit <- function(backend) sdmTMB(y ~ x, data = d, spatial = "off",
    family = gaussian(), control = sdmTMBcontrol(backend = backend))
  cpp <- fit("tmb")
  rt <- fit("rtmb")
  sampled_eta <- predict(rt)$est + 0.1
  actual <- residuals(rt, type = "mle-mcmc", mcmc_samples = sampled_eta)
  expected <- (d$y - sampled_eta) / exp(get_pars(rt)[["ln_phi"]])
  expect_equal(actual, expected, tolerance = 1e-6)
  expect_equal(actual, residuals(cpp, type = "mle-mcmc",
    mcmc_samples = sampled_eta), tolerance = 1e-5)
})
