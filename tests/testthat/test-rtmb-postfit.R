test_that("RTMB sdmTMB_simulate() output matches its moments", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(X = runif(300), Y = runif(300), a1 = rnorm(300),
    year = rep(1:6, each = 50))
  mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 0.1)
  sim <- function(seed) {
    sdmTMB_simulate(~ 1 + a1, data = d, time = "year", mesh = mesh,
      family = gaussian(), range = 0.5, sigma_E = 0.1, phi = 0.1,
      sigma_O = 0.2, seed = seed, B = c(0.2, -0.4))
  }
  rt <- sim(1)
  expect_equal(rt$mu - rt$omega_s - rt$epsilon_st, 0.2 - 0.4 * d$a1,
    tolerance = 1e-10)
  draws <- vapply(seq_len(100), function(i) sim(i)$observed,
    numeric(nrow(d)))
  expect_equal(mean(draws), mean(0.2 - 0.4 * d$a1), tolerance = 0.03)
  expect_equal(mean(apply(draws, 1, sd)), sqrt(0.2^2 + 0.1^2 + 0.1^2),
    tolerance = 0.1)
})

test_that("RTMB project() forecasts grow more uncertain", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  grid <- qcs_grid[seq(1, nrow(qcs_grid), by = 20), ]
  nd <- replicate_df(grid, "year", c(2011, 2013, 2015, 2017, 2019, 2021))
  moments <- function() {
    fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
      family = tweedie(), time = "year", spatiotemporal = "ar1")
    set.seed(1)
    p <- project(fit, newdata = nd, nsim = 300, silent = TRUE)
    rbind(
      mean = tapply(rowMeans(p$est), nd$year, mean),
      sd = tapply(apply(p$est, 1, sd), nd$year, mean)
    )
  }
  rt <- moments()
  # The forecast years are more uncertain than the fitted ones.
  expect_true(all(rt["sd", c("2019", "2021")] > rt["sd", "2017"]))
})

# The tmbstan check is in tests/optional/, since tmbstan is not in Suggests.

test_that("RTMB mle-mcmc residuals use supplied linear predictors", {
  set.seed(105)
  d <- data.frame(x = seq(-1, 1, length.out = 40L))
  d$y <- rnorm(nrow(d), 0.4 + 0.7 * d$x, 0.5)
  rt <- sdmTMB(y ~ x, data = d, spatial = "off", family = gaussian())
  sampled_eta <- predict(rt)$est + 0.1
  actual <- residuals(rt, type = "mle-mcmc", mcmc_samples = sampled_eta)
  expected <- (d$y - sampled_eta) / exp(get_pars(rt)[["ln_phi"]])
  expect_equal(actual, expected, tolerance = 1e-6)
})
