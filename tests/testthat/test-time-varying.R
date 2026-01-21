test_that("AR1 time-varying works", {
  skip_on_cran()
  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(4000), Y = runif(4000),
    year = rep(seq_len(800), each = 5)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)

  sigma_V_hats <- numeric(10L)
  rho_hats <- numeric(10L)

  set.seed(1234)
  for (j in seq_along(rho_hats)) {
    print(j)
    # https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html
    ar1_scale <- 0.8
    rho <- 0.7
    sigma <- sqrt(1 - rho^2)
    devs <- rnorm(length(unique(predictor_dat$year))) * ar1_scale
    x0 <- rnorm(1) * ar1_scale
    B <- numeric(length(unique(predictor_dat$year)))
    B[1] <- rho * x0 + sigma * devs[1]
    for (i in seq(2, length(B))) {
      B[i] <- rho * B[i-1] + sigma * devs[i]
    }

    sim_dat <- sdmTMB_simulate(
      formula = ~ 0 + as.factor(year),
      data = predictor_dat,
      time = NULL,
      mesh = mesh,
      family = gaussian(),
      range = 0.5,
      sigma_E = NULL,
      phi = 0.1,
      sigma_O = 0,
      seed = j,
      B = B
    )

    sim_dat$year <- predictor_dat$year
    m <- sdmTMB(observed ~ 1, data = sim_dat, mesh = mesh,
      time = "year", spatial = 'off', time_varying_type = "ar1",
      spatiotemporal = 'off', time_varying = ~ 1)

    s <- as.list(m$sd_report, "Estimate")

    sigma_V_hats[j] <- exp(s$ln_tau_V)[1,1]
    m121 <- function(x) 2 * plogis(x) - 1
    rho_hats[j] <- m121(s$rho_time_unscaled)[1,1]
  }

  plot(s$b_rw_t[,,1], B)
  abline(0, 1)
  expect_gt(cor(B, s$b_rw_t[,,1]), 0.99)

  expect_equal(round(mean(sigma_V_hats), 2L), ar1_scale)
  expect_equal(round(mean(rho_hats), 2L), rho)

  ss <- simulate(m)
  sim_dat$obs2 <- ss[,1]
  m <- sdmTMB(obs2 ~ 0, data = sim_dat, mesh = mesh,
    time = "year", spatial = 'off', time_varying_type = "ar1",
    spatiotemporal = 'off', time_varying = ~ 1)

  s <- as.list(m$sd_report, "Estimate")
  expect_equal(dim(ss), c(4000L, 1L))
  expect_gt(exp(s$ln_tau_V)[1,1], 0.75)
  expect_gt(m121(s$rho_time_unscaled)[1,1], 0.65)

  # test that tidy works
  fit <- sdmTMB(density ~ 1, time = "year",
                time_varying = ~ 1,
                time_varying_type = "ar1",
                data = pcod_2011, spatial="off",
                spatiotemporal = "off",
                family = tweedie())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(3L, 5L))
  expect_equal(s$term, c("phi","rho_time","tweedie_p"))

  # test that tidy works -- RW
  fit <- sdmTMB(density ~ 1, time = "year",
                time_varying = ~ -1 + depth_scaled,
                data = pcod_2011, spatial="off",
                spatiotemporal = "off",
                family = tweedie())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(2L, 5L))
  expect_equal(s$term, c("phi","tweedie_p"))

  # test that tidy works -- AR1
  fit <- sdmTMB(density ~ 1, time = "year",
                time_varying = ~ -1 + depth_scaled,
                data = pcod_2011, spatial="off",
                time_varying_type = "ar1",
                spatiotemporal = "off",
                family = tweedie())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(3L, 5L))
  expect_equal(s$term, c("phi", "rho_time", "tweedie_p"))

  # test that tidy works -- no time varying
  fit <- sdmTMB(density ~ as.factor(year), time = "year",
                data = pcod_2011, spatial="off",
                spatiotemporal = "off",
                family = tweedie())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(2L, 5L))
  expect_equal(s$term, c("phi","tweedie_p"))

  # test that tidy works with 2 parameters
  fit <- sdmTMB(density ~ 1, time = "year",
                time_varying = ~ -1+depth_scaled + I(depth_scaled^2),
                time_varying_type = "ar1",
                data = pcod_2011, spatial="off",
                spatiotemporal = "off",
                family = tweedie())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(4L, 5L))
  expect_equal(s$term, c("phi","rho_time","rho_time","tweedie_p"))

  # test that tidy works with delta
  fit <- sdmTMB(density ~ 1, time = "year",
                time_varying = ~ -1+depth_scaled + I(depth_scaled^2),
                time_varying_type = "ar1",
                data = pcod_2011, spatial="off",
                spatiotemporal = "off",
                family = delta_gamma())
  s <- tidy(fit, "ran_pars")
  expect_equal(dim(s), c(2L, 6L))
  expect_equal(s$term, c("rho_time","rho_time"))

  s <- tidy(fit, "ran_pars", model = 1)
  expect_equal(dim(s), c(2L, 6L))
  expect_equal(s$term, c("rho_time","rho_time"))
  expect_equal(s$estimate, c(1.0, 0.9582931))

  s <- tidy(fit, "ran_pars", model = 2)
  expect_equal(s$term, c("phi", "rho_time","rho_time"))
  expect_equal(s$estimate, c(0.652, 1.0, 0.868), tolerance = 0.001)

})


test_that("RW with mean-zero (rw0) time-varying works", {
  skip_on_cran()
  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(4000), Y = runif(4000),
    year = rep(seq_len(800), each = 5)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)

  sigma_V_hats <- numeric(12L)
  sigma_true <- 0.3

  set.seed(1234)
  for (j in seq_along(sigma_V_hats)) {
    print(j)
    devs <- rnorm(length(unique(predictor_dat$year)))
    B <- numeric(length(unique(predictor_dat$year)))
    B[1] <- sigma_true * devs[1]
    for (i in seq(2, length(B))) {
      B[i] <- B[i-1] + sigma_true * devs[i]
    }

    sim_dat <- sdmTMB_simulate(
      formula = ~ 0 + as.factor(year),
      data = predictor_dat,
      time = NULL,
      mesh = mesh,
      family = gaussian(),
      range = 0.5,
      sigma_E = NULL,
      phi = 0.1,
      sigma_O = 0,
      seed = j,
      B = B
    )

    sim_dat$year <- predictor_dat$year
    m <- sdmTMB(observed ~ 1, data = sim_dat, mesh = mesh,
      time = "year", spatial = 'off', time_varying_type = "rw0",
      spatiotemporal = 'off', time_varying = ~ 1)

    s <- as.list(m$sd_report, "Estimate")

    sigma_V_hats[j] <- exp(s$ln_tau_V)[1,1]
  }

  plot(s$b_rw_t[,,1], B)
  abline(0, 1)
  expect_gt(cor(B, s$b_rw_t[,,1]), 0.99)
  expect_equal(round(mean(sigma_V_hats), 2L), sigma_true)
})
