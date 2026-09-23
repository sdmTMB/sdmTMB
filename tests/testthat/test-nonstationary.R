# basic model fitting and prediction tests

test_that("removed epsilon_model options and vector sigma_E error", {
  d <- data.frame(X = runif(20), Y = runif(20), year = rep(1:2, each = 10))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.2)
  d$y <- rnorm(20)
  for (x in c("re", "trend-re")) {
    expect_error(sdmTMB(y ~ 1, data = d, mesh = mesh, time = "year",
      experimental = list(epsilon_model = x), do_fit = FALSE),
      regexp = "epsilon_model")
  }
  expect_error(sdmTMB_simulate(~ 1, data = d, mesh = mesh, time = "year",
    range = 0.5, sigma_E = c(0.1, 0.2), phi = 0.1, B = 0),
    regexp = "sigma_E")
})

test_that("Test that non-stationary model works without spatial field and epsilon trend works", {
  local_edition(2)
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "ar1",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.05852822, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(1.0534572, 1.0409799, 1.0285026, 1.0035480, 0.9785934, 0.9536388, 0.9286842, 0.9037296, 0.8787750), tolerance = 0.002)

  # fit non-stationary model - iid
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "iid",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.04915406, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(1.1262184, 1.1157395, 1.1052607, 1.0843029, 1.0633452, 1.0423874, 1.0214297, 1.0004719, 0.9795142), tolerance = 0.002)
})

test_that("Test that non-stationary model works with spatial field and epsilon trend works", {
  local_edition(2)
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="on",
    time = "year",
    spatiotemporal = "ar1",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.04435818, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(0.8882275, 0.8787710, 0.8693146, 0.8504016, 0.8314887, 0.8125758, 0.7936628, 0.7747499, 0.7558370), tolerance = 0.002)

  # fit non-stationary model - iid
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="on",
    time = "year",
    spatiotemporal = "iid",
    family = tweedie(link = "log"),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )
  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), -0.0457674, tolerance = 0.002)

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="log_sigma_E")]
  expect_equal(as.numeric(par), c(0.8916701, 0.8819133, 0.8721564, 0.8526426, 0.8331288, 0.8136150, 0.7941013, 0.7745875, 0.7550737), tolerance = 0.002)
})


test_that("Test that non-stationary model works with epsilon trend and delta model", {
  local_edition(2)
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  pcod$fyear <- as.factor(pcod$year)
  pcod$time <- pcod$year - min(pcod$year) + 1
  pcod$time = scale(pcod$year)
  fit <- sdmTMB(
    density ~ s(depth),
    data = pcod, mesh = mesh,
    spatial="off",
    time = "year",
    spatiotemporal = "ar1",
    family = delta_gamma(),
    experimental = list(epsilon_model = "trend", epsilon_predictor = "time"),
    control = sdmTMBcontrol(lower = list(b_epsilon = -1),
                            upper = list(b_epsilon = 1))
  )

  par <- fit$sd_report$value[which(names(fit$sd_report$value)=="b_epsilon")]
  expect_equal(as.numeric(par), c(-0.07908264, -0.09297464), tolerance = 0.002)

})
