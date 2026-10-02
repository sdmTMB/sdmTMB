test_that("removed epsilon_model options and vector sigma_E error", {
  d <- data.frame(X = runif(20), Y = runif(20), year = rep(1:2, each = 10))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.2)
  d$y <- rnorm(20)
  d$cov <- d$year
  removed <- list(
    list(epsilon_model = "re"),
    list(epsilon_model = "trend-re"),
    list(epsilon_model = "trend", epsilon_predictor = "cov"),
    list(epsilon_predictor = "cov")
  )
  for (x in removed) {
    expect_error(sdmTMB(y ~ 1, data = d, mesh = mesh, time = "year",
      experimental = x, do_fit = FALSE),
      regexp = "have been removed")
  }
  expect_error(sdmTMB_simulate(~ 1, data = d, mesh = mesh, time = "year",
    range = 0.5, sigma_E = c(0.1, 0.2), phi = 0.1, B = 0),
    regexp = "sigma_E")
})
