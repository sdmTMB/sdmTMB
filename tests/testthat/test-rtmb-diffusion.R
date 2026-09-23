rtmb_diffusion_data <- function() {
  set.seed(74)
  d <- expand.grid(X = 1:6, Y = 1:2, year = 1:5)
  d$x1 <- as.numeric(scale(sin(d$year / 2) + d$X / 6))
  d$x2 <- as.numeric(scale(cos(d$year / 3) + d$Y / 2))
  d$y <- 0.2 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(nrow(d), sd = 0.2)
  d$response <- ifelse(d$x2 > 0, exp(d$y), 0)
  d
}

test_that("RTMB covariate diffusion matches every TMB report and sdreport row", {
  d <- rtmb_diffusion_data()
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.5)
  grid <- make_nl_covariate_grid(mesh, 1:5, c("x1", "x2"))
  newdata <- d[c(1, 20, 60), ]
  cases <- list(
    spatial_time = list(y ~ x1 + x2, family = gaussian(),
      nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
      spatiotemporal = "iid"),
    joint = list(y ~ x1, family = gaussian(),
      nonlocal_formula = ~ diffusion(x1) + time_lag(x1), spatial = "off"),
    delta = list(response ~ x1, family = delta_gamma(),
      nonlocal_formula = ~ diffusion(x1) + time_lag(x2), spatial = "off")
  )
  for (name in names(cases)) {
    args <- c(cases[[name]], list(data = d, mesh = mesh, time = "year",
      nonlocal_data = grid, do_fit = FALSE))
    fit <- do.call(sdmTMB, args)
    expect_rtmb_fit_data_matches(fit, newdata, info = name)
  }
})
