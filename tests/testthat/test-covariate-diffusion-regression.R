simulate_covariate_diffusion_regression_data <- function(
  n_obs = 700,
  n_t = 6,
  cutoff = 0.08,
  n_hotspots = 40,
  hotspot_sd = 0.04,
  seed = 321,
  kappaS = 6
) {
  set.seed(seed)

  dat <- data.frame(
    X = runif(n_obs),
    Y = runif(n_obs),
    year = sample(seq_len(n_t), size = n_obs, replace = TRUE)
  )

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  centers <- cbind(runif(n_hotspots), runif(n_hotspots))
  weights <- rnorm(n_hotspots)
  sq_dist <- outer(dat$X, centers[, 1], "-")^2 +
    outer(dat$Y, centers[, 2], "-")^2
  hotspot_field <- as.numeric(exp(-sq_dist / (2 * hotspot_sd^2)) %*% weights)

  A_st <- mesh$A_st
  M0 <- mesh$spde$c0
  M1 <- mesh$spde$g1

  # Vignette-style spatial diffusion simulator.
  x_fine <- as.numeric(scale(hotspot_field + rnorm(n_obs, sd = 0.25)))
  numerator <- as.vector(Matrix::crossprod(A_st, x_fine))
  denominator <- as.vector(Matrix::crossprod(A_st, rep(1, n_obs)))
  x_vertex <- numeric(ncol(A_st))
  keep <- denominator > 0
  x_vertex[keep] <- numerator[keep] / denominator[keep]
  diffusion_matrix <- M0 + (1 / kappaS^2) * M1
  x_vertex_diffused <- as.numeric(Matrix::solve(diffusion_matrix, M0 %*% x_vertex))
  x_diffused <- as.numeric(A_st %*% x_vertex_diffused)

  t_scaled <- as.numeric(scale(dat$year))

  dat$x_space <- x_fine
  dat$x_time <- as.numeric(scale(sin(dat$year) + 0.3 * dat$year + rnorm(n_obs, sd = 0.3)))
  dat$x_st <- as.numeric(scale(hotspot_field * t_scaled + rnorm(n_obs, sd = 0.25)))

  dat$y_space <- 0.5 + 0.8 * x_diffused + rnorm(n_obs, sd = 0.2)
  dat$y_time <- -0.2 + 0.7 * dat$x_time + rnorm(n_obs, sd = 0.2)
  dat$y_st <- 0.1 + 0.6 * dat$x_st + rnorm(n_obs, sd = 0.2)

  list(data = dat, mesh = mesh)
}

test_that("covariate diffusion regression estimates and logLik stay stable", {
  skip_on_cran()

  sim <- simulate_covariate_diffusion_regression_data()
  dat <- sim$data
  mesh <- sim$mesh

  ctrl <- sdmTMBcontrol(newton_loops = 0, getsd = FALSE)

  fit_space <- sdmTMB(
    y_space ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ space(x_space),
    control = ctrl
  )

  fit_time <- sdmTMB(
    y_time ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ time(x_time),
    control = ctrl
  )

  fit_st <- suppressWarnings(sdmTMB(
    y_st ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ spacetime(x_st),
    control = ctrl
  ))

  expect_equal(as.numeric(logLik(fit_space)), 92.0763, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit_time)), -240.2167, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit_st)), -278.3500, tolerance = 1e-4)

  get_beta <- function(fit) {
    b <- fit$tmb_obj$env$parList()$b_j
    names(b) <- colnames(fit$tmb_data$X_ij[[1]])
    b
  }
  beta_space <- get_beta(fit_space)
  beta_time <- get_beta(fit_time)
  beta_st <- get_beta(fit_st)

  expect_equal(beta_space[["(Intercept)"]], 0.4815081, tolerance = 1e-4)
  expect_equal(beta_space[["cov_diff_space_x_space"]], 0.3733943, tolerance = 1e-4)

  expect_equal(beta_time[["(Intercept)"]], -0.1992140, tolerance = 1e-4)
  expect_equal(beta_time[["cov_diff_time_x_time"]], 0.7224500, tolerance = 1e-4)

  expect_equal(beta_st[["(Intercept)"]], 0.1152084, tolerance = 1e-4)
  expect_equal(beta_st[["cov_diff_spacetime_x_st"]], 0.7260176, tolerance = 1e-4)

  rep_space <- fit_space$tmb_obj$report()
  rep_time <- fit_time$tmb_obj$report()
  rep_st <- fit_st$tmb_obj$report()

  expect_equal(rep_space$kappaS_dl[1], 10.206683946, tolerance = 1e-4)
  expect_equal(rep_time$kappaT_dl[1], 0.001683076, tolerance = 1e-4)
  expect_equal(rep_time$rhoT[1], 0.001680248, tolerance = 1e-4)
  expect_equal(rep_st$kappaST_dl[1], 25514.0, tolerance = 1e-4)
})
