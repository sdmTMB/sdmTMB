simulate_nonlocal_regression_data <- function(
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

  vertex_grid <- expand.grid(
    vertex_id = seq_len(ncol(A_st)),
    year = seq_len(n_t)
  )
  loc <- as.data.frame(mesh$mesh$loc[, 1:2, drop = FALSE])
  names(loc) <- c("X", "Y")
  vertex_grid$X <- loc$X[vertex_grid$vertex_id]
  vertex_grid$Y <- loc$Y[vertex_grid$vertex_id]
  vertex_grid$x_space <- x_vertex[vertex_grid$vertex_id]
  x_time_vertex <- matrix(0, nrow = ncol(A_st), ncol = n_t)
  for (tt in seq_len(n_t)) {
    obs_t <- which(dat$year == tt)
    A_t <- A_st[obs_t, , drop = FALSE]
    numerator_t <- as.vector(Matrix::crossprod(A_t, dat$x_time[obs_t]))
    denominator_t <- as.vector(Matrix::crossprod(A_t, rep(1, length(obs_t))))
    keep_t <- denominator_t > 0
    x_time_vertex[keep_t, tt] <- numerator_t[keep_t] / denominator_t[keep_t]
  }
  vertex_grid$x_time <- x_time_vertex[cbind(vertex_grid$vertex_id, vertex_grid$year)]
  vertex_grid$vertex_id <- NULL

  dat$y_space <- 0.5 + 0.8 * x_diffused + rnorm(n_obs, sd = 0.2)
  dat$y_time <- -0.2 + 0.7 * dat$x_time + rnorm(n_obs, sd = 0.2)
  dat$y_st <- 0.1 + 0.6 * dat$x_st + rnorm(n_obs, sd = 0.2)

  list(data = dat, mesh = mesh, nonlocal_data = vertex_grid)
}

test_that("covariate diffusion regression estimates and logLik stay stable", {
  skip_on_cran()
  skip_on_ci()

  sim <- simulate_nonlocal_regression_data()
  dat <- sim$data
  mesh <- sim$mesh
  grid <- sim$nonlocal_data

  ctrl <- sdmTMBcontrol(newton_loops = 2L, getsd = FALSE)

  fit_space <- sdmTMB(
    y_space ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x_space),
    nonlocal_data = grid,
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
    nonlocal_formula = ~ time_lag(x_time),
    nonlocal_data = grid,
    control = ctrl
  )

  expect_equal(as.numeric(logLik(fit_space)), 154.3423, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit_time)), -240.2167, tolerance = 1e-4)

  get_beta <- function(fit) {
    b <- fit$tmb_obj$env$parList()$b_j
    names(b) <- colnames(fit$tmb_data$X_ij[[1]])
    b
  }
  beta_space <- get_beta(fit_space)
  beta_time <- get_beta(fit_time)

  expect_equal(beta_space[["(Intercept)"]], 0.4901745, tolerance = 1e-4)
  expect_equal(beta_space[["nl_diffusion_x_space"]], 0.6815684, tolerance = 1e-4)

  expect_equal(beta_time[["(Intercept)"]], -0.1992140, tolerance = 1e-4)
  expect_equal(beta_time[["nl_time_lag_x_time"]], 0.7224500, tolerance = 1e-4)

  rep_space <- fit_space$tmb_obj$report()
  rep_time <- fit_time$tmb_obj$report()

  expect_equal(rep_space$kappaS_nl[1], 7.894961005, tolerance = 1e-3)
  expect_equal(rep_time$kappaT_nl[1], 0.001682999, tolerance = 1e-3)
  expect_equal(rep_time$rhoT[1], 0.001680172, tolerance = 1e-3)
})
