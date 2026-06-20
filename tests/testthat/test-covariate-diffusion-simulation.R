make_dl_sim_data <- function() {
  grid <- expand.grid(
    X = seq(0, 1, length.out = 4),
    Y = seq(0, 1, length.out = 4),
    year = 1:3
  )
  grid$x_s <- sin(grid$X * pi) + cos(grid$Y * pi)
  grid$x_t <- grid$year
  grid$x_st <- grid$x_s + grid$year / 10
  grid
}

make_dl_sim_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.6)
}

test_that("simulate_new supports space covariate diffusion", {
  skip_on_cran()

  dat <- make_dl_sim_data()
  mesh <- make_dl_sim_mesh(dat)

  sim <- simulate_new(
    formula = ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    range = 0.5,
    sigma_O = 0,
    phi = 0.1,
    B = c(0.2, 0.7),
    covariate_diffusion = ~ space(x_s),
    diffusion_kappaS = 1.3,
    seed = 1
  )

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), nrow(dat))
  expect_true(all(c("observed", "eta") %in% names(sim)))
  expect_true("diffusion_truth_space_x_s" %in% names(sim))
  expect_false(any(grepl("^cov_diff_", names(sim))))
})

test_that("simulate_new supports time covariate diffusion", {
  skip_on_cran()

  dat <- make_dl_sim_data()
  mesh <- make_dl_sim_mesh(dat)

  sim <- simulate_new(
    formula = ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    range = 0.5,
    sigma_O = 0,
    phi = 0.1,
    B = c(0.2, 0.7),
    covariate_diffusion = ~ time(x_t),
    diffusion_rhoT = 0.4,
    seed = 2
  )

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), nrow(dat))
  expect_true(all(c("observed", "eta") %in% names(sim)))
  expect_true("diffusion_truth_time_x_t" %in% names(sim))
  expect_false(any(grepl("^cov_diff_", names(sim))))
})

test_that("simulate_new supports combined covariate diffusion terms", {
  skip_on_cran()

  dat <- make_dl_sim_data()
  mesh <- make_dl_sim_mesh(dat)

  sim <- simulate_new(
    formula = ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    range = 0.5,
    sigma_O = 0,
    phi = 0.1,
    B = c(0.2, 0.5, -0.3),
    covariate_diffusion = ~ space(x_s) + time(x_t),
    diffusion_kappaS = 1.3,
    diffusion_rhoT = 0.4,
    seed = 3
  )

  expect_s3_class(sim, "data.frame")
  expect_true(all(c(
    "diffusion_truth_space_x_s",
    "diffusion_truth_time_x_t"
  ) %in% names(sim)))
  expect_false(any(grepl("^cov_diff_", names(sim))))
})

test_that("simulate_new errors for missing covariate diffusion parameters", {
  dat <- make_dl_sim_data()
  mesh <- make_dl_sim_mesh(dat)

  expect_error(
    simulate_new(
      formula = ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      family = gaussian(),
      spatial = "off",
      spatiotemporal = "off",
      range = 0.5,
      sigma_O = 0,
      phi = 0.1,
      B = c(0.2, 0.7),
      covariate_diffusion = ~ space(x_s)
    ),
    "diffusion_kappaS"
  )

  expect_error(
    simulate_new(
      formula = ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      family = gaussian(),
      spatial = "off",
      spatiotemporal = "off",
      range = 0.5,
      sigma_O = 0,
      phi = 0.1,
      B = c(0.2, 0.7),
      covariate_diffusion = ~ time(x_t)
    ),
    "diffusion_rhoT"
  )

  expect_error(
    simulate_new(
      formula = ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      family = gaussian(),
      spatial = "off",
      spatiotemporal = "off",
      range = 0.5,
      sigma_O = 0,
      phi = 0.1,
      B = c(0.2, 0.7),
      covariate_diffusion = ~ spacetime(x_st),
      diffusion_kappaS = 1.3
    ),
    "Unsupported wrapper"
  )
})

test_that("simulate_new B length validation includes diffusion coefficient columns", {
  dat <- make_dl_sim_data()
  mesh <- make_dl_sim_mesh(dat)

  expect_error(
    simulate_new(
      formula = ~ 1,
      data = dat,
      mesh = mesh,
      time = "year",
      family = gaussian(),
      spatial = "off",
      spatiotemporal = "off",
      range = 0.5,
      sigma_O = 0,
      phi = 0.1,
      B = 0.2,
      covariate_diffusion = ~ space(x_s),
      diffusion_kappaS = 1.3
    ),
    "covariate_diffusion"
  )
})

test_that("simulated spatial diffusion recovers kappaS reasonably well", {
  skip_on_cran()

  dat <- expand.grid(
    X = seq(0, 1, length.out = 8),
    Y = seq(0, 1, length.out = 8),
    year = 1:8
  )
  dat$x_s <- as.numeric(scale(
    sin(dat$X * 2 * pi) + cos(dat$Y * 2 * pi) + rnorm(nrow(dat), sd = 0.15)
  ))
  mesh <- make_dl_sim_mesh(dat)

  true_kappa <- 1.3
  sim <- simulate_new(
    formula = ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    range = 0.5,
    sigma_O = 0,
    phi = 0.02,
    B = c(0.3, 1.1),
    covariate_diffusion = ~ space(x_s),
    diffusion_kappaS = true_kappa,
    seed = 123
  )

  dat$y <- sim$observed
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    covariate_diffusion = ~ space(x_s),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  est_kappa <- fit$tmb_obj$report()$kappaS_dl[1]
  expect_equal(est_kappa, true_kappa, tolerance = 0.01)
})

test_that("simulated time diffusion recovers kappaT reasonably well", {
  skip_on_cran()

  dat <- expand.grid(
    X = seq(0, 1, length.out = 8),
    Y = seq(0, 1, length.out = 8),
    year = 1:10
  )
  dat$x_t <- as.numeric(scale(
    sin(dat$year / 2) + 0.2 * dat$year + rnorm(nrow(dat), sd = 0.2)
  ))
  mesh <- make_dl_sim_mesh(dat)

  true_kappa <- 0.45
  sim <- simulate_new(
    formula = ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    range = 0.5,
    sigma_O = 0,
    phi = 0.02,
    B = c(0.3, 1.1),
    covariate_diffusion = ~ time(x_t),
    diffusion_rhoT = true_kappa,
    seed = 124
  )

  dat$y <- sim$observed
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    covariate_diffusion = ~ time(x_t),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  est_rho <- fit$tmb_obj$report()$rhoT[1]
  expect_equal(est_rho, true_kappa, tolerance = 0.01)
})
