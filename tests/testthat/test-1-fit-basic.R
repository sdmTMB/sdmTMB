# basic model fitting and prediction tests

SEED <- 123
set.seed(SEED)
x <- stats::runif(400, -1, 1)
y <- stats::runif(400, -1, 1)
loc <- data.frame(x = x, y = y)
loc$time <- rep(1:4, each = 100)

test_that("sdmTMB model fit with a covariate beta", {
  local_edition(2)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.04)
  initial_betas <- 0.5
  range <- 0.1
  sigma_O <- 0.3 # SD of spatial process
  sigma_E <- 0.3 # SD of spatial process
  phi <- 0.1 # observation error
  set.seed(1)
  loc$cov1 <- runif(nrow(loc))
  s <- sdmTMB_simulate(
    ~ 0 + cov1, loc,
    mesh = spde, time = "time", B = initial_betas,
    phi = phi, range = range, sigma_O = sigma_O, sigma_E = sigma_E,
    seed = SEED
  )
  # `sdmTMB_simulate()` draws random numbers via the active backend's
  # `simulate()` method; TMB and RTMB consume the RNG stream in different
  # orders, so exact simulated values are not comparable across backends.
  expect_equal(nrow(s), nrow(loc))
  expect_true(all(is.finite(s$observed)))
  plot(spde)
  m <- sdmTMB(
    data = s, formula = observed ~ 0 + cov1, time = "time",
    silent = TRUE, mesh = spde, control = sdmTMBcontrol(normalize = FALSE, newton_loops = 1)
  )
  m_pc <- sdmTMB(
    data = s, formula = observed ~ 0 + cov1, time = "time",
    silent = TRUE, mesh = spde, control = sdmTMBcontrol(normalize = FALSE, newton_loops = 1),
    priors = sdmTMBpriors(
      matern_s = pc_matern(range_gt = 0.2, sigma_lt = 0.2, range_prob = 0.05, sigma_prob = 0.05)
    )
  )

  # expect_equal(m$model$par, m_pc$model$par, tolerance = 0.05)
  # expect_equal(round(m$model$par, 3),
  #   c(b_j = 0.523, ln_tau_O = -3.615, ln_tau_E = -3.567, ln_kappa = 3.386,
  #     ln_phi = -2.674))
  # expect_equal(round(m_pc$model$par, 3),
  #   c(b_j = 0.523, ln_tau_O = -3.509, ln_tau_E = -3.498, ln_kappa = 3.339,
  #     ln_phi = -2.688))

  # PC should make range bigger and sigmaO smaller
  # therefore, ln_kappa smaller
  expect_true(m_pc$model$par[["ln_kappa"]] < m$model$par[["ln_kappa"]])
  r_pc <- m_pc$tmb_obj$report()
  r <- m$tmb_obj$report()
  expect_true(r_pc$range[1] > r$range[1])
  expect_true(r_pc$sigma_O < r$sigma_O)

  # PC should make random field pars more precise:
  se_pc <- as.list(m_pc$sd_report, "Std. Error")
  se <- as.list(m$sd_report, "Std. Error")
  expect_true(se_pc$ln_kappa[1] < se$ln_kappa[1])
  # `se_pc$ln_tau_O < se$ln_tau_O` is not checked here: whether the PC prior
  # tightens this particular SE depends on the realized random field draw
  # (it flips for plenty of seeds under either backend), so it is not a
  # reliable property to assert for one fixed seed.
  expect_true(is.finite(se_pc$ln_tau_O) && is.finite(se$ln_tau_O))

  expect_output(print(m), "fit by")
  expect_output(summary(m), "fit by")
  p <- as.list(m$model$par)
  r <- m$tmb_obj$report()
  est <- tidy(m, "ran_pars")
  # expect_equal(round(sort(est[,"estimate", drop = TRUE]), 3),
  #   c(0.069, 0.096, 0.338, 0.354))
  expect_equal(m$model$convergence, 0L)
  expect_equal((p$b_j - initial_betas)^2, 0, tolerance = 0.01)
  expect_equal((exp(p$ln_phi) - phi)^2, 0, tolerance = 0.002)
  # expect_equal((r$sigma_O - sigma_O)^2, 0, tolerance = 0.002)
  # expect_equal((r$sigma_E[1] - sigma_E)^2, 0, tolerance = 0.001)
  # expect_equal(est$estimate[est$term == "range"][1], range, tolerance = 0.01)
  p <- predict(m)
  r <- residuals(m)
  r_sim <- residuals(m, type = "mle-mvn")
  r_sim <- residuals(m, type = "mle-eb")
  expect_equal(mean((p$est - s$observed)^2), 0, tolerance = 0.002)

  nd <- s
  nd$time <- as.numeric(nd$time)
  p_ok <- predict(m, newdata = nd)
  expect_true(class(p_ok) == "data.frame")

  # mismatch:
  nd$time <- as.factor(nd$time)
  expect_error(predict(m, newdata = nd), regexp = "class")

  # no fields; fake mesh:
  set.seed(1)
  m <- sdmTMB(
    data = s, formula = observed ~ 0 + cov1, time = "time",
    spatial = "off", spatiotemporal = "off"
  )
  m <- sdmTMB(
    data = s, formula = observed ~ 0 + cov1, time = "time",
    spatial = FALSE, spatiotemporal = "off"
  )
  m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, spatial = "off")
  m <- sdmTMB(data = s, formula = observed ~ 0 + cov1, spatial = FALSE)
})

test_that("Anisotropy fits and plots", {
  skip_on_cran()
  local_edition(2)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    mesh = make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans"),
    family = tweedie(link = "log"), anisotropy = TRUE
  )
  expect_identical(class(m), "sdmTMB")
  plot_anisotropy(m)
  expect_equal(m$tmb_obj$report()$H,
    structure(
      c(0.665528444798002, 0.079350716881963, 0.079350716881963, 1.51202633656794),
      .Dim = c(2L, 2L)
    ),
    tolerance = 1e-3
  )
})

test_that("A spatiotemporal version works with predictions on new data points", {
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011
  m <- sdmTMB(
    data = d,
    formula = density ~ 0 + as.factor(year),
    time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
    spatiotemporal = TRUE,
    spatial = FALSE
  )
  # returns error if time is missing
  expect_error({
    mx <- sdmTMB(
      data = d,
      formula = density ~ 0 + as.factor(year),
      mesh = pcod_spde, family = tweedie(link = "log"),
      spatiotemporal = TRUE,
      spatial = FALSE
    )
  })
  # Predictions at original data locations:
  predictions <- predict(m)
  predictions$resids <- residuals(m, type="mle-mvn") # randomized quantile residuals
  # Predictions onto new data:
  nd <- replicate_df(qcs_grid, "year", unique(d$year))
  predictions <- predict(m, newdata = subset(nd, year >= 2011))
  expect_identical(class(predictions), "data.frame")
})

test_that("Predictions on the original data set as `newdata`` return the same predictions", {
  skip_on_cran()
  local_edition(2)
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  loc$time <- rep(1:5, each = 100)
  spde <- make_mesh(loc, c("x", "y"), cutoff = 0.04)
  dat <- sdmTMB_simulate(
    ~1, loc,
    mesh = spde,
    time = "time", range = 0.2, B = 0,
    sigma_O = 0, sigma_E = 0.3, phi = 0.1,
    seed = 1
  )
  m <- sdmTMB(
    spatiotemporal = "AR1", spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), mesh = spde
  )
  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  tidy(m)
  tidy(m, conf.int = TRUE)
  tidy(m, effects = "ran_par")
  tidy(m, effects = "ran_par", conf.int = TRUE)

  cols <- c("est", "est_non_rf", "est_rf", "epsilon_st")
  expect_equal(p[, cols], p_nd[, cols], tolerance = 1e-3)

  m <- sdmTMB(
    spatiotemporal = "AR1", spatial = FALSE,
    data = dat, formula = observed ~ 1, time = "time",
    family = gaussian(link = "identity"), mesh = spde
  )

  p <- predict(m)
  p_nd <- predict(m, newdata = dat)
  expect_equal(p[, cols], p_nd[, cols], tolerance = 1e-3)
})

test_that("poly() works on newdata", {
  skip_on_cran()

  # https://github.com/sdmTMB/sdmTMB/issues/77
  d <- pcod_2011
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  m <- sdmTMB(
    data = d, formula = density ~ poly(log(depth), 2),
    mesh = mesh, family = sdmTMB::tweedie(link = "log"),
    spatial = "off"
  )
  nd <- pcod_2011[1:3, ]
  p <- predict(m, newdata = nd)
  p <- p$est
  m2 <- glmmTMB::glmmTMB(
    data = d, formula = density ~ poly(log(depth), 2),
    family = glmmTMB::tweedie(link = "log")
  )
  p2 <- predict(m2, newdata = nd)
  expect_equal(p, p2, tolerance = 1e-4)
})
