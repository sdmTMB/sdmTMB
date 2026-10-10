test_that("Basic prior parsing works", {
  expect_equal(normal(0, 1), matrix(c(0, 1), ncol = 2L), ignore_attr = TRUE)
  expect_equal(halfnormal(0, 1), matrix(c(0, 1), ncol = 2L), ignore_attr = TRUE)
  expect_equal(pc_matern(5, 5, 0.05, 0.05), c(5, 5, 0.05, 0.05), ignore_attr = TRUE)

  expect_error(normal(NA, 1))
  expect_error(normal(0, -1))
  expect_error(normal(1, NA))
  expect_error(normal(c(1, 1), 1))

  expect_error(pc_matern(NA, 1))
  expect_error(pc_matern(1, NA))
  expect_error(pc_matern(-1, 1))
  expect_error(pc_matern(1, -1))
  expect_error(pc_matern(1, 1, 1.5, 0.05))
  expect_error(pc_matern(1, 1, -1.5, 0.05))
  expect_error(pc_matern(1, 1, 0.05, -0.05))
  expect_error(pc_matern(1, 1, 0.05, NA))
})

test_that("Fixed-effect normal priors use SDs and expand scalars", {
  get_b_prior <- function(b) {
    fit <- sdmTMB(y ~ x, data = data.frame(y = c(1, 2, 3), x = c(-1, 0, 1)),
      spatial = "off", do_fit = FALSE, priors = sdmTMBpriors(b = b))
    fit$tmb_data[c("priors_b_n", "priors_b_mean", "priors_b_Sigma")]
  }
  p <- get_b_prior(normal(c(0, 1), c(2, 3)))
  expect_equal(p$priors_b_Sigma, diag(c(4, 9)))
  expect_message(p <- get_b_prior(normal(1, 2)), "Expanding")
  expect_equal(p$priors_b_n, 2L)
  expect_equal(p$priors_b_mean, c(1, 1))
  expect_equal(p$priors_b_Sigma, diag(c(4, 4)))
  expect_message(p <- get_b_prior(mvnormal(1, matrix(4))), "Expanding")
  expect_equal(p$priors_b_Sigma, diag(c(4, 4)))
  expect_equal(get_b_prior(normal(NA, NA))$priors_b_n, 0L)
})

test_that("PC prior Jacobian is log(range * sigma)", {
  prior <- pc_matern(5, 1)
  log_tau <- 0.7
  log_kappa <- 0.3
  range <- sqrt(8) / exp(log_kappa)
  sigma <- exp(-log_tau - log_kappa) / sqrt(4 * pi)
  for (share_range in c(FALSE, TRUE)) {
    jacobian <- rtmb_pc_matern(log_tau, log_kappa, prior,
      include_range = !share_range, stan = TRUE) -
      rtmb_pc_matern(log_tau, log_kappa, prior,
        include_range = !share_range, stan = FALSE)
    expect_equal(jacobian, log(sigma) + if (share_range) 0 else log(range))
  }
})

test_that("Prior fitting works", {
  skip_on_cran()
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011

  # no priors
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
    spatiotemporal = "AR1",
    share_range = FALSE
  )

  # population-effects missing a prior; should error
  expect_error(
    {
      mp <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
        data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
        share_range = FALSE, spatiotemporal = "AR1",
        priors = sdmTMBpriors(
          b = normal(c(0, 0, NA, NA, NA), c(2, 2, NA, NA, NA))
        )
      )
    },
    regexp = "prior"
  )

  # all the priors
  mp <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2,
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
    share_range = FALSE, spatiotemporal = "AR1",
    priors = sdmTMBpriors(
      # b = normal(c(0, 0, NA, NA, NA, NA), c(2, 2, NA, NA, NA, NA)),
      phi = halfnormal(0, 10),
      # tweedie_p = normal(1.5, 2),
      ar1_rho = normal(0, 1),
      matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
      matern_st = pc_matern(range_gt = 5, sigma_lt = 1)
    )
  )

  expect_lt(abs(mp$model$par[1]), abs(m$model$par[1]))
  expect_gt(mp$model$par[["ln_tau_O"]], m$model$par[["ln_tau_O"]])
  expect_gt(mp$model$par[["ln_tau_E"]], m$model$par[["ln_tau_E"]])
})

# FIXME: random-slopes: get these priors working again?
# test_that("Priors on random intercept SDs work", {
#   skip_on_ci()
#   skip_on_cran()
#
#   pcod$fyear <- as.factor(pcod$year)
#   fit0 <- sdmTMB(
#     density ~ 1 + (1 | fyear), family = tweedie(),
#     data = pcod, spatial = "off"
#   )
#   fit1 <- sdmTMB(
#     density ~ 1 + (1 | fyear), family = tweedie(),
#     data = pcod, spatial = "off",
#     priors = sdmTMBpriors(sigma_G = halfnormal(0, 0.1))
#   )
#   t0 <- tidy(fit0, "ran_pars")
#   t1 <- tidy(fit1, "ran_pars")
#   G0 <- t0$estimate[t0$term == "sigma_G"]
#   G1 <- t1$estimate[t0$term == "sigma_G"]
#   expect_lt(G1, G0) # prior reduces SD
#
#   fit0 <- sdmTMB(
#     density ~ 1 + (1 | fyear), family = poisson(),
#     data = pcod, spatial = "off"
#   )
#   pcod$fake_count <- round(pcod$density)
#   pcod$obs_id <- as.factor(seq_len(nrow(pcod)))
#   fit0 <- sdmTMB(
#     fake_count ~ 1 + (1 | fyear) + (1 | obs_id), family = poisson(),
#     data = pcod, spatial = "off"
#   )
#   fit1 <- sdmTMB(
#     fake_count ~ 1 + (1 | fyear) + (1 | obs_id), family = poisson(),
#     data = pcod, spatial = "off",
#     priors = sdmTMBpriors(sigma_G = halfnormal(c(0, 0), c(0.1, 0.1)))
#   )
#   fit2 <- sdmTMB(
#     fake_count ~ 1 + (1 | fyear) + (1 | obs_id), family = poisson(),
#     data = pcod, spatial = "off",
#     priors = sdmTMBpriors(sigma_G = halfnormal(c(1, 0), c(0.1, 0.5)))
#   )
#   t0 <- tidy(fit0, "ran_pars")
#   t1 <- tidy(fit1, "ran_pars")
#   t2 <- tidy(fit2, "ran_pars")
#   G0 <- t0$estimate[t0$term == "sigma_G"][2]
#   G1 <- t1$estimate[t0$term == "sigma_G"][2]
#   G2 <- t2$estimate[t0$term == "sigma_G"][2]
#   expect_lt(G1, G0) # prior reduces SD
#   expect_lt(G1, G2) # stronger prior reduces SD more
#
#   # high mean prior keeps SD away from zero:
#   G1 <- t1$estimate[t0$term == "sigma_G"][1]
#   G2 <- t2$estimate[t0$term == "sigma_G"][1]
#   expect_gt(G2, G1)
# })

test_that("Additional priors work", {
  skip_on_cran()
  d <- pcod_2011
  pcod_spde <- pcod_mesh_2011

  # priors with one covariate/intercept only and the normal scale not 1:
  m_norm <- sdmTMB(density ~ 1,
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(0, 10)),
    spatial = "off", spatiotemporal = "off"
  )

  # univariate normal priors
  m_norm <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(rep(0, 6), rep(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_norm), "sdmTMB")
  m_mvn <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = mvnormal(rep(0, 6), diag(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_mvn), "sdmTMB")
  m_mvn_na <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = mvnormal(c(NA, 0, 0, 0, 0, 0), diag(1, 6))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_mvn_na), "sdmTMB")

  m_norm_na <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, mesh = pcod_spde, family = tweedie(link = "log"),
    priors = sdmTMBpriors(b = normal(c(NA, rep(0, 5)), c(NA, rep(1, 5)))),
    spatial = "off", spatiotemporal = "off"
  )
  expect_identical(class(m_norm_na), "sdmTMB")

  expect_equal(tidy(m_norm), tidy(m_mvn), tolerance = 0.0001)
  expect_equal(tidy(m_norm_na), tidy(m_mvn_na), tolerance = 0.0001)
  expect_true(
    abs(tidy(m_mvn)$estimate[tidy(m_mvn)$term == "depth_scaled"]) <
      abs(tidy(m_mvn_na)$estimate[tidy(m_mvn_na)$term == "depth_scaled"])
  )
})

test_that("Threshold priors work", {
  skip_on_cran()
  d <- pcod_2011

  m <- sdmTMB(density ~ 0 + logistic(depth_scaled),
    data = d,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(threshold_logistic_s50 = normal(-1, 0.005)),
    spatial = "off", spatiotemporal = "off"
  )
  x <- tidy(m)
  expect_equal(x$estimate[x$term == "depth_scaled-s50"], -1, tolerance = 0.01)

  m <- sdmTMB(density ~ 0 + logistic(depth_scaled),
    data = d,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(threshold_logistic_s95 = normal(-1, 0.005)),
    spatial = "off", spatiotemporal = "off"
  )
  x <- tidy(m)
  expect_equal(x$estimate[x$term == "depth_scaled-s95"], -1, tolerance = 0.01)

  m <- sdmTMB(density ~ 0 + logistic(depth_scaled),
    data = d,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(threshold_logistic_smax = normal(4.0, 0.005)),
    spatial = "off", spatiotemporal = "off"
  )
  x <- tidy(m)
  expect_equal(x$estimate[x$term == "depth_scaled-smax"], 4.0, tolerance = 0.01)

  m <- sdmTMB(density ~ 0 + breakpt(depth_scaled),
    data = d,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(threshold_breakpt_slope = normal(-4, 0.005)),
    spatial = "off", spatiotemporal = "off"
  )
  x <- tidy(m)
  expect_equal(x$estimate[x$term == "depth_scaled-slope"], -4, tolerance = 0.01)

  m <- sdmTMB(density ~ 0 + breakpt(depth_scaled),
    data = d,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(threshold_breakpt_cut = normal(-1, 0.005)),
    spatial = "off", spatiotemporal = "off"
  )
  x <- tidy(m)
  expect_equal(x$estimate[x$term == "depth_scaled-breakpt"], -1, tolerance = 0.01)
})

test_that("the Stan Jacobian gives a proper PC Matern prior on (log_tau, log_kappa)", {
  # With the correct Jacobian, the implied density integrates to 1
  prior <- pc_matern(range_gt = 5, sigma_lt = 2)
  h <- 0.025
  lt <- seq(-15, 15, by = h)
  lk <- seq(-15, 10, by = h)
  density <- outer(lt, lk, function(lt, lk) {
    exp(rtmb_pc_matern(lt, lk, prior, stan = TRUE))
  })
  expect_equal(sum(density) * h^2, 1, tolerance = 1e-4)
})

test_that("Priors apply the right terms and Jacobians in both backends", {
  make_obj <- function(backend, priors, bayesian = FALSE, ...) {
    sdmTMB(..., priors = priors, bayesian = bayesian, do_fit = FALSE,
      control = sdmTMBcontrol(backend = backend, multiphase = FALSE))$tmb_obj
  }
  for (backend in c("tmb", "rtmb")) {
    # A spatiotemporal PC prior keeps its range term when no spatial PC prior
    # supplies the shared range
    objs <- lapply(list(sdmTMBpriors(), sdmTMBpriors(matern_st = pc_matern(5, 1))),
      make_obj, backend = backend, formula = density ~ 1, data = pcod_2011,
      mesh = pcod_mesh_2011, time = "year", family = tweedie())
    p <- objs[[1]]$par
    p[names(p) == "ln_tau_E"] <- 1
    p[names(p) == "ln_kappa"] <- -1
    range <- sqrt(8) / exp(-1)
    sigma <- 1 / sqrt(4 * pi)
    l_range <- -log(0.05) * 5
    l_sigma <- -log(0.05)
    expect_equal(as.numeric(objs[[2]]$fn(p) - objs[[1]]$fn(p)),
      -(log(l_range) - 2 * log(range) - l_range / range +
        log(l_sigma) - l_sigma * sigma), label = backend)

    # Logistic-threshold s95 Jacobian
    objs <- lapply(c(FALSE, TRUE), make_obj, backend = backend,
      priors = sdmTMBpriors(threshold_logistic_s50 = normal(0, 1),
        threshold_logistic_s95 = normal(1, 1)),
      formula = y ~ 1 + logistic(x), spatial = "off",
      data = data.frame(x = seq(-2, 2, length.out = 10),
        y = seq(-1, 1, length.out = 10)))
    p <- objs[[1]]$par
    p[names(p) == "b_threshold"] <- c(0, log(2), 1)
    expect_equal(as.numeric(objs[[2]]$fn(p) - objs[[1]]$fn(p)), -log(2),
      label = backend)
  }
})

test_that("Lognormal sigma_V priors work", {
  expect_equal(lognormal_prior(log(0.2), 0.5), matrix(c(log(0.2), 0.5), ncol = 2L),
    ignore_attr = TRUE)
  expect_error(lognormal_prior(0, -1))
  expect_error(lognormal_prior(NA, 1))
  expect_error(sdmTMBpriors(sigma_V = normal(0, 1)))

  make_obj <- function(priors, bayesian = FALSE) {
    sdmTMB(density ~ 0, time = "year", time_varying = ~ 1 + depth_scaled,
      data = pcod_2011, spatial = "off", spatiotemporal = "off",
      family = tweedie(), priors = priors, bayesian = bayesian, do_fit = FALSE,
      control = sdmTMBcontrol(multiphase = FALSE))$tmb_obj
  }
  sigma <- c(0.3, 0.6)
  for (bayesian in c(FALSE, TRUE)) {
    # vector of priors with NA = no prior on the second SD
    objs <- lapply(list(sdmTMBpriors(),
      sdmTMBpriors(sigma_V = lognormal_prior(c(log(0.2), NA), c(0.5, NA)))),
      make_obj, bayesian = bayesian)
    p <- objs[[1]]$par
    p[names(p) == "ln_tau_V"] <- log(sigma)
    expected <- -dlnorm(sigma[1], log(0.2), 0.5, log = TRUE)
    if (bayesian) expected <- expected - log(sigma[1])
    expect_equal(as.numeric(objs[[2]]$fn(p) - objs[[1]]$fn(p)), expected,
      tolerance = 1e-6, label = paste("bayesian =", bayesian))
  }
  # gamma priors are unchanged
  objs <- lapply(list(sdmTMBpriors(), sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5))),
    make_obj)
  p <- objs[[1]]$par
  p[names(p) == "ln_tau_V"] <- log(sigma)
  expect_equal(as.numeric(objs[[2]]$fn(p) - objs[[1]]$fn(p)),
    -sum(dgamma(sigma, shape = 1 / 0.5^2, scale = 0.5^2 * 0.2, log = TRUE)),
    tolerance = 1e-6)
})
