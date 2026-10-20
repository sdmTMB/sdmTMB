test_that("TMB IID simulation works", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:10, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.1,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  fit <- sdmTMB(observed ~ a1, sim_dat, mesh = mesh, time = "year")
  b <- tidy(fit)
  b
  expect_equal(b$estimate[b$term == "a1"], -0.4, tolerance = 0.1)
  expect_equal(b$estimate[b$term == "(Intercept)"], 0.2, tolerance = 0.2)
  b <- tidy(fit, "ran_pars")
  b
})

test_that("TMB AR1 simulation works", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:10, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.1,
    phi = 0.1,
    sigma_O = 0,
    seed = 42,
    rho = 0.8, #<
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  fit <- sdmTMB(observed ~ a1, sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "ar1", spatial = "off"
  )
  b <- tidy(fit)
  b
  b <- tidy(fit, "ran_pars")
  b
  rho_hat <- b$estimate[b$term == "rho"]
  expect_true(rho_hat > 0.7 && rho_hat < 0.9)
  sigma_E_hat <- b$estimate[b$term == "sigma_E"]
  expect_true(sigma_E_hat > 0.07 && sigma_E_hat < 0.13)
})

test_that("TMB RW simulation works", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(2000), Y = runif(2000),
    a1 = rnorm(2000), year = rep(1:20, each = 200)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.3,
    sigma_E = 0.1,
    phi = 0.05,
    sigma_O = 0,
    seed = 42,
    rho = 1, #<
    B = 0
  )
  fit_ar1 <- sdmTMB(observed ~ 0,
    sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "ar1", #<
    spatial = "off"
  )
  b_ar1 <- tidy(fit_ar1, "ran_pars")
  b_ar1
  rho_hat <- b_ar1$estimate[b_ar1$term == "rho"]
  expect_true(rho_hat > 0.9)
  sigma_E_hat_ar1 <- b_ar1$estimate[b_ar1$term == "sigma_E"]

  fit_rw <- sdmTMB(observed ~ 0,
    sim_dat,
    mesh = mesh, time = "year",
    spatiotemporal = "rw", #<
    spatial = "off"
  )
  b_rw <- tidy(fit_rw, "ran_pars")
  b_rw
  sigma_E_hat_rw <- b_rw$estimate[b_rw$term == "sigma_E"]
  sigma_E_hat_rw
  expect_true(sigma_E_hat_rw > 0.07 && sigma_E_hat_rw < 1.13)
  expect_true(sigma_E_hat_rw < sigma_E_hat_ar1)
})

test_that("sdmTMB_simulate can generate time-varying RW effects", {
  skip_on_cran()

  set.seed(99)
  predictor_dat <- data.frame(
    X = runif(600), Y = runif(600),
    year = rep(seq_len(60), each = 10)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.3)

  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.4,
    sigma_E = 0,
    phi = 0.1,
    sigma_O = 0,
    seed = 101,
    B = 0,
    time_varying = ~1,
    time_varying_type = "rw0",
    sigma_V = 0.25
  )

  fit <- sdmTMB(observed ~ 1,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off", spatiotemporal = "off",
    time_varying = ~1, time_varying_type = "rw0"
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_equal(exp(s$ln_tau_V)[1, 1], 0.25, tolerance = 0.1)
})

test_that("sdmTMB_simulate can generate time-varying AR1 effects", {
  skip_on_cran()

  set.seed(77)
  predictor_dat <- data.frame(
    X = runif(800), Y = runif(800),
    year = rep(seq_len(80), each = 10)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.3)

  rho_true <- 0.6
  sigma_v_true <- 0.2

  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.4,
    sigma_E = 0,
    phi = 0.1,
    sigma_O = 0,
    seed = 202,
    B = 0,
    time_varying = ~1,
    time_varying_type = "ar1",
    sigma_V = sigma_v_true,
    rho_time = rho_true
  )

  fit <- sdmTMB(observed ~ 1,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off", spatiotemporal = "off",
    time_varying = ~1, time_varying_type = "ar1"
  )
  s <- as.list(fit$sd_report, "Estimate")
  m121 <- function(x) 2 * plogis(x) - 1
  expect_equal(exp(s$ln_tau_V)[1, 1], sigma_v_true, tolerance = 0.1)
  expect_equal(m121(s$rho_time_unscaled)[1, 1], rho_true, tolerance = 0.15)
})

test_that("sdmTMB_simulate supports multiple AR1 time-varying coefficients", {
  skip_on_cran()

  set.seed(54)
  predictor_dat <- data.frame(
    X = runif(900), Y = runif(900),
    year = rep(seq_len(90), each = 10),
    cov1 = rnorm(900),
    cov2 = rnorm(900)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.3)

  sigma_true <- c(0.15, 0.3)
  rho_true <- c(0.4, -0.5)

  sim_dat <- sdmTMB_simulate(
    formula = ~1 + cov1 + cov2,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.4,
    sigma_E = 0,
    phi = 0.1,
    sigma_O = 0,
    seed = 909,
    B = c(0, 0.5, -0.25),
    time_varying = ~0 + cov1 + cov2,
    time_varying_type = "ar1",
    sigma_V = sigma_true,
    rho_time = rho_true
  )

  fit <- sdmTMB(observed ~ cov1 + cov2,
    data = sim_dat, mesh = mesh, time = "year",
    spatial = "off", spatiotemporal = "off",
    time_varying = ~0 + cov1 + cov2,
    time_varying_type = "ar1"
  )
  s <- as.list(fit$sd_report, "Estimate")
  m121 <- function(x) 2 * plogis(x) - 1

  expect_equal(exp(s$ln_tau_V)[, 1], sigma_true, tolerance = 0.1)
  expect_equal(m121(s$rho_time_unscaled)[, 1], rho_true, tolerance = 0.2)
})

test_that("TMB breakpt sims work", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    phi = 0.02,
    sigma_O = 0.001,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.3, 0)
  )
  expect_lt(max(sim_dat$observed), 0.1)
  sim_dat$a1 <- predictor_dat$a1
  fit <- sdmTMB(observed ~ 1 + breakpt(a1), sim_dat, mesh = mesh, spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1L))
  b <- tidy(fit)
  expect_gt(b$estimate[b$term == "a1-breakpt"], -0.02)
  expect_lt(b$estimate[b$term == "a1-breakpt"], 0.02)
  expect_gt(b$estimate[b$term == "a1-slope"], 0.28)
  expect_lt(b$estimate[b$term == "a1-slope"], 0.32)

  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + logistic(a1),
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    phi = 0.001,
    sigma_O = 0.001,
    seed = 42,
    B = 0,
    threshold_coefs = c(0.2, 0.4, 0.5)
  )
  expect_lt(max(sim_dat$observed), 0.53)
  expect_gt(min(sim_dat$observed), -0.05)
})

test_that("simulate.sdmTMB returns the right length", {
  skip_on_cran()
  pcod$os <- rep(log(0.01), nrow(pcod)) # offset
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0,
    time_varying = ~ 1,
    offset = pcod$os,
    family = tweedie(link = "log"),
    spatial = "off",
    time = "year",
    extra_time = c(2006, 2008, 2010, 2012, 2014, 2016),
    spatiotemporal = "off"
  )
  s <- simulate(m, nsim = 2)
  expect_equal(nrow(s), nrow(pcod))
})

test_that("simulate.sdmTMB works with sizes in binomial GLMs #465", {
  skip_on_cran()
  set.seed(1)
  w <- sample(1:10, size = 200, replace = TRUE)
  x <- rnorm(200)
  dat <- data.frame(y = stats::rbinom(length(w), size = w, prob = plogis(x * 0.5)))
  dat$prop <- dat$y / w
  dat$x <- x
  fit <- sdmTMB(prop ~ x, data = dat, weights = w, family = binomial(), spatial = "off")
  dat$X <- dat$Y <- NA # FIXME
  set.seed(1)
  snd <- simulate(fit, nsim = 500, newdata = dat, size = w)
  set.seed(1)
  s <- simulate(fit, nsim = 500)
  expect_equal(as.matrix(s), as.matrix(snd), ignore_attr = TRUE)
  expect_true(nrow(snd) == 200L)
  expect_true(ncol(snd) == 500L)
  sim_means <- rowMeans(snd)
  expect_gt(cor(sim_means, dat$y), 0.7)
})

test_that("coarse meshes with zeros in simulation still return fields #370", {
  set.seed(123)
  predictor_dat <- data.frame(
    X = runif(100), Y = runif(100),
    a1 = rnorm(100), year = rep(1:2, each = 50))
  mesh <- sdmTMB::make_mesh(predictor_dat, xy_cols = c("X", "Y"), n_knots = 30)
  sim_dat <- sdmTMB::sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.5,
    sigma_E = 0.1,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 42,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  nm <- names(sim_dat)
  expect_true("omega_s" %in% nm)
  expect_true("epsilon_st" %in% nm)
  expect_false("zeta_s" %in% nm)
})

test_that("simulate without observation error works for binomial likelihoods #431", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit.dg <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh, family = delta_gamma(type="standard")
  )
  s.dg <- simulate(
    fit.dg,
    newdata = qcs_grid_small,
    type = "mle-mvn", # fixed effects at MLE values and random effect MVN draws
    mle_mvn_samples = "multiple", # take an MVN draw for each sample
    nsim = 50, # increase this for more stable results
    observation_error = FALSE, # do not include observation error
    seed = 23859
  )
  expect_gt(min(s.dg), 0)
  m <- apply(s.dg, 1, mean)
  p <- predict(fit.dg, newdata = qcs_grid_small)
  expect_gt(cor(plogis(p$est1) * exp(p$est2), m), 0.98)

  fit.b <- sdmTMB(present ~ 1,
    data = pcod, mesh = mesh, family = binomial()
  )
  s.b <- simulate(
    fit.b,
    newdata = qcs_grid_small,
    type = "mle-mvn", # fixed effects at MLE values and random effect MVN draws
    mle_mvn_samples = "multiple", # take an MVN draw for each sample
    nsim = 50, # increase this for more stable results
    observation_error = FALSE, # do not include observation error
    seed = 23859
  )
  expect_gt(min(s.dg), 0)
  m <- apply(s.b, 1, mean)
  p <- predict(fit.b, newdata = qcs_grid_small)
  expect_gt(cor(plogis(p$est), m), 0.95)

  # with size specified (but wrong length at first)
  expect_error({simulate(
    fit.b,
    newdata = qcs_grid_small,
    nsim = 1,
    observation_error = FALSE,
    size = c(1, 2, 3)
  )}, regexp = "size")

  set.seed(1)
  w <- sample(1:9, size = nrow(qcs_grid_small), replace = TRUE)
  s.b1 <- simulate(
    fit.b,
    newdata = qcs_grid_small,
    type = "mle-mvn",
    mle_mvn_samples = "multiple",
    nsim = 50,
    observation_error = FALSE,
    seed = 23859,
    size = w
  )
  expect_true(max(s.b1) > 1)
  expect_equal(mean(s.b1[1,]), m[1] * w[1])
  expect_equal(mean(s.b1[51,]), m[51] * w[51])
})

test_that("simulate without observation error works for binomial likelihoods and Poisson-link delta", {
  skip_on_cran()
  skip_on_ci()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit.dg <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh, family = delta_gamma(type="poisson-link")
  )
  s.dg <- simulate(
    fit.dg,
    newdata = qcs_grid_small,
    type = "mle-mvn", # fixed effects at MLE values and random effect MVN draws
    mle_mvn_samples = "multiple", # take an MVN draw for each sample
    nsim = 200, # increase this for more stable results
    observation_error = FALSE, # do not include observation error
    seed = 23859
  )
  expect_gt(min(s.dg), 0)
  m <- apply(s.dg, 1, mean)
  p <- predict(fit.dg, newdata = qcs_grid_small)
  expect_gt(cor(exp(p$est1) * exp(p$est2), m), 0.98)
})

test_that("simulate_re redraws only the named random effects", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  for (backend in c("rtmb", "tmb")) {
    fit <- sdmTMB(density ~ 0 + as.factor(year), data = pcod, mesh = mesh,
      family = tweedie(), time = "year", spatiotemporal = "iid",
      control = sdmTMBcontrol(backend = backend))
    sim <- function(...) {
      simulate(fit, return_tmb_report = TRUE, seed = 1, silent = TRUE, ...)
    }

    # mle-eb: held omega equals the EB estimates
    cond <- sim(nsim = 1)[[1]]
    r <- sim(nsim = 2, simulate_re = "spatiotemporal")
    expect_equal(r[[1]]$omega_s, cond$omega_s)
    expect_identical(r[[1]]$omega_s, r[[2]]$omega_s)
    expect_false(identical(r[[1]]$epsilon_st, r[[2]]$epsilon_st))

    # mle-mvn single: one shared posterior draw, different from EB
    r <- sim(nsim = 3, type = "mle-mvn", mle_mvn_samples = "single",
      simulate_re = "spatiotemporal")
    expect_identical(r[[1]]$omega_s, r[[3]]$omega_s)
    expect_false(isTRUE(all.equal(r[[1]]$omega_s, cond$omega_s)))
    expect_false(identical(r[[1]]$epsilon_st, r[[2]]$epsilon_st))

    # mle-mvn multiple: held omega differs across nsim
    r <- sim(nsim = 2, type = "mle-mvn", mle_mvn_samples = "multiple",
      simulate_re = "spatiotemporal")
    expect_false(identical(r[[1]]$omega_s, r[[2]]$omega_s))

    # naming both redraws both
    r <- sim(nsim = 2, simulate_re = c("spatial", "spatiotemporal"))
    expect_false(identical(r[[1]]$omega_s, r[[2]]$omega_s))

    # newdata
    nd <- replicate_df(qcs_grid_small, "year", unique(pcod$year))
    r <- sim(nsim = 2, newdata = nd, simulate_re = "spatiotemporal")
    expect_equal(r[[1]]$omega_s, cond$omega_s)
    expect_identical(r[[1]]$omega_s, r[[2]]$omega_s)
    expect_false(identical(r[[1]]$epsilon_st, r[[2]]$epsilon_st))

    # matrix output and reproducibility
    y <- simulate(fit, nsim = 2, simulate_re = "spatiotemporal", seed = 1, silent = TRUE)
    expect_equal(dim(y), c(nrow(pcod), 2L))
    expect_identical(y, simulate(fit, nsim = 2, simulate_re = "spatiotemporal",
      seed = 1, silent = TRUE))

    # all named is the same as re_form = NA
    expect_identical(
      simulate(fit, nsim = 2, seed = 1, silent = TRUE, re_form = NA),
      simulate(fit, nsim = 2, seed = 1, silent = TRUE,
        simulate_re = c("spatial", "spatiotemporal", "spatial_varying",
          "group_re", "time_varying"))
    )
  }

  expect_error(simulate(fit, simulate_re = "spatiotemporal", re_form = ~0), "only one")
  expect_error(simulate(fit, simulate_re = "spatio"), "Unknown")
  expect_error(simulate(fit, simulate_re = "omega_s"), "Unknown")
})

test_that("simulate_re works for delta models and time-varying effects", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit <- sdmTMB(density ~ 1, data = pcod, mesh = mesh,
    family = delta_gamma(type = "poisson-link"), time = "year",
    spatiotemporal = "iid")
  cond <- simulate(fit, nsim = 1, return_tmb_report = TRUE, seed = 1, silent = TRUE)[[1]]
  r <- simulate(fit, nsim = 2, simulate_re = "spatiotemporal",
    return_tmb_report = TRUE, seed = 1, silent = TRUE)
  expect_equal(r[[1]]$omega_s, cond$omega_s)
  expect_equal(ncol(r[[1]]$omega_s), 2L)
  expect_false(identical(r[[1]]$epsilon_st, r[[2]]$epsilon_st))

  fit <- sdmTMB(density ~ 0, time_varying = ~1, data = pcod, mesh = mesh,
    family = tweedie(), time = "year", spatiotemporal = "iid")
  r <- simulate(fit, nsim = 2, simulate_re = "time_varying",
    return_tmb_report = TRUE, seed = 1, silent = TRUE)
  expect_identical(r[[1]]$omega_s, r[[2]]$omega_s)
  expect_identical(r[[1]]$epsilon_st, r[[2]]$epsilon_st)
  expect_false(identical(r[[1]]$b_rw_t, r[[2]]$b_rw_t))
})

test_that("simulate_re = 'group_re' simulates correlated random slopes", {
  skip_on_cran()
  set.seed(1)
  ng <- 100
  g <- rep(seq_len(ng), each = 10)
  b <- matrix(rnorm(2 * ng, 0, c(0.5, 0.3)), ncol = 2, byrow = TRUE)
  d <- data.frame(g = factor(g), x = rnorm(length(g)))
  d$y <- 1 + b[g, 1] + (0.5 + b[g, 2]) * d$x + rnorm(length(g), 0, 0.3)
  for (backend in c("rtmb", "tmb")) {
    fit <- sdmTMB(y ~ x + (1 + x | g), data = d, spatial = "off",
      control = sdmTMBcontrol(backend = backend))
    r <- simulate(fit, nsim = 200, simulate_re = "group_re",
      return_tmb_report = TRUE, seed = 1, silent = TRUE)
    b_sim <- do.call(rbind, lapply(r, function(.x) {
      matrix(.x$re_b_pars[, 1], ncol = 2, byrow = TRUE)
    }))
    sds <- tidy(fit, "ran_pars")
    sds <- sds$estimate[grepl("^sd__", sds$term)]
    expect_equal(apply(b_sim, 2, sd), sds, tolerance = 0.05)
    expect_false(identical(r[[1]]$re_b_pars, r[[2]]$re_b_pars))
  }
})
