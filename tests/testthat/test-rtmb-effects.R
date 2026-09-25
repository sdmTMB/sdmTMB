test_that("RTMB spatially varying coefficients match TMB", {
  set.seed(46)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  d$response <- 0.3 + 0.4 * d$z + rnorm(n, sd = 0.3)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 18L, type = "kmeans")
  for (spatial in c("on", "off")) {
    fit <- sdmTMB(response ~ z, data = d, mesh = mesh,
      spatial = spatial, spatiotemporal = "off",
      spatial_varying = ~1 + z, family = gaussian(),
      control = sdmTMBcontrol(multiphase = FALSE), do_fit = FALSE)
    p <- fit$tmb_params
    p$b_j[] <- c(0.3, 0.4)
    p$ln_tau_Z[] <- 0.2
    p$ln_kappa[] <- 0.1
    p$zeta_s[] <- seq(-0.1, 0.1, length.out = length(p$zeta_s))
    cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random)
    rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")
    expect_equal(rt$fn(rt$par), cpp$fn(cpp$par), tolerance = 1e-7,
      info = spatial)
    expect_equal(rt$gr(rt$par), cpp$gr(cpp$par),
      ignore_attr = TRUE, tolerance = 1e-6, info = spatial)
    expected <- cpp$report()
    actual <- rt$report()
    for (name in c("sigma_O", "sigma_Z", "zeta_s_A", "eta_i")) {
      expect_equal(actual[[name]], expected[[name]], tolerance = 1e-6,
        info = paste(spatial, name))
    }
    if (spatial == "on") {
      expect_equal(
        sort(rownames(suppressWarnings(summary(sdreport_sdmTMB(rt), "report")))),
        sort(rownames(suppressWarnings(summary(sdreport_sdmTMB(cpp), "report")))))
    }
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    expect_equal(ncol(projection$proj_z_i), ncol(fit$tmb_data$z_i))
    proj_cpp <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random)$report()
    proj_rt <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")$report()
    expect_equal(proj_rt$proj_eta, proj_cpp$proj_eta, tolerance = 1e-6)
    expect_equal(proj_rt$proj_zeta_s_A, proj_cpp$proj_zeta_s_A,
      tolerance = 1e-6)
    expect_setequal(names(proj_rt), names(proj_cpp))
    simulated_data <- fit$tmb_data
    simulated_data$sim_re[[3L]] <- 1L
    sim <- make_sdmTMB_adfun(simulated_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")$simulate()
    expect_equal(dim(sim$zeta_s), dim(p$zeta_s))
    expect_false(isTRUE(all.equal(sim$zeta_s, p$zeta_s)))
    expect_equal(as.numeric(sim$eta_i),
      as.numeric(simulated_data$X_ij[[1L]] %*% p$b_j +
        simulated_data$offset_i +
        rowSums(sim$zeta_s_A[, , 1L] * simulated_data$z_i)),
      tolerance = 1e-6)
  }
})

test_that("RTMB time-varying coefficients match TMB", {
  set.seed(50)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    w = rnorm(n), time = rep(1:3, each = 20L))
  d$response <- rnorm(n)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 18L, type = "kmeans")
  for (time_type in c("rw", "rw0", "ar1")) {
    fit <- sdmTMB(response ~ z, data = d, mesh = mesh, time = "time",
      spatial = "off", spatiotemporal = "off",
      time_varying = ~0 + z + w, time_varying_type = time_type,
      family = gaussian(), control = sdmTMBcontrol(multiphase = FALSE),
      do_fit = FALSE)
    p <- fit$tmb_params
    p$b_rw_t[] <- seq(-0.2, 0.2, length.out = length(p$b_rw_t))
    p$ln_tau_V[] <- log(0.6)
    p$rho_time_unscaled[] <- 0.4
    cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random)
    rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")
    expect_equal(rt$fn(rt$par), cpp$fn(cpp$par), tolerance = 1e-7,
      info = time_type)
    expect_equal(rt$gr(rt$par), cpp$gr(cpp$par),
      ignore_attr = TRUE, tolerance = 1e-6, info = time_type)
    for (name in c("sigma_V", "eta_rw_i", "eta_i")) {
      expect_equal(rt$report()[[name]], cpp$report()[[name]],
        tolerance = 1e-6, info = paste(time_type, name))
    }
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    proj_cpp <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random)$report()
    proj_rt <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")$report()
    expect_equal(proj_rt$proj_eta, proj_cpp$proj_eta, tolerance = 1e-6)
    expect_equal(proj_rt$proj_rw_i, proj_cpp$proj_rw_i, tolerance = 1e-6)
    simulated_data <- fit$tmb_data
    simulated_data$sim_re[[5L]] <- 1L
    simulated_data$simulate_t <- c(0L, 1L, 1L)
    sim <- make_sdmTMB_adfun(simulated_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")$simulate()
    expect_equal(sim$b_rw_t[1L, , 1L], p$b_rw_t[1L, , 1L])
    expect_false(isTRUE(all.equal(sim$b_rw_t[2L, , 1L],
      p$b_rw_t[2L, , 1L])))
  }
})

test_that("RTMB correlated IID effects preserve TMB parameter order", {
  set.seed(51)
  n <- 40L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    w = rnorm(n), v = rnorm(n), g = factor(rep(1:5, each = 8L)))
  d$response <- rnorm(n)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 10L, type = "kmeans")
  fit <- sdmTMB(response ~ z + (1 + z + w + v | g), data = d,
    mesh = mesh, spatial = "off", family = gaussian(), do_fit = FALSE)
  p <- fit$tmb_params
  p$re_cov_pars[] <- c(log(0.6), 0.1, 0.2, -0.15, log(0.4),
    0.3, -0.1, log(0.5), 0.25, log(0.7))
  p$re_b_pars[] <- seq(-0.1, 0.1, length.out = length(p$re_b_pars))
  joint_cpp <- TMB::MakeADFun(fit$tmb_data, p, map = fit$tmb_map,
    DLL = "sdmTMB", silent = TRUE)
  joint_rt <- RTMB::MakeADFun(
    rtmb_make_objective(rtmb_prepare(fit$tmb_data)), p,
    map = fit$tmb_map, silent = TRUE)
  expect_equal(joint_rt$fn(), joint_cpp$fn(), tolerance = 1e-7)
  expect_equal(joint_rt$gr(), joint_cpp$gr(),
    ignore_attr = TRUE, tolerance = 1e-6)
  cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random,
    backend = "rtmb")
  expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7)
  expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(rt$report()$eta_iid_re_i, cpp$report()$eta_iid_re_i,
    tolerance = 1e-6)
  projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
  proj_cpp <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
    fit$tmb_random)$report()
  proj_rt <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
    fit$tmb_random, backend = "rtmb")$report()
  expect_equal(proj_rt$proj_eta, proj_cpp$proj_eta, tolerance = 1e-6)
  expect_equal(proj_rt$proj_iid_re_i, proj_cpp$proj_iid_re_i,
    tolerance = 1e-6)
})

test_that("RTMB thresholds and their priors match TMB", {
  set.seed(54)
  d <- data.frame(z = seq(-2, 2, length.out = 60L),
    response = rnorm(60L))
  cases <- list(
    list(term = "breakpt", coefficients = c(0.7, 0.2),
      priors = sdmTMBpriors(
        threshold_breakpt_slope = normal(0, 1),
        threshold_breakpt_cut = normal(0, 1))),
    list(term = "logistic", coefficients = c(-0.1, log(0.8), 1.2),
      priors = sdmTMBpriors(
        threshold_logistic_s50 = normal(0, 1),
        threshold_logistic_s95 = normal(1, 1),
        threshold_logistic_smax = normal(0, 2)))
  )
  for (case in cases) {
    formula <- stats::as.formula(paste0("response ~ 1 + ",
      case$term, "(z)"))
    fit <- sdmTMB(formula, data = d, spatial = "off",
      priors = case$priors, do_fit = FALSE)
    p <- fit$tmb_params
    p$b_threshold[] <- case$coefficients
    cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random)
    rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")
    expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7)
    expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE,
      tolerance = 1e-6)
    if (case$term == "breakpt") {
      off <- cpp$par
      off[tail(which(names(off) == "b_threshold"), 1L)] <- -0.6
      expect_equal(rt$fn(off), cpp$fn(off), tolerance = 1e-7)
      expect_equal(rt$gr(off), cpp$gr(off), ignore_attr = TRUE,
        tolerance = 1e-6)
      expect_equal(as.numeric(cpp$report(off)$eta_i[, 1L]),
        cpp$env$parList(off)$b_j[[1L]] + case$coefficients[[1L]] *
          pmin(d$z, -0.6), tolerance = 1e-7)
    }
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    proj_cpp <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random)$report()
    proj_rt <- make_sdmTMB_adfun(projection, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")$report()
    expect_equal(proj_rt$proj_eta, proj_cpp$proj_eta, tolerance = 1e-6)
    expect_setequal(names(proj_rt), names(proj_cpp))
  }
})

test_that("RTMB Matérn and transformed-parameter priors match TMB", {
  set.seed(55)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    time = rep(1:3, each = 20L), response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 18L, type = "kmeans")
  priors <- sdmTMBpriors(
    matern_s = pc_matern(0.5, 1),
    matern_st = pc_matern(0.5, 1),
    phi = halfnormal(0, 1), ar1_rho = normal(0, 0.7),
    b = normal(c(0, 0), c(2, 2)))
  for (bayesian in c(FALSE, TRUE)) {
    fit <- sdmTMB(response ~ z, data = d, mesh = mesh, time = "time",
      spatial = "on", spatiotemporal = "ar1", priors = priors,
      bayesian = bayesian, do_fit = FALSE)
    p <- fit$tmb_params
    p$b_j[] <- c(0.2, -0.1)
    p$ln_tau_O[] <- 0.2
    p$ln_tau_E[] <- 0.3
    p$ln_kappa[] <- 0.1
    p$ar1_phi[] <- 0.4
    p$ln_phi[] <- log(0.7)
    cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random)
    rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")
    expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7,
      info = paste("bayesian", bayesian))
    expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE,
      tolerance = 1e-6, info = paste("bayesian", bayesian))
  }
})

test_that("RTMB nonstationary epsilon variance matches TMB", {
  set.seed(57)
  n <- 30L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    time = rep(1:3, each = 10L),
    cov = rep(c(-1, 0, 1), each = 10L), response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 10L, type = "kmeans")
  fit <- sdmTMB(response ~ z, data = d, mesh = mesh, time = "time",
    spatial = "on", spatiotemporal = "ar1",
    experimental = list(epsilon_model = "trend", epsilon_predictor = "cov"),
    do_fit = FALSE)
  p <- fit$tmb_params
  p$b_epsilon[] <- 0.3
  p$ln_tau_E[] <- 0.2
  p$ar1_phi[] <- 0.4
  cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random,
    backend = "rtmb")
  expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7)
  expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(rt$report()$sigma_E, cpp$report()$sigma_E, tolerance = 1e-6)
})

test_that("RTMB REML and profiled fixed effects match TMB", {
  set.seed(58)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  d$response <- 0.4 + 0.3 * d$z + rnorm(n, sd = 0.3)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  common <- list(formula = response ~ z, data = d, mesh = mesh,
    spatial = "on", spatiotemporal = "off", reml = TRUE)
  cpp <- do.call(sdmTMB, c(common,
    list(control = sdmTMBcontrol(multiphase = FALSE, backend = "tmb"))))
  rt <- do.call(sdmTMB, c(common,
    list(control = sdmTMBcontrol(multiphase = FALSE, backend = "rtmb"))))
  expect_equal(as.numeric(rt$model$objective),
    as.numeric(cpp$model$objective), tolerance = 1e-6)
  expect_equal(rt$model$par, cpp$model$par, tolerance = 1e-5)
  expect_true(rt$sd_report$pdHess)
  expect_sdreport_matches(rt, cpp)
  expected <- predict(cpp, newdata = d[1:4, ])
  actual <- predict(rt, newdata = d[1:4, ])
  expect_equal(actual$est, expected$est, tolerance = 1e-5)

  prof_cpp <- sdmTMB(response ~ z, data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE, profile = TRUE,
      backend = "tmb"))
  prof_rt <- sdmTMB(response ~ z, data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE, profile = TRUE,
      backend = "rtmb"))
  expect_equal(as.numeric(prof_rt$model$objective),
    as.numeric(prof_cpp$model$objective), tolerance = 1e-6)
  expect_equal(prof_rt$model$par, prof_cpp$model$par, tolerance = 1e-5)
})

test_that("RTMB dispersion formulas match TMB likelihoods and reports", {
  set.seed(60)
  n <- 60L
  d <- data.frame(z = runif(n), w = rnorm(n), response = rnorm(n))
  for (family in list(gaussian(), nbinom1(), nbinom2())) {
    if (family$family != "gaussian") {
      d$response <- rnbinom(n, mu = exp(0.2 + 0.4 * d$z), size = 2)
    }
    fit <- sdmTMB(response ~ z, data = d, spatial = "off",
      family = family, dispformula = ~w, do_fit = FALSE)
    p <- fit$tmb_params
    p$b_disp_k[] <- c(-0.2, 0.3)
    cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random)
    rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map,
      fit$tmb_random, backend = "rtmb")
    expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7)
    expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE,
      tolerance = 1e-6)
    expect_equal(suppressWarnings(rt$report())$phi_i,
      suppressWarnings(cpp$report())$phi_i,
      tolerance = 1e-6)
  }
  single <- sdmTMB(response ~ 1, data = d, spatial = "off",
    family = gaussian(), dispformula = ~0 + z, do_fit = FALSE)
  cpp <- make_sdmTMB_adfun(single$tmb_data, single$tmb_params,
    single$tmb_map, single$tmb_random)
  rt <- make_sdmTMB_adfun(single$tmb_data, single$tmb_params,
    single$tmb_map, single$tmb_random, backend = "rtmb")
  for (obj in list(rt, cpp)) {
    expect_true(all(c("ln_phi_i(0)", "phi_i(0)") %in%
      rownames(summary(sdreport_sdmTMB(obj), "report"))))
  }
})

test_that("RTMB gamma prior on time-varying SD matches TMB", {
  set.seed(56)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n),
    z = seq(-2, 2, length.out = n), time = rep(1:3, each = 20L),
    response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  priors <- sdmTMBpriors(phi = halfnormal(0, 1),
    sigma_V = gamma_cv(0.5, 0.5),
    threshold_breakpt_slope = normal(0, 1),
    threshold_breakpt_cut = normal(0, 1))
  fit <- sdmTMB(response ~ 1 + breakpt(z), data = d, mesh = mesh,
    time = "time", spatial = "off", time_varying = ~1,
    time_varying_type = "rw0", priors = priors, bayesian = TRUE,
    do_fit = FALSE)
  p <- fit$tmb_params
  p$b_threshold[] <- c(0.7, 0.2)
  p$b_rw_t[] <- c(-0.1, 0.2, 0.3)
  p$ln_tau_V[] <- log(0.4)
  cpp <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  rt <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, fit$tmb_random,
    backend = "rtmb")
  expect_equal(rt$fn(), cpp$fn(), tolerance = 1e-7)
  expect_equal(rt$gr(), cpp$gr(), ignore_attr = TRUE,
    tolerance = 1e-6)
})

test_that("RTMB combined model matches every TMB report and sdreport row", {
  set.seed(61)
  n <- 90L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n), w = runif(n),
    v = seq(-2, 2, length.out = n), time = rep(1:3, each = 30L),
    g = factor(rep(1:6, 15L)), response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  fit <- sdmTMB(response ~ 0 + z + s(w, k = 4) + (1 | g) + breakpt(v),
    data = d, mesh = mesh, time = "time", spatial = "on",
    spatiotemporal = "ar1", spatial_varying = ~ 0 + z,
    time_varying = ~1, time_varying_type = "rw0", do_fit = FALSE)
  expect_rtmb_fit_data_matches(fit, d[c(1, 40, 90), ], info = "combined")
})
