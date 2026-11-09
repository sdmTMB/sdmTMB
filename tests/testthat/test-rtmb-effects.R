test_that("RTMB spatially varying coefficients match frozen TMB values", {
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
    expect_frozen_objective(paste("effects svc spatial", spatial),
      fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    expect_equal(ncol(projection$proj_z_i), ncol(fit$tmb_data$z_i))
    expect_frozen_objective(paste("effects svc projection spatial", spatial),
      projection, p, fit$tmb_map, fit$tmb_random, reports = "proj_eta",
      objective = FALSE)
    simulated_data <- fit$tmb_data
    simulated_data$sim_re[[3L]] <- 1L
    sim <- make_sdmTMB_adfun(simulated_data, p, fit$tmb_map,
      fit$tmb_random)$simulate()
    expect_equal(dim(sim$zeta_s), dim(p$zeta_s))
    expect_false(isTRUE(all.equal(sim$zeta_s, p$zeta_s)))
    expect_equal(as.numeric(sim$eta_i),
      as.numeric(simulated_data$X_ij[[1L]] %*% p$b_j +
        simulated_data$offset_i +
        rowSums(sim$zeta_s_A[, , 1L] * simulated_data$z_i)),
      tolerance = 1e-6)
  }
})

test_that("RTMB time-varying coefficients match frozen TMB values", {
  set.seed(50)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    w = rnorm(n), time = rep(1:3, each = 20L))
  d$response <- rnorm(n)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 18L, type = "kmeans")
  for (time_type in c("rw", "rw0", "ar1")) {
    fit <- sdmTMB(response ~ 1, data = d, mesh = mesh, time = "time",
      spatial = "off", spatiotemporal = "off",
      time_varying = ~0 + z + w, time_varying_type = time_type,
      family = gaussian(), control = sdmTMBcontrol(multiphase = FALSE),
      do_fit = FALSE)
    p <- fit$tmb_params
    p$b_rw_t[] <- seq(-0.2, 0.2, length.out = length(p$b_rw_t))
    p$ln_tau_V[] <- log(0.6)
    p$rho_time_unscaled[] <- 0.4
    expect_frozen_objective(paste("effects time-varying", time_type),
      fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    expect_frozen_objective(paste("effects time-varying projection",
      time_type), projection, p, fit$tmb_map, fit$tmb_random,
      reports = c("proj_eta", "proj_rw_i"), objective = FALSE)
    simulated_data <- fit$tmb_data
    simulated_data$sim_re[[5L]] <- 1L
    simulated_data$simulate_t <- c(0L, 1L, 1L)
    sim <- make_sdmTMB_adfun(simulated_data, p, fit$tmb_map,
      fit$tmb_random)$simulate()
    expect_equal(sim$b_rw_t[1L, , 1L], p$b_rw_t[1L, , 1L])
    expect_false(isTRUE(all.equal(sim$b_rw_t[2L, , 1L],
      p$b_rw_t[2L, , 1L])))
  }
})

test_that("RTMB correlated IID effects match frozen TMB values", {
  set.seed(51)
  n <- 40L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    w = rnorm(n), v = rnorm(n), g = factor(rep(1:5, each = 8L)))
  d$response <- rnorm(n)
  mesh <- make_mesh(d, c("x", "y"), cutoff = 0.1)
  fit <- sdmTMB(response ~ z + (1 + z + w + v | g), data = d,
    mesh = mesh, spatial = "off", family = gaussian(), do_fit = FALSE)
  p <- fit$tmb_params
  p$re_cov_pars[] <- c(log(0.6), 0.1, 0.2, -0.15, log(0.4),
    0.3, -0.1, log(0.5), 0.25, log(0.7))
  p$re_b_pars[] <- seq(-0.1, 0.1, length.out = length(p$re_b_pars))
  expect_frozen_objective("effects correlated iid", fit$tmb_data, p,
    fit$tmb_map, fit$tmb_random)
  projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
  expect_frozen_objective("effects correlated iid projection", projection,
    p, fit$tmb_map, fit$tmb_random, reports = c("proj_eta", "proj_iid_re_i"),
    objective = FALSE)
})

test_that("RTMB thresholds and their priors match frozen TMB values", {
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
    rt <- expect_frozen_objective(paste("effects threshold", case$term),
      fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
    if (case$term == "breakpt") {
      off <- rt$par
      off[tail(which(names(off) == "b_threshold"), 1L)] <- -0.6
      expect_equal(as.numeric(rt$report(off)$eta_i[, 1L]),
        rt$env$parList(off)$b_j[[1L]] + case$coefficients[[1L]] *
          pmin(d$z, -0.6), tolerance = 1e-7)
    }
    projection <- predict(fit, newdata = d[1:5, ], return_tmb_data = TRUE)
    expect_frozen_objective(paste("effects threshold projection",
      case$term), projection, p, fit$tmb_map, fit$tmb_random,
      reports = "proj_eta", objective = FALSE)
  }
})

test_that("RTMB Matérn and transformed-parameter priors match frozen TMB values", {
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
    expect_frozen_objective(paste("effects matern priors bayesian",
      bayesian), fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  }
})

test_that("RTMB REML and profiled fixed effects fit", {
  set.seed(58)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  d$response <- 0.4 + 0.3 * d$z + rnorm(n, sd = 0.3)
  mesh <- suppressWarnings(
    make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans"),
    classes = "sdmTMB_mesh_outside_warning"
  )
  rt <- sdmTMB(response ~ z, data = d, mesh = mesh, spatial = "on",
    spatiotemporal = "off", reml = TRUE,
    control = sdmTMBcontrol(multiphase = FALSE))
  expect_true(rt$sd_report$pdHess)
  expect_frozen_objective("effects reml", rt$tmb_data,
    rtmb_test_parameters(rt$tmb_params, rt$tmb_map), rt$tmb_map,
    rt$tmb_random)
  ml <- sdmTMB(response ~ z, data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE))
  prof <- sdmTMB(response ~ z, data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE, profile = TRUE))
  expect_equal(as.numeric(prof$model$objective),
    as.numeric(ml$model$objective), tolerance = 1e-6)
  expect_equal(tidy(prof)$estimate, tidy(ml)$estimate, tolerance = 1e-5)
})

test_that("RTMB dispersion formulas match frozen TMB values", {
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
    suppressWarnings(expect_frozen_objective(paste("effects dispformula",
      family$family), fit$tmb_data, p, fit$tmb_map, fit$tmb_random))
  }
  single <- sdmTMB(response ~ 1, data = d, spatial = "off",
    family = gaussian(), dispformula = ~0 + z, do_fit = FALSE)
  rt <- make_sdmTMB_adfun(single$tmb_data, single$tmb_params,
    single$tmb_map, single$tmb_random)
  expect_true(all(c("ln_phi_i(0)", "phi_i(0)") %in%
    rownames(summary(RTMB::sdreport(rt), "report"))))
})

test_that("RTMB gamma prior on time-varying SD matches frozen TMB values", {
  set.seed(56)
  n <- 60L
  d <- data.frame(x = runif(n), y = runif(n),
    z = seq(-2, 2, length.out = n), time = rep(1:3, each = 20L),
    response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
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
  expect_frozen_objective("effects gamma prior sigma_V", fit$tmb_data, p,
    fit$tmb_map, fit$tmb_random)
})

test_that("RTMB combined model matches every TMB report and sdreport row", {
  set.seed(61)
  n <- 90L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n), w = runif(n),
    v = seq(-2, 2, length.out = n), time = rep(1:3, each = 30L),
    g = factor(rep(1:6, 15L)), response = rnorm(n))
  mesh <- make_mesh(d, c("x", "y"), cutoff = 0.1)
  fit <- sdmTMB(response ~ 0 + z + s(w, k = 4) + (1 | g) + breakpt(v),
    data = d, mesh = mesh, time = "time", spatial = "on",
    spatiotemporal = "ar1", spatial_varying = ~ 0 + z,
    time_varying = ~1, time_varying_type = "rw0", do_fit = FALSE)
  expect_frozen_fit_data("effects", fit, d[c(1, 40, 90), ], info = "combined")
})
