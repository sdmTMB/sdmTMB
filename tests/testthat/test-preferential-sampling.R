pref_dat <- function() {
  set.seed(1)
  n <- 300
  dat <- data.frame(
    x = runif(n, 0, 10), y = runif(n, 0, 10),
    year = rep(2018:2020, length.out = n), depth = runif(n, 10, 100),
    gear = factor(sample(c("a", "b", "c"), n, TRUE)),
    vessel = factor(sample(letters[1:5], n, TRUE)), effort = runif(n, 0.5, 2)
  )
  dat$catch <- rpois(n, dat$effort * exp(0.5 + 0.01 * dat$depth +
    c(0, 0.4, -0.4)[dat$gear] + c(-0.6, -0.3, 0, 0.3, 0.6)[dat$vessel] +
    sin(dat$x / 2) * cos(dat$y / 3) + 0.3 * sin(dat$x + dat$year)))
  dat
}

# Eligible frame: every cell in every year, with one standardized gear.
pref_grid <- function(seed = 2) {
  set.seed(seed)
  grid <- expand.grid(x = seq(1.5, 8.5, 1), y = seq(1.5, 8.5, 1),
    year = 2018:2020)
  grid$depth <- runif(nrow(grid), 10, 100)
  grid$gear <- factor("b", levels = c("a", "b", "c"))
  grid$sampled <- rbinom(nrow(grid), 1, 0.3)
  grid
}

pref_mesh <- function(dat) make_mesh(dat, c("x", "y"), cutoff = 1.5)

pref_formula <- catch ~ poly(depth, 2) + gear + (1 | vessel)

pref_fit <- fit_once(function() {
  dat <- pref_dat()
  sdmTMB(pref_formula,
    data = dat, mesh = pref_mesh(dat), time = "year", family = poisson(),
    offset = log(dat$effort), control = sdmTMBcontrol(backend = "rtmb")
  )
})

# A coupled simulation: the sampling indicators depend on the same latent
# catch surface that generates catches at the sampled cells. `xi_sd > 0` adds
# a sampling-only field. Catches exist only at sampled cells.
pref_sim <- function(b = 0.8, xi_sd = 0, seed = 1) {
  old <- options(sdmTMB.backend = "tmb")
  on.exit(options(old))
  set.seed(seed)
  grid <- expand.grid(x = seq(0.25, 9.75, 0.5), y = seq(0.25, 9.75, 0.5),
    year = 1:3)
  mesh <- make_mesh(grid, c("x", "y"), cutoff = 1)
  sim <- sdmTMB_simulate(~1, data = grid, mesh = mesh, time = "year",
    family = poisson(), range = 4, sigma_O = 0.7, sigma_E = 0.3,
    B = log(3), seed = seed, spatiotemporal = "iid")
  xi <- 0
  if (xi_sd > 0) {
    cells <- grid[grid$year == 1L, ]
    xi_sim <- sdmTMB_simulate(~1, data = cells,
      mesh = make_mesh(cells, c("x", "y"), mesh = mesh$mesh),
      family = gaussian(), range = 3, sigma_O = xi_sd, phi = 0.1, B = 0,
      seed = seed + 1L)
    xi <- rep(xi_sim$omega_s, 3L)
  }
  gamma <- c(-2, -1.5, -1)
  grid$sampled <- stats::rbinom(nrow(grid), 1,
    stats::plogis(gamma[grid$year] + b * sim$eta + xi))
  dat <- grid[grid$sampled == 1, ]
  dat$catch <- stats::rpois(nrow(dat), exp(sim$eta[grid$sampled == 1]))
  list(grid = grid, dat = dat,
    mesh = make_mesh(dat, c("x", "y"), mesh = mesh$mesh))
}

pref_joint <- function(sim, spatial = "off", control = list(), ...) {
  sdmTMB(catch ~ 1, data = sim$dat, mesh = sim$mesh, time = "year",
    family = poisson(),
    preferential = preferential_sampling(sampled ~ 0 + factor(year),
      data = sim$grid, spatial = spatial),
    control = do.call(sdmTMBcontrol, c(list(backend = "rtmb"), control)), ...
  )
}

pref_joint_fit <- fit_once(function() pref_joint(pref_sim()))

pref_names <- c("gamma_pref", "b_pref", "ln_tau_xi", "ln_kappa_xi", "xi_s")

prepare_for <- function(fit, spec) {
  .prepare_preferential(spec, fit$terms, fit$xlevels, fit$contrasts,
    fit$tmb_data$X_ij, fit$spde, fit$time, fit$time_lu)
}

# Shared catch predictor at the frame rows, from the RTMB predictor code at
# the fitted parameters and latent effects.
shared_predictor <- function(fit, prep) {
  data <- fit$tmb_data
  data$preferential <- prep$data
  prepared <- rtmb_prepare(data)
  par <- fit$tmb_obj$env$parList(par = fit$tmb_obj$env$last.par.best)
  theta <- rtmb_transform(par, prepared)
  effects <- rtmb_latent_effects(par, theta, prepared, character(0))
  rtmb_linear_predictors(par, theta, effects, prepared,
    prepared$preferential$rows)
}

test_that("preferential_sampling() validates its specification", {
  grid <- pref_grid()
  spec <- preferential_sampling(sampled ~ 0 + factor(year), data = grid)
  expect_s3_class(spec, "sdmTMB_preferential")
  expect_identical(spec$response, "sampled")
  expect_false(spec$include_iid)
  expect_identical(spec$offset, rep(0, nrow(grid)))
  expect_identical(spec$spatial, "off")
  expect_true(preferential_sampling(sampled ~ 1, grid, re_form_iid = NULL)$include_iid)

  expect_error(preferential_sampling(~1, grid), "two-sided")
  expect_error(preferential_sampling(log(sampled) ~ 1, grid), "left side")
  expect_error(preferential_sampling(missing ~ 1, grid), "response column")
  expect_error(preferential_sampling(sampled ~ s(depth), grid), "Smoothers")
  expect_error(preferential_sampling(sampled ~ (1 | gear), grid), "Random effects")
  expect_error(preferential_sampling(sampled ~ offset(depth), grid), "offset")
  expect_error(preferential_sampling(sampled ~ 1, as.list(grid)), "data frame")
  expect_error(preferential_sampling(sampled ~ 1, grid[0, ]), "no rows")
  bad <- grid
  bad$sampled[1] <- 2
  expect_error(preferential_sampling(sampled ~ 1, bad), "0, 1, or `NA`")
  expect_error(preferential_sampling(sampled ~ 1, grid, re_form_iid = ~0), "re_form_iid")
  expect_error(preferential_sampling(sampled ~ 1, grid, offset = c(1, 2)), "offset")
  expect_error(preferential_sampling(sampled ~ 1, grid, offset = NA_real_), "offset")
  expect_error(preferential_sampling(sampled ~ 1, grid, coefficient = "rw"), "constant")
  expect_error(preferential_sampling(sampled ~ 1, grid, spatial = "maybe"))
})

test_that("preferential sampling requires RTMB and rejects unsupported features", {
  skip_on_cran()
  dat <- pref_dat()
  mesh <- pref_mesh(dat)
  spec <- preferential_sampling(sampled ~ 0 + factor(year), data = pref_grid())
  build <- function(formula = catch ~ depth, family = poisson(),
                    backend = "rtmb", preferential = spec, ...) {
    sdmTMB(formula,
      data = dat, mesh = mesh, time = "year", family = family,
      preferential = preferential,
      control = sdmTMBcontrol(backend = backend), do_fit = FALSE, ...
    )
  }
  expect_error(build(backend = "tmb"),
    "requires backend = \"rtmb\"; set control = sdmTMBcontrol(backend = \"rtmb\")",
    fixed = TRUE)
  expect_error(build(preferential = list()), "preferential_sampling()", fixed = TRUE)
  expect_error(build(family = delta_gamma()), "delta models")
  expect_error(build(catch ~ s(depth)), "smoothers")
  expect_error(build(catch ~ breakpt(depth)), "threshold")
  expect_error(build(time_varying = ~depth), "time_varying")
  expect_error(build(spatial_varying = ~depth), "spatial_varying")
  expect_error(build(anisotropy = TRUE), "anisotropy")
  expect_error(build(family = gaussian()), "families other than")
  expect_error(build(family = Gamma(link = "inverse")), "families other than")
  expect_error(build(family = binomial()), "families other than")
  expect_error(build(spatial = "off", spatiotemporal = "off"), "spatial or spatiotemporal field")
  spec_iid <- preferential_sampling(sampled ~ 0 + factor(year), data = pref_grid(),
    re_form_iid = NULL)
  expect_error(build(catch ~ depth + (1 | vessel), preferential = spec_iid), "IID")
  expect_error(
    sdmTMB(catch ~ depth, data = dat, time = "year", family = poisson(),
      spatial = "off", preferential = spec,
      control = sdmTMBcontrol(backend = "rtmb"), do_fit = FALSE),
    "mesh"
  )

  # A supported specification builds with its sampling parameters. The
  # preference coefficient waits for the catch fields in multiphase fits.
  obj <- build(pref_formula, offset = log(dat$effort))
  expect_length(obj$tmb_params$gamma_pref, 3L)
  expect_identical(obj$tmb_params$b_pref, 0)
  expect_false(any(c("ln_tau_xi", "ln_kappa_xi", "xi_s") %in%
    names(obj$tmb_params)))
  expect_false("xi_s" %in% obj$tmb_random)
  expect_false(any(c("gamma_pref", "b_pref") %in% names(obj$tmb_map)))
  spec_xi <- preferential_sampling(sampled ~ 0 + factor(year),
    data = pref_grid(), spatial = "on")
  obj_xi <- build(pref_formula, offset = log(dat$effort),
    preferential = spec_xi)
  expect_true("xi_s" %in% obj_xi$tmb_random)
  expect_length(obj_xi$tmb_params$xi_s, mesh$mesh$n)
})

test_that("the TMB backend refuses preferential model data", {
  skip_on_cran()
  fit <- pref_fit()
  prep <- prepare_for(fit, preferential_sampling(sampled ~ 1, pref_grid()))
  data <- fit$tmb_data
  data$preferential <- prep$data
  expect_error(
    make_sdmTMB_adfun(data, fit$tmb_params, fit$tmb_map, fit$tmb_random,
      backend = "tmb"),
    "requires backend = \"rtmb\"", fixed = TRUE
  )
  obj <- make_sdmTMB_adfun(data, c(fit$tmb_params, prep$parameters),
    fit$tmb_map, fit$tmb_random, backend = "rtmb")
  expect_true(all(c("gamma_pref", "b_pref") %in% names(obj$par)))
})

test_that("update() keeps the preferential specification and rejects TMB", {
  skip_on_cran()
  fit <- pref_joint_fit()
  spec <- fit$preferential$spec
  call <- update(fit, evaluate = FALSE)
  expect_identical(call$preferential, spec)
  expect_identical(call$control$backend, "rtmb")
  expect_error(update(fit, control = sdmTMBcontrol(backend = "tmb")),
    "requires backend = \"rtmb\"", fixed = TRUE)
  refit <- update(fit)
  expect_equal(refit$model$objective, fit$model$objective, tolerance = 1e-6)
})

test_that("the shared surface matches ordinary RTMB prediction on the frame", {
  skip_on_cran()
  fit <- pref_fit()
  grid <- pref_grid()
  grid$sampled[grid$year == 2020][1:10] <- NA
  spec <- preferential_sampling(sampled ~ 0 + factor(year), data = grid,
    offset = 0.3)
  prep <- prepare_for(fit, spec)
  lp <- shared_predictor(fit, prep)

  # The frame needs no `vessel` column: IID effects are excluded. predict()
  # still wants one.
  expect_false("vessel" %in% names(grid))
  nd <- grid
  nd$vessel <- factor("a", levels = levels(fit$data$vessel))
  p <- predict(fit, newdata = nd, re_form_iid = NA, offset = rep(0.3, nrow(nd)))
  expect_equal(lp$eta[, 1], p$est, tolerance = 1e-10)
  expect_equal(lp$fixed[, 1] + 0.3, p$est_non_rf, tolerance = 1e-10)
  expect_equal(lp$omega[, 1] + lp$epsilon[, 1], p$est_rf, tolerance = 1e-10)
  expect_equal(unname(prep$data$X_ij[[1]]),
    unname(predict(fit, newdata = nd, offset = rep(0.3, nrow(nd)),
      return_tmb_data = TRUE)$proj_X_ij[[1]]))

  # Row order is preserved in every derived input.
  perm <- sample(nrow(grid))
  prep_perm <- prepare_for(fit, preferential_sampling(sampled ~ 0 + factor(year),
    data = grid[perm, ], offset = 0.3))
  expect_equal(shared_predictor(fit, prep_perm)$eta[, 1], lp$eta[perm, 1],
    tolerance = 1e-10)
  expect_identical(prep_perm$data$R_i, prep$data$R_i[perm])
  expect_equal(prep_perm$data$Z_ij, prep$data$Z_ij[perm, ], ignore_attr = TRUE)
  expect_identical(prep_perm$data$year_i, prep$data$year_i[perm])

  # Changing a main-only covariate changes the shared design only.
  grid2 <- grid
  grid2$gear[] <- "c"
  prep2 <- prepare_for(fit, preferential_sampling(sampled ~ 0 + factor(year),
    data = grid2, offset = 0.3))
  expect_identical(prep2$data$Z_ij, prep$data$Z_ij)
  expect_equal(prep2$data$X_ij[[1]][, "gearc"], rep(1, nrow(grid)),
    ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(shared_predictor(fit, prep2)$eta, lp$eta)))
})

test_that("preferential preparation records the frame and missing history", {
  skip_on_cran()
  fit <- pref_fit()
  grid <- pref_grid()
  grid$sampled[grid$year == 2020][1:10] <- NA
  prep <- prepare_for(fit, preferential_sampling(sampled ~ 0 + factor(year),
    data = grid, spatial = "on"))
  d <- prep$data
  expect_identical(d$n_pref, nrow(grid))
  expect_identical(d$year_i, as.integer(grid$year - 2018L))
  expect_identical(d$spatial_xi, 1L)
  expect_identical(d$include_iid, 0L)
  # Unique locations are projected once.
  expect_identical(nrow(d$A_station), 64L)
  expect_identical(ncol(d$A_station), fit$spde$mesh$n)
  expect_identical(max(d$station_i), 63L)
  expect_identical(prep$info$n_unknown, 10L)
  expect_identical(prep$info$n_observed, nrow(grid) - 10L)
  expect_identical(prep$info$n_sampled, sum(grid$sampled, na.rm = TRUE))
  expect_identical(colnames(d$Z_ij), paste0("factor(year)", 2018:2020))
  expect_s3_class(prep$info$sampling_terms, "terms")

  data <- fit$tmb_data
  data$preferential <- d
  pref <- rtmb_prepare(data)$preferential
  expect_identical(pref$observed, which(!is.na(grid$sampled)))
  expect_identical(pref$rows$time, d$year_i + 1L)
  expect_true(pref$xi)

  # A missing sampling year is fine: biological time steps are unaffected.
  prep_missing <- prepare_for(fit, preferential_sampling(
    sampled ~ 0 + factor(year), data = grid[grid$year != 2019, ]))
  expect_identical(sort(unique(prep_missing$data$year_i)), c(0L, 2L))

  # So is an unobserved `extra_time` step.
  dat <- pref_dat()
  fit_extra <- sdmTMB(catch ~ depth, data = dat, mesh = pref_mesh(dat),
    time = "year", extra_time = 2021L, family = poisson(), do_fit = FALSE,
    control = sdmTMBcontrol(backend = "rtmb"))
  grid_extra <- rbind(grid, transform(grid[grid$year == 2020, ], year = 2021L,
    sampled = NA))
  prep_extra <- prepare_for(fit_extra,
    preferential_sampling(sampled ~ 1, data = grid_extra))
  expect_identical(max(prep_extra$data$year_i), 3L)
  expect_identical(sum(is.na(prep_extra$data$R_i)), 74L)
})

test_that("preferential preparation rejects invalid frames", {
  skip_on_cran()
  fit <- pref_fit()
  grid <- pref_grid()
  prep <- function(data, formula = sampled ~ 0 + factor(year), ...) {
    prepare_for(fit, preferential_sampling(formula, data = data, ...))
  }

  bad <- grid
  bad$sampled[bad$year == 2019] <- NA
  expect_error(prep(bad), "not full rank")
  expect_error(prep(bad), "factor(year)2019", fixed = TRUE)
  expect_silent(prep(bad, sampled ~ 1))

  bad <- grid
  bad$sampled[bad$year == 2019] <- 0
  expect_error(prep(bad), "all 0 wherever column `factor(year)2019` is 1",
    fixed = TRUE)
  bad$sampled[] <- 0
  expect_error(prep(bad, sampled ~ 1), "both 0s and 1s")

  expect_error(prep(grid[, names(grid) != "depth"]), "Missing: depth")
  expect_error(prep(grid[, names(grid) != "year"]), "time")
  bad <- grid
  bad$depth[3] <- NA
  expect_error(prep(bad), "can't contain `NA`")
  bad$depth[3] <- Inf
  expect_error(prep(bad), "finite")
  bad <- grid
  bad$year[1] <- 2030L
  expect_error(prep(bad), "outside the fitted time steps")
  expect_error(prep(rbind(grid, grid[1, ])), "one row per cell")
  bad <- grid
  bad$x[5] <- 50
  expect_error(prep(bad), "outside the mesh")
  bad <- grid
  bad$gear <- factor("d")
  expect_error(prep(bad), "new level")
  bad <- grid
  bad$dist <- 1
  expect_error(prep(bad, sampled ~ 0 + factor(year) + dist), "not full rank")
  expect_error(prep(grid, sampled ~ missing_covariate), "missing_covariate")
})

test_that("the shared design reuses fitted bases and rejects recomputed ones", {
  skip_on_cran()
  # The retired builder recomputed poly() on the frame; the shared design
  # uses the fitted basis in `predvars`.
  dat <- data.frame(x = seq(0, 9), y = seq(0, 9), year = 1L, depth = 1:10,
    catch = rpois(10, 5))
  grid <- data.frame(x = c(1, 2, 3, 4), y = c(1, 2, 3, 4), year = 1L,
    depth = c(2, 4, 7, 9), sampled = c(0, 1, 0, 1))
  mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 1)
  fit <- sdmTMB(catch ~ poly(depth, 2), data = dat, mesh = mesh,
    time = "year", family = poisson(), do_fit = FALSE,
    control = sdmTMBcontrol(backend = "rtmb"))
  prep <- prepare_for(fit, preferential_sampling(sampled ~ 1, grid))
  fitted_terms <- stats::terms(stats::model.frame(catch ~ poly(depth, 2), dat))
  expected <- stats::model.matrix(stats::delete.response(fitted_terms), grid)
  expect_equal(prep$data$X_ij[[1]], expected, ignore_attr = TRUE)
  recomputed <- stats::model.matrix(~ poly(depth, 2), grid)
  expect_gt(max(abs(prep$data$X_ij[[1]] - recomputed)), 0.1)

  fit_scale <- sdmTMB(catch ~ scale(depth), data = dat, mesh = mesh,
    time = "year", family = poisson(), do_fit = FALSE,
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_error(prepare_for(fit_scale, preferential_sampling(sampled ~ 1, grid)),
    "Precompute")
})

# Objective with every latent effect held fixed (nothing integrated), so the
# joint and catch-only objectives differ by exactly the sampling terms.
fixed_latent_objectives <- function(fit, par) {
  keep <- setdiff(names(par), pref_names)
  map <- fit$tmb_map
  joint <- make_sdmTMB_adfun(fit$tmb_data, par, map, random = NULL,
    backend = "rtmb")
  data <- fit$tmb_data
  data$preferential <- NULL
  catch <- make_sdmTMB_adfun(data, par[keep], map[intersect(names(map), keep)],
    random = NULL, backend = "rtmb")
  list(joint = joint, catch = catch)
}

# Negative log density of a zero-mean GMRF with precision Q, computed densely.
dense_gmrf_nll <- function(x, Q) {
  Q <- as.matrix(Q)
  0.5 * length(x) * log(2 * pi) -
    0.5 * as.numeric(determinant(Q)$modulus) + 0.5 * sum(x * (Q %*% x))
}

test_that("the joint objective adds the Bernoulli and sampling-field terms", {
  skip_on_cran()
  sim <- pref_sim(xi_sd = 0.5)
  fit <- pref_joint(sim, spatial = "on", do_fit = FALSE)
  par <- fit$tmb_params
  set.seed(3)
  par$b_j <- 1.1
  par$omega_s[] <- stats::rnorm(length(par$omega_s), 0, 0.5)
  par$epsilon_st[] <- stats::rnorm(length(par$epsilon_st), 0, 0.3)
  par$gamma_pref <- c(-2, -1.4, -0.9)
  par$ln_tau_xi <- 0.4
  par$ln_kappa_xi <- -0.2
  par$xi_s <- stats::rnorm(length(par$xi_s), 0, 0.4)
  par$ln_tau_O <- -0.5
  par$ln_tau_E <- 0.2
  par$ln_kappa[] <- -0.3

  # Independent shared surface and sampling logits at the frame rows.
  grid <- sim$grid
  A <- fmesher::fm_basis(sim$mesh$mesh, loc = as.matrix(grid[, c("x", "y")]))
  eps <- as.matrix(A %*% par$epsilon_st[, , 1])
  h <- par$b_j + as.vector(A %*% par$omega_s[, 1]) +
    eps[cbind(seq_len(nrow(grid)), grid$year)]
  xi <- as.vector(A %*% par$xi_s)
  Q <- with(fit$tmb_data$spde, exp(par$ln_kappa_xi)^4 * M0 +
    2 * exp(par$ln_kappa_xi)^2 * M1 + M2)
  xi_nll <- dense_gmrf_nll(par$xi_s, exp(2 * par$ln_tau_xi) * Q)

  # Include extreme logits: a naive log(1 - plogis(eta)) underflows there.
  for (b in c(0.7, 40)) {
    par$b_pref <- b
    eta <- par$gamma_pref[grid$year] + b * h + xi
    bern_nll <- -sum(ifelse(grid$sampled == 1,
      stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE)))
    obj <- fixed_latent_objectives(fit, par)
    expect_equal(obj$joint$fn(obj$joint$par) - obj$catch$fn(obj$catch$par),
      bern_nll + xi_nll, tolerance = 1e-8)
    r <- obj$joint$report(obj$joint$par)
    expect_equal(r$sampling_target_i, h, tolerance = 1e-10)
    expect_equal(r$sampling_eta_i, eta, tolerance = 1e-10)
    expect_equal(r$sampling_field_i, xi, tolerance = 1e-10)
    expect_equal(r$sampling_fixed_i + r$sampling_preference_i +
      r$sampling_field_i, r$sampling_eta_i, tolerance = 1e-10)
  }
  expect_true(any(abs(eta) > 50))

  # Unknown indicators contribute nothing.
  par$b_pref <- 0.7
  eta <- par$gamma_pref[grid$year] + 0.7 * h + xi
  known <- grid$year != 2L | seq_len(nrow(grid)) %% 3 != 0
  grid_na <- grid
  grid_na$sampled[!known] <- NA
  fit_na <- pref_joint(list(grid = grid_na, dat = sim$dat, mesh = sim$mesh),
    spatial = "on", do_fit = FALSE)
  obj <- fixed_latent_objectives(fit_na, par)
  bern_nll <- -sum(stats::dbinom(grid$sampled[known], 1,
    stats::plogis(eta[known]), log = TRUE))
  expect_equal(obj$joint$fn(obj$joint$par) - obj$catch$fn(obj$catch$par),
    bern_nll + xi_nll, tolerance = 1e-8)

  # AD gradients of the sampling parameters and shared catch parameters
  # agree with central differences.
  obj <- fixed_latent_objectives(fit, par)$joint
  x <- obj$par
  check <- c(which(names(x) %in% c("b_j", "ln_tau_O", "ln_kappa",
    "gamma_pref", "b_pref", "ln_tau_xi", "ln_kappa_xi")),
    which(names(x) == "omega_s")[c(1, 10)], which(names(x) == "xi_s")[c(2, 20)])
  numeric_gr <- vapply(check, function(i) {
    step <- 1e-5
    up <- down <- x
    up[i] <- x[i] + step
    down[i] <- x[i] - step
    (obj$fn(up) - obj$fn(down)) / (2 * step)
  }, numeric(1))
  expect_equal(as.vector(obj$gr(x))[check], numeric_gr, tolerance = 1e-5)
})

test_that("ordinary models get no preferential parameters on either backend", {
  skip_on_cran()
  sim <- pref_sim()
  for (backend in c("tmb", "rtmb")) {
    fit <- sdmTMB(catch ~ 1, data = sim$dat, mesh = sim$mesh, time = "year",
      family = poisson(), do_fit = FALSE,
      control = sdmTMBcontrol(backend = backend))
    expect_false(any(pref_names %in% names(fit$tmb_params)))
    expect_null(fit$tmb_data$preferential)
  }
})

test_that("a coupled simulation recovers the preference coefficient", {
  skip_on_cran()
  fit <- pref_joint_fit()
  expect_true(fit$pos_def_hessian)
  expect_lt(max(abs(fit$gradients)), 1e-4)
  sdr <- summary(fit$sd_report, "fixed")
  b <- sdr[rownames(sdr) == "b_pref", ]
  expect_lt(abs(b[["Estimate"]] - 0.8), 3 * b[["Std. Error"]])
  expect_gt(b[["Estimate"]] / b[["Std. Error"]], 3)
  expect_identical(sum(rownames(sdr) == "gamma_pref"), 3L)
  expect_false(any(c("ln_tau_xi", "xi_s") %in% names(fit$tmb_obj$env$par)))

  # Reports are by frame row, in the frame's order.
  r <- fit$tmb_obj$report(fit$tmb_obj$env$last.par.best)
  expect_length(r$sampling_p_i, nrow(fit$preferential$spec$data))
  expect_equal(r$sampling_p_i, stats::plogis(r$sampling_eta_i))
  expect_identical(fit$preferential$n_sampled,
    sum(fit$preferential$spec$data$sampled))

  # A prediction retape keeps the joint likelihood.
  lifecycle::expect_deprecated(
    p <- predict(fit, newdata = pref_sim()$grid, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  expect_equal(p$obj$fn(fit$model$par), fit$model$objective, tolerance = 1e-6)

  # Permuting the frame rows doesn't change the likelihood.
  grid <- fit$preferential$spec$data
  sim <- list(grid = grid[sample(nrow(grid)), ], dat = fit$data,
    mesh = fit$spde)
  perm <- pref_joint(sim, do_fit = FALSE)
  expect_equal(perm$tmb_obj$fn(fit$model$par), fit$model$objective,
    tolerance = 1e-6)
})

test_that("a sampling-only field can be estimated", {
  skip_on_cran()
  fit <- pref_joint(pref_sim(xi_sd = 0.8), spatial = "on")
  expect_true(fit$pos_def_hessian)
  expect_true("xi_s" %in% fit$tmb_random)
  sdr <- summary(fit$sd_report, "report")
  expect_true(all(c("sigma_xi", "range_xi") %in% rownames(sdr)))
  expect_gt(sdr["sigma_xi", "Estimate"], 0.3)
})

test_that("with the preference coefficient fixed at 0 the models decouple", {
  skip_on_cran()
  sim <- pref_sim()
  joint <- pref_joint(sim,
    control = list(map = list(b_pref = factor(NA)), start = list(b_pref = 0)))
  catch <- sdmTMB(catch ~ 1, data = sim$dat, mesh = sim$mesh, time = "year",
    family = poisson(), control = sdmTMBcontrol(backend = "rtmb"))
  shared <- names(catch$model$par)
  expect_equal(joint$model$par[names(joint$model$par) %in% shared],
    catch$model$par, tolerance = 1e-4)
  sampling <- stats::glm(sampled ~ 0 + factor(year), family = binomial(),
    data = sim$grid)
  expect_equal(unname(joint$model$par[names(joint$model$par) == "gamma_pref"]),
    unname(stats::coef(sampling)), tolerance = 1e-4)
})

test_that("a spatial-only model without time fits", {
  skip_on_cran()
  sim <- pref_sim()
  grid <- sim$grid[sim$grid$year == 1L, names(sim$grid) != "year"]
  dat <- sim$dat[sim$dat$year == 1L, ]
  fit <- sdmTMB(catch ~ 1, data = dat,
    mesh = make_mesh(dat, c("x", "y"), mesh = sim$mesh$mesh),
    family = poisson(),
    preferential = preferential_sampling(sampled ~ 1, data = grid),
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_true(fit$pos_def_hessian)
  expect_true("b_pref" %in% names(fit$model$par))
})
