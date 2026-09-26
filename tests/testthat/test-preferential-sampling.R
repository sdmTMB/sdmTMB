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

  # A supported specification is prepared, then stopped before any objective
  # is taped without the sampling likelihood.
  expect_error(build(pref_formula, offset = log(dat$effort)), "not implemented yet")
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
  expect_error(
    make_sdmTMB_adfun(data, fit$tmb_params, fit$tmb_map, fit$tmb_random,
      backend = "rtmb"),
    "not implemented yet"
  )
})

test_that("update() keeps the preferential specification and rejects TMB", {
  skip_on_cran()
  fit <- pref_fit()
  spec <- preferential_sampling(sampled ~ 1, pref_grid())
  # Stand-in for a fitted preferential model until the likelihood exists.
  fit$preferential <- list(spec = spec)
  call <- update(fit, evaluate = FALSE)
  expect_identical(call$preferential, spec)
  expect_identical(call$control$backend, "rtmb")
  expect_error(
    update(fit, offset = log(fit$data$effort),
      control = sdmTMBcontrol(backend = "tmb")),
    "requires backend = \"rtmb\"", fixed = TRUE
  )
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
