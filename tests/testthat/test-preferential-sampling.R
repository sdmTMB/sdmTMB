test_that("preferential-sampling R-side wiring validates and assembles tmb_data", {
  skip_on_cran()
  set.seed(1)

  n <- 60
  dat <- data.frame(
    x = runif(n, 0, 10),
    y = runif(n, 0, 10),
    year = rep(2018:2020, length.out = n),
    region = factor(sample(c("A", "B"), n, replace = TRUE)),
    catch = rpois(n, 5)
  )
  mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 2)

  grid_xy <- expand.grid(x = seq(0, 10, length.out = 5), y = seq(0, 10, length.out = 3))
  pref_grid <- do.call(rbind, lapply(2018:2020, function(yr) {
    g <- grid_xy
    g$year <- yr
    g$region <- factor(sample(c("A", "B"), nrow(g), replace = TRUE), levels = c("A", "B"))
    g$sampled <- rbinom(nrow(g), 1, 0.5)
    g
  }))

  # happy path: preferential_formula explicitly asks for the full `region`
  # term, matching the pre-`preferential_formula` behavior of reusing the
  # whole fixed-effect design matrix.
  fit <- sdmTMB(
    catch ~ region,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(
      backend = "tmb",
      preferential_grid = pref_grid,
      preferential_response = "sampled",
      preferential_formula = ~region
    ),
    do_fit = FALSE
  )
  pref <- fit$tmb_data$preferential
  expect_identical(pref$n_pref, nrow(pref_grid))
  expect_length(pref$R_i, nrow(pref_grid))
  expect_true(all(pref$R_i %in% c(0, 1)))
  expect_identical(colnames(pref$X_pref_ij), colnames(fit$tmb_data$X_ij[[1]]))
  expect_identical(nrow(pref$A_pref), nrow(pref_grid))
  expect_identical(ncol(pref$A_pref), mesh$mesh$n)
  expect_length(pref$year_i_pref, nrow(pref_grid))
  expect_identical(pref$b_pref_type, 0L) # "constant"
  # `preferential_formula = ~region` selected both columns, so neither is
  # zeroed out.
  expect_true(all(pref$X_pref_ij[, "(Intercept)"] == 1))
  expect_identical(unname(pref$X_pref_ij[, "regionB"]), as.numeric(pref_grid$region == "B"))

  # default `preferential_formula` (`NULL`, treated as `~1`): only the
  # intercept column is reused, `regionB` is zeroed out even though
  # `pref_grid` still has a `region` column.
  fit_default <- sdmTMB(
    catch ~ region,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(
      backend = "tmb",
      preferential_grid = pref_grid,
      preferential_response = "sampled"
    ),
    do_fit = FALSE
  )
  pref_default <- fit_default$tmb_data$preferential
  expect_identical(colnames(pref_default$X_pref_ij), colnames(fit_default$tmb_data$X_ij[[1]]))
  expect_true(all(pref_default$X_pref_ij[, "(Intercept)"] == 1))
  expect_true(all(pref_default$X_pref_ij[, "regionB"] == 0))

  # default `preferential_formula` doesn't need `region` in the grid at all
  pref_grid_no_region <- pref_grid
  pref_grid_no_region$region <- NULL
  fit_default2 <- sdmTMB(
    catch ~ region,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(
      backend = "tmb",
      preferential_grid = pref_grid_no_region,
      preferential_response = "sampled"
    ),
    do_fit = FALSE
  )
  expect_true(all(fit_default2$tmb_data$preferential$X_pref_ij[, "(Intercept)"] == 1))

  # preferential_grid defaults to `data` when NULL
  dat2 <- dat
  dat2$sampled <- rbinom(n, 1, 0.5)
  fit2 <- sdmTMB(
    catch ~ region,
    data = dat2, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(backend = "tmb", preferential_response = "sampled"),
    do_fit = FALSE
  )
  expect_identical(fit2$tmb_data$preferential$n_pref, nrow(dat2))

  # feature off -> harmless zero-length placeholder
  fit3 <- sdmTMB(
    catch ~ region,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    do_fit = FALSE
  )
  expect_identical(fit3$tmb_data$preferential$n_pref, 0L)
  expect_length(fit3$tmb_data$preferential$R_i, 0L)

  # hard error: novel factor level in a factor preferential_formula uses
  bad_grid <- pref_grid
  bad_grid$region <- factor(sample(c("A", "B", "C"), nrow(bad_grid), replace = TRUE))
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb",
        preferential_grid = bad_grid,
        preferential_response = "sampled",
        preferential_formula = ~region
      ),
      do_fit = FALSE
    ),
    regexp = "new levels"
  )

  # hard error: preferential_grid missing a column required by an
  # explicit `preferential_formula` (the default `~1` needs no covariate
  # columns at all, see above, so this only errors once `region` is
  # actually requested)
  bad_grid2 <- pref_grid
  bad_grid2$region <- NULL
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb",
        preferential_grid = bad_grid2,
        preferential_response = "sampled",
        preferential_formula = ~region
      ),
      do_fit = FALSE
    ),
    regexp = "Missing.*region"
  )

  # hard error: preferential_formula must be one-sided
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb",
        preferential_grid = pref_grid,
        preferential_response = "sampled",
        preferential_formula = catch ~ region
      ),
      do_fit = FALSE
    ),
    regexp = "one-sided"
  )

  # hard error: preferential_formula references a term not in `formula`'s
  # own fixed-effect design matrix -- the sub-model reuses `b_j`, it can't
  # invent a coefficient for a term the main model never fit
  bad_grid5 <- pref_grid
  bad_grid5$not_in_formula <- runif(nrow(bad_grid5))
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb",
        preferential_grid = bad_grid5,
        preferential_response = "sampled",
        preferential_formula = ~not_in_formula
      ),
      do_fit = FALSE
    ),
    regexp = "not present"
  )

  # hard error: preferential_grid missing a time slice
  bad_grid3 <- pref_grid[pref_grid$year != 2020, ]
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb",
        preferential_grid = bad_grid3,
        preferential_response = "sampled"
      ),
      do_fit = FALSE
    ),
    regexp = "time slice"
  )
})

test_that("preferential sampling fits, converges, and supports tidy()/simulate()", {
  skip_on_cran()
  set.seed(42)

  # 40 "years" is unrealistic; it's just replication so xi_s is identifiable
  # rather than sitting at a zero-variance boundary (see
  # scratch/tests/test-preferential-xi-diagnostic.R).
  years <- 1:40
  n <- 1500
  dat <- data.frame(
    x = runif(n, 0, 10),
    y = runif(n, 0, 10),
    year = rep(years, length.out = n)
  )
  dat$catch <- rpois(n, exp(1 + 0.15 * dat$y))
  mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 2.5)

  # 7x7 grid, finer than the fitting mesh's cutoff, for spatial replication.
  grid_xy <- expand.grid(x = seq(0.5, 9.5, length.out = 7), y = seq(0.5, 9.5, length.out = 7))

  # Simulate a genuine spatial field (standing in for `xi_s`) via sdmTMB's own
  # simulate_new(), so the sampling submodel has real residual signal rather
  # than pure noise. simulate_new() requires the mesh built from the exact
  # data it simulates over, hence a separate small mesh here (not `mesh`,
  # which the fit below uses).
  mesh_grid <- make_mesh(grid_xy, xy_cols = c("x", "y"), cutoff = 1.4)
  sim_xi <- simulate_new(
    formula = ~1, data = grid_xy, mesh = mesh_grid, family = gaussian(),
    range = 3, sigma_O = 0.8, phi = 0.01, B = 0, seed = 99,
    control = sdmTMBcontrol(backend = "tmb")
  )
  xi_true <- sim_xi$omega_s

  pref_grid <- do.call(rbind, lapply(years, function(yr) {
    g <- grid_xy
    g$year <- yr
    g$sampled <- rbinom(nrow(g), 1, plogis(-1 + 0.3 * g$y + xi_true))
    g
  }))

  # Legacy C++ assembly smoke test only: catch and sampling are simulated
  # independently, so this is not evidence of bias correction.
  fit <- sdmTMB(
    catch ~ 1,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(
      backend = "tmb",
      newton_loops = 0,
      preferential_grid = pref_grid,
      preferential_response = "sampled"
    )
  )
  expect_lt(max(abs(fit$gradients)), 1e-2)
  expect_identical(fit$model$convergence, 0L)
  expect_true(fit$sd_report$pdHess)
  expect_true(all(c("gamma_0", "b_pref", "range_xi", "sigma_xi") %in%
    names(fit$sd_report$value)))

  td <- tidy(fit, "ran_pars")
  expect_true(all(c("gamma_0", "b_pref", "range_xi", "sigma_xi") %in% td$term))

  sims <- simulate(fit, nsim = 3, seed = 1)
  expect_true(is.matrix(sims))
  expect_identical(dim(sims), c(nrow(dat), 3L))
})

test_that("preferential sampling rejects unsupported backends and features", {
  skip_on_cran()
  set.seed(1)
  dat <- data.frame(
    x = runif(60, 0, 10), y = runif(60, 0, 10),
    year = rep(2018:2020, length.out = 60), depth = runif(60)
  )
  dat$catch <- rpois(60, 5)
  dat$sampled <- rbinom(60, 1, 0.5)
  mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 2)
  build <- function(formula = catch ~ 1, family = poisson(), backend = "tmb",
                    b_type = "constant", ...) {
    sdmTMB(
      formula,
      data = dat, mesh = mesh, time = "year", family = family,
      control = sdmTMBcontrol(
        backend = backend,
        preferential_response = "sampled",
        preferential_b_type = b_type
      ),
      do_fit = FALSE, ...
    )
  }
  expect_error(build(backend = "rtmb"), "not yet implemented for the RTMB")
  expect_error(build(family = delta_gamma()), "delta models")
  expect_error(build(b_type = "rw"), "other than")
  expect_error(build(b_type = "iid"), "other than")
  expect_error(build(catch ~ s(depth)), "smoothers")
  expect_error(build(catch ~ breakpt(depth)), "threshold")
  expect_error(build(time_varying = ~depth), "time_varying")
  expect_error(build(spatial_varying = ~depth), "spatial_varying")
  expect_error(build(anisotropy = TRUE), "anisotropy")
  expect_error(
    sdmTMB(
      catch ~ 1,
      data = dat, mesh = mesh, time = "year", family = poisson(),
      control = sdmTMBcontrol(
        backend = "tmb", normalize = TRUE,
        preferential_response = "sampled"
      ),
      do_fit = FALSE
    ),
    "normalize"
  )

  # The RTMB objective must also refuse preferential data from internal
  # callers that bypass sdmTMB(), rather than omitting the likelihood.
  fit <- build()
  expect_error(
    make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
      fit$tmb_random, backend = "rtmb"),
    "not yet implemented for the RTMB"
  )
})

test_that("preferential shared design reuses the fitted poly() basis", {
  # Regression case for the legacy builder, which rebuilds terms from
  # `preferential_formula` and so recomputes poly() on the grid values.
  # Enable once the shared design is built from the fitted terms (phase 1).
  skip("Legacy preferential design ignores fitted predvars")
  dat <- data.frame(
    x = seq(0, 9), y = seq(0, 9), year = 1L, depth = 1:10,
    catch = rpois(10, 5)
  )
  grid <- data.frame(x = 0, y = 0, year = 1L, depth = c(2, 4, 7, 9),
    sampled = c(0, 1, 0, 1))
  mesh <- make_mesh(dat, xy_cols = c("x", "y"), n_knots = 4, type = "kmeans")
  fit <- sdmTMB(
    catch ~ poly(depth, 2),
    data = dat, mesh = mesh, time = "year", family = poisson(),
    spatial = "off",
    control = sdmTMBcontrol(
      backend = "tmb",
      preferential_grid = grid,
      preferential_response = "sampled",
      preferential_formula = ~ poly(depth, 2)
    ),
    do_fit = FALSE
  )
  # The fitted terms carry poly()'s basis in `predvars`.
  fitted_terms <- stats::terms(stats::model.frame(catch ~ poly(depth, 2), dat))
  expected <- stats::model.matrix(stats::delete.response(fitted_terms), grid)
  expect_equal(fit$tmb_data$preferential$X_pref_ij, expected,
    ignore_attr = TRUE)
})
