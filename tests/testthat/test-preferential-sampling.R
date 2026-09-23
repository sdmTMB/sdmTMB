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

  # happy path
  fit <- sdmTMB(
    catch ~ region,
    data = dat, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(
      preferential_grid = pref_grid,
      preferential_response = "sampled",
      preferential_b_type = "rw"
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
  expect_identical(pref$b_pref_type, 2L) # "rw"

  # preferential_grid defaults to `data` when NULL
  dat2 <- dat
  dat2$sampled <- rbinom(n, 1, 0.5)
  fit2 <- sdmTMB(
    catch ~ region,
    data = dat2, mesh = mesh, time = "year",
    family = poisson(),
    control = sdmTMBcontrol(preferential_response = "sampled"),
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

  # hard error: novel factor level in preferential_grid
  bad_grid <- pref_grid
  bad_grid$region <- factor(sample(c("A", "B", "C"), nrow(bad_grid), replace = TRUE))
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        preferential_grid = bad_grid,
        preferential_response = "sampled"
      ),
      do_fit = FALSE
    ),
    regexp = "new levels"
  )

  # hard error: preferential_grid missing a fixed-effect column
  bad_grid2 <- pref_grid
  bad_grid2$region <- NULL
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        preferential_grid = bad_grid2,
        preferential_response = "sampled"
      ),
      do_fit = FALSE
    ),
    regexp = "Missing.*region"
  )

  # hard error: preferential_grid missing a time slice
  bad_grid3 <- pref_grid[pref_grid$year != 2020, ]
  expect_error(
    sdmTMB(
      catch ~ region,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
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
    range = 3, sigma_O = 0.8, phi = 0.01, B = 0, seed = 99
  )
  xi_true <- sim_xi$omega_s

  pref_grid <- do.call(rbind, lapply(years, function(yr) {
    g <- grid_xy
    g$year <- yr
    g$sampled <- rbinom(nrow(g), 1, plogis(-1 + 0.3 * g$y + xi_true))
    g
  }))

  fit_one <- function(preferential_b_type) {
    sdmTMB(
      catch ~ 1,
      data = dat, mesh = mesh, time = "year",
      family = poisson(),
      control = sdmTMBcontrol(
        newton_loops = 0,
        preferential_grid = pref_grid,
        preferential_response = "sampled",
        preferential_b_type = preferential_b_type
      ),
      do_fit = TRUE
    )
  }

  for (b_type in c("constant", "rw", "iid")) {
    fit <- suppressWarnings(fit_one(b_type)) # rw/iid: see NaN-SE note below
    expect_lt(max(abs(fit$gradients)), 1e-2)
    nms <- names(fit$sd_report$value)
    expect_true(all(c("gamma_0", "b_pref", "range_xi", "sigma_xi") %in% nms))
    if (b_type == "constant") {
      # no temporal-variance parameter for b_pref here, so this reliably
      # reaches a genuine interior optimum given xi_s's real signal above.
      expect_identical(fit$model$convergence, 0L)
      expect_true(fit$sd_report$pdHess)
    } else {
      # rw/iid add a b_pref temporal-variance parameter with no true
      # year-to-year signal to detect in this data, so it can sit at a
      # zero-variance boundary (same phenomenon as xi_s above, just for a
      # different parameter -- see scratch/tests/test-preferential-rw-
      # diagnostic*.R). convergence/pdHess aren't asserted here as a result.
      expect_true(fit$model$convergence %in% c(0L, 1L))
    }
  }

  fit_rw <- suppressWarnings(fit_one("rw"))
  fit_constant <- fit_one("constant")

  # fit_rw's b_pref variance sits at the ~0 boundary (see above), producing
  # harmless NaN-SE warnings from TMB's delta method.
  td <- suppressWarnings(tidy(fit_rw, "ran_pars"))
  expect_true(any(grepl("^gamma_0$", td$term)))
  expect_true(any(grepl("^b_pref:", td$term))) # rw -> one row per time slice
  expect_true(any(grepl("^range_xi$", td$term)))
  expect_true(any(grepl("^sigma_xi$", td$term)))

  td_const <- tidy(fit_constant, "ran_pars")
  expect_true("b_pref" %in% td_const$term) # constant -> single row, no time suffix

  sims <- simulate(fit_rw, nsim = 3, seed = 1)
  expect_true(is.matrix(sims))
  expect_identical(ncol(sims), 3L)
  expect_identical(nrow(sims), nrow(dat))
})
