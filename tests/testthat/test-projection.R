test_that("projection helpers preserve parameter order", {
  x <- setNames(1:10, c("a", "a", rep("epsilon_st", 6), "z", "z"))
  expect_equal(
    unname(insert_pars(x, "epsilon_st", 2)),
    c(1:8, 0, 0, 9:10)
  )
  expect_equal(
    unname(insert_pars(x, "epsilon_st", 4, n_groups = 2)),
    c(1:2, 3:5, 0, 0, 6:8, 0, 0, 9:10)
  )

  at_start <- setNames(1:4, c("epsilon_st", "epsilon_st", "z", "z"))
  at_end <- setNames(1:4, c("a", "a", "epsilon_st", "epsilon_st"))
  expect_equal(unname(insert_pars(at_start, "epsilon_st", 2)), c(1:2, 0, 0, 3:4))
  expect_equal(unname(insert_pars(at_end, "epsilon_st", 2)), c(1:4, 0, 0))
  expect_equal(
    unname(insert_pars(x, "epsilon_st", 3, n_groups = 3)),
    c(1:2, 3:4, 0, 5:6, 0, 7:8, 0, 9:10)
  )
  expect_error(insert_pars(x, "epsilon_st", 1, n_groups = 2), "divide evenly")
  expect_error(insert_pars(x, "epsilon_st", 2, n_groups = 0), "positive")
  expect_error(insert_pars(x, "epsilon_st", 1.5), "whole number")
  expect_error(insert_pars(x, "missing", 2), "not found")
})

test_that("project supports all independent sampling combinations", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  historical_years <- sort(unique(pcod_2011$year))
  grid <- replicate_df(
    pcod_2011[seq(1, nrow(pcod_2011), by = 100), , drop = FALSE],
    "year", c(historical_years, max(historical_years) + 2)
  )
  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011,
    mesh = mesh,
    time = "year",
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "ar1",
    family = tweedie()
  )

  controls <- expand.grid(
    sample_fe = c(FALSE, TRUE),
    sample_historical_re = c(FALSE, TRUE),
    sample_future_re = c(FALSE, TRUE)
  )
  reports <- lapply(seq_len(nrow(controls)), function(i) {
    project(
      fit,
      grid,
      nsim = 2,
      sample_fe = controls$sample_fe[[i]],
      sample_historical_re = controls$sample_historical_re[[i]],
      sample_future_re = controls$sample_future_re[[i]],
      return_tmb_report = TRUE,
      silent = TRUE
    )
  })
  expect_length(reports, 8L)
  expect_true(all(vapply(reports, length, integer(1L)) == 2L))

  n_historical <- nrow(fit$time_lu)
  for (i in which(!controls$sample_future_re)) {
    for (report in reports[[i]]) {
      rho <- report$rho[[1L]]
      expect_equal(
        report$epsilon_st[, n_historical + 1L, 1L],
        rho * report$epsilon_st[, n_historical, 1L]
      )
    }
  }
  for (i in which(!controls$sample_historical_re)) {
    expect_equal(
      reports[[i]][[1]]$epsilon_st[, seq_len(n_historical), 1L],
      reports[[i]][[2]]$epsilon_st[, seq_len(n_historical), 1L]
    )
  }

  stochastic <- which(
    !controls$sample_fe & !controls$sample_historical_re &
      controls$sample_future_re
  )
  future <- grid$year > max(historical_years)
  expect_false(all(
    reports[[stochastic]][[1]]$proj_epsilon_st_A_vec[future, 1L] == 0
  ))
  expect_false(isTRUE(all.equal(
    reports[[stochastic]][[1]]$proj_epsilon_st_A_vec[future, 1L],
    reports[[stochastic]][[2]]$proj_epsilon_st_A_vec[future, 1L]
  )))
  expect_false(isTRUE(all.equal(
    reports[[stochastic]][[1]]$proj_eta[future, 1L],
    reports[[stochastic]][[2]]$proj_eta[future, 1L]
  )))
})

test_that("project treats REML regression coefficients as parameters", {
  skip_on_cran()

  data <- pcod_2011
  data$year_f <- factor(data$year)
  mesh <- make_mesh(data, c("X", "Y"), cutoff = 20)
  historical_years <- sort(unique(data$year))
  future_year <- max(historical_years) + 2
  grid <- replicate_df(
    data[seq(1, nrow(data), by = 100), , drop = FALSE],
    "year", c(historical_years, future_year)
  )
  grid$year_f <- factor(grid$year, levels = c(historical_years, future_year))
  fit <- sdmTMB(
    log1p(density) ~ 1 + (1 | year_f),
    data = data,
    mesh = mesh,
    time = "year",
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    reml = TRUE,
    family = gaussian()
  )

  parameter_only <- project(
    fit, grid, nsim = 2, sample_fe = TRUE,
    sample_historical_re = FALSE, sample_future_re = FALSE,
    return_tmb_report = TRUE, allow_new_levels = TRUE, silent = TRUE
  )
  historical_only <- project(
    fit, grid, nsim = 2, sample_fe = FALSE,
    sample_historical_re = TRUE, sample_future_re = FALSE,
    return_tmb_report = TRUE, allow_new_levels = TRUE, silent = TRUE
  )

  expect_false(isTRUE(all.equal(
    parameter_only[[1]]$eta_fixed_i, parameter_only[[2]]$eta_fixed_i
  )))
  expect_equal(
    historical_only[[1]]$eta_fixed_i, historical_only[[2]]$eta_fixed_i
  )
  expect_equal(parameter_only[[1]]$re_b_pars, parameter_only[[2]]$re_b_pars)
  expect_false(isTRUE(all.equal(
    historical_only[[1]]$re_b_pars, historical_only[[2]]$re_b_pars
  )))
})

test_that("project propagates conditional means for future process effects", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  historical_years <- sort(unique(pcod_2011$year))
  grid <- replicate_df(
    pcod_2011[seq(1, nrow(pcod_2011), by = 100), , drop = FALSE],
    "year", c(historical_years, max(historical_years) + 2)
  )
  for (spatiotemporal in c("ar1", "rw", "iid")) {
    fit <- sdmTMB(
      density ~ 1,
      data = pcod_2011,
      mesh = mesh,
      time = "year",
      extra_time = historical_years,
      spatial = "off",
      spatiotemporal = spatiotemporal,
      family = tweedie()
    )
    report <- project(
      fit, grid, nsim = 1, sample_fe = FALSE,
      sample_historical_re = FALSE, sample_future_re = FALSE,
      return_tmb_report = TRUE, silent = TRUE
    )[[1L]]
    n_historical <- nrow(fit$time_lu)
    if (spatiotemporal == "ar1") {
      rho <- report$rho[[1L]]
      expect_equal(
        report$epsilon_st[, n_historical + 1L, 1L],
        rho * report$epsilon_st[, n_historical, 1L]
      )
    } else if (spatiotemporal == "rw") {
      expect_equal(
        report$epsilon_st[, n_historical + 1L, 1L],
        report$epsilon_st[, n_historical, 1L]
      )
    } else {
      expect_equal(
        as.numeric(report$epsilon_st[, n_historical + 1L, 1L]),
        rep(0, nrow(report$epsilon_st))
      )
    }
  }
})

test_that("project can zero or fix future spatiotemporal effects", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  historical_years <- sort(unique(pcod_2011$year))
  grid <- replicate_df(
    pcod_2011[seq(1, nrow(pcod_2011), by = 100), , drop = FALSE],
    "year", c(historical_years, max(historical_years) + 2)
  )
  fit <- sdmTMB(
    density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
    extra_time = historical_years, spatial = "off",
    spatiotemporal = "ar1", family = tweedie()
  )
  n_historical <- nrow(fit$time_lu)

  zero <- project(
    fit, grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, future_re = "zero",
    return_tmb_report = TRUE, silent = TRUE
  )[[1L]]
  future_index <- seq.int(n_historical + 1L, dim(zero$epsilon_st)[2L])
  expect_true(all(zero$epsilon_st[, future_index, 1L] == 0))

  fixed <- project(
    fit, grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, future_re = "fix",
    return_tmb_report = TRUE, silent = TRUE
  )[[1L]]
  for (tt in future_index) {
    expect_equal(
      fixed$epsilon_st[, tt, 1L],
      fixed$epsilon_st[, n_historical, 1L]
    )
  }
})

test_that("projection time extensions use the fitted cadence", {
  annual <- data.frame(time_from_data = 2020:2022)
  expect_equal(
    project_time_extension(annual, c(2020:2023, 2025))$time_from_data,
    2023:2025
  )

  five_year <- data.frame(time_from_data = seq(2000, 2010, by = 5))
  expect_equal(project_time_extension(five_year, 2015)$time_from_data, 2015)

  one_time <- data.frame(time_from_data = 2010)
  expect_equal(
    project_time_extension(one_time, c(2015, 2020))$time_from_data,
    c(2015, 2020)
  )
  expect_error(project_time_extension(one_time, 2015), "at least two")
  expect_error(project_time_extension(one_time, c(2015, 2025)), "align")

  decimal <- project_time_extension(
    data.frame(time_from_data = c(0, 0.1, 0.2)),
    c(0, 0.1 + 0.2)
  )
  expect_equal(decimal$time_from_data, 0.3, tolerance = 1e-12)
  expect_equal(decimal$canonical_time, c(0, 0.3), tolerance = 1e-12)

  expect_error(
    project_time_extension(data.frame(time_from_data = c(2000, 2001, 2003)), 2004),
    "constant cadence"
  )
  expect_error(project_time_extension(five_year, 2012), "align")
  expect_error(project_time_extension(annual, 2021), "No new time")
  expect_error(project_time_extension(annual, c(2021.5, 2023)), "neither fitted nor future")
  expect_error(project_time_extension(annual, NA_real_), "finite")
  expect_error(project_time_extension(annual, Inf), "finite")
  expect_error(project_time_extension(annual, "2023"), "numeric")
})

test_that("project simulation report extraction has stable dimensions", {
  reports <- list(
    list(vector = 1:3, matrix = matrix(1:6, 3, 2)),
    list(vector = 4:6, matrix = matrix(7:12, 3, 2))
  )
  expect_equal(dim(extract_project_sims(reports, "vector")), c(3L, 2L))
  expect_equal(dim(extract_project_sims(reports, "matrix")), c(3L, 2L, 2L))
  expect_equal(extract_project_sims(reports, "vector")[, 2], 4:6)
  delta_reports <- list(
    list(delta = matrix(1:6, nrow = 3, ncol = 2)),
    list(delta = matrix(7:12, nrow = 3, ncol = 2))
  )
  expect_equal(dim(extract_project_sims(delta_reports, "delta")), c(3L, 2L, 2L))
  expect_error(
    extract_project_sims(
      c(delta_reports, list(list(delta = matrix(1:4, nrow = 2, ncol = 2)))),
      "changed dimensions"
    )
  )
  expect_error(extract_project_sims(reports, "missing"), "not found")
})

test_that("project() handles decimal time cadence and preserves rows", {
  skip_on_cran()

  decimal_data <- do.call(
    rbind,
    lapply(split(pcod_2011, pcod_2011$year), head, 30)
  )
  decimal_data <- decimal_data[decimal_data$year < 2017, , drop = FALSE]
  decimal_data$time_decimal <- (decimal_data$year - 2011) / 20
  mesh <- make_mesh(decimal_data, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ 1,
    data = decimal_data,
    time = "time_decimal",
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = tweedie()
  )

  newdata <- decimal_data[seq(1, nrow(decimal_data), by = 3), , drop = FALSE]
  future <- newdata[seq_len(5), , drop = FALSE]
  future$time_decimal <- 0.1 + 0.2
  newdata <- rbind(newdata, future)
  newdata <- newdata[c(nrow(newdata), seq_len(nrow(newdata) - 1L)), , drop = FALSE]

  out <- project(
    fit, newdata, nsim = 2, sample_fe = FALSE,
    sample_historical_re = FALSE, silent = TRUE
  )
  expect_identical(names(out), "est")
  expect_equal(dim(out$est), c(nrow(newdata), 2L))
})

test_that("project() works with delta models", {
  skip_on_cran()

  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)

  # we could fit our model like this, but for long projections, this becomes slow:
  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = all_years, #< note that all years here
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )
  p <- predict(fit, newdata = proj_grid)

  # instead, we could fit our model like this and then take simulation draws
  # from the projection time period:
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )
  set.seed(1)
  out <- project(
    fit2, newdata = proj_grid, nsim = 25, sample_fe = FALSE,
    sample_historical_re = FALSE
  )
  none <- out
  expect_identical(names(out), c("est1", "est2", "epsilon_st1", "epsilon_st2"))

  proj_grid$est_mean1 <- apply(out$est1, 1, mean)
  proj_grid$est_mean2 <- apply(out$est2, 1, mean)
  proj_grid$eps_mean1 <- apply(out$epsilon_st1, 1, mean)
  proj_grid$eps_mean2 <- apply(out$epsilon_st2, 1, mean)

  i <- p$year == 2023

  OK_COR <- 0.98
  expect_gt(cor(p$est1[i], proj_grid$est_mean1[i]), OK_COR)
  expect_gt(cor(p$est2[i], proj_grid$est_mean2[i]), OK_COR)
  expect_gt(cor(p$epsilon_st1[i], proj_grid$eps_mean1[i]), OK_COR)
  expect_gt(cor(p$epsilon_st2[i], proj_grid$eps_mean2[i]), OK_COR)

  # Future simulation must not shift the fitted delta-component fields.
  historical <- proj_grid$year %in% historical_years
  p_historical <- predict(fit2, newdata = proj_grid[historical, ])
  expect_equal(out$epsilon_st1[historical, 1], p_historical$epsilon_st1)
  expect_equal(out$epsilon_st2[historical, 1], p_historical$epsilon_st2)

  # test return_tmb_report:

  # if instead we wanted to grab, say, the spatiotemporal random field values,
  # we can return the report and work with the raw output ourselves:
  set.seed(1)
  out <- project(
    fit2,
    newdata = proj_grid,
    nsim = 25,
    sample_fe = FALSE,
    sample_historical_re = FALSE,
    return_tmb_report = TRUE #< difference from above example
  )

  expect_true(all(vapply(out, \(x) "proj_eta" %in% names(x), logical(1L))))
  eps <- lapply(out, \(x) x[["proj_epsilon_st_A_vec"]][, 1])
  eps <- do.call(cbind, eps)
  eps_mean <- apply(eps, 1, mean)
  expect_gt(cor(eps_mean, proj_grid$eps_mean1), OK_COR)

  eps <- lapply(out, \(x) x[["proj_epsilon_st_A_vec"]][, 2])
  eps <- do.call(cbind, eps)
  eps_mean <- apply(eps, 1, mean)
  expect_gt(cor(eps_mean, proj_grid$eps_mean2), OK_COR)

  # test the types of uncertainty
  set.seed(1)
  both <- project(fit2, newdata = proj_grid, nsim = 20)
  set.seed(1)
  suppressWarnings({
    random <- project(
      fit2, newdata = proj_grid, nsim = 20,
      sample_fe = FALSE, sample_historical_re = TRUE
    )
  })

  sd_both <- mean(apply(both$est1, 1, sd)[i])
  sd_none <- mean(apply(none$est1, 1, sd)[i])
  sd_random <- mean(apply(random$est1, 1, sd)[i], na.rm = TRUE)
  # expect_gt(sd_both, sd_random) # !!?
  expect_gt(sd_both, sd_none)
  expect_gt(sd_random, sd_none)
})

test_that("project() pads delta spatiotemporal fields by active component", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  historical_years <- sort(unique(pcod_2011$year))
  grid <- replicate_df(
    pcod_2011[seq(1, nrow(pcod_2011), by = 50), , drop = FALSE],
    "year", c(historical_years, max(historical_years) + 2)
  )
  historical <- grid$year %in% historical_years

  for (spatiotemporal in list(list("ar1", "off"), list("off", "ar1"))) {
    fit <- sdmTMB(
      density ~ 1,
      data = pcod_2011,
      mesh = mesh,
      spatial = "off",
      time = "year",
      spatiotemporal = spatiotemporal,
      family = delta_gamma(),
      control = sdmTMBcontrol(newton_loops = 0)
    )
    out <- project(
      fit, grid, nsim = 1, sample_fe = FALSE,
      sample_historical_re = FALSE, silent = TRUE
    )
    fitted <- predict(fit, grid[historical, , drop = FALSE])
    expect_equal(out$est1[historical, 1], fitted$est1)
    expect_equal(out$est2[historical, 1], fitted$est2)
  }
})

test_that("project() works with non-delta models", {
  skip_on_cran()

  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)

  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = all_years, #< note that all years here
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = tweedie()
  )
  p <- predict(fit, newdata = proj_grid)
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off", # speed
    spatiotemporal = "ar1",
    data = dogfish,
    mesh = mesh,
    family = tweedie()
  )
  set.seed(1)
  out <- project(
    fit2, newdata = proj_grid, nsim = 25, sample_fe = FALSE,
    sample_historical_re = FALSE
  )
  expect_identical(names(out), c("est", "epsilon_st"))

  i <- p$year == 2023
  proj_grid$est_mean <- apply(out$est, 1, mean)
  proj_grid$eps_mean <- apply(out$epsilon_st, 1, mean)
  expect_gt(cor(p$est[i], proj_grid$est_mean[i]), 0.98)
  expect_gt(cor(p$epsilon_st[i], proj_grid$eps_mean[i]), 0.98)

  eta <- project(
    fit2, newdata = proj_grid, nsim = 2, sample_fe = FALSE,
    sample_historical_re = FALSE,
    sims_var = "proj_eta", silent = TRUE
  )
  expect_equal(dim(eta), c(nrow(proj_grid), 1L, 2L))

  old_fit <- fit2
  old_fit$tmb_data$sim_obs <- NULL
  old_fit_out <- project(
    old_fit, newdata = proj_grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, silent = TRUE
  )
  expect_equal(dim(old_fit_out$est), c(nrow(proj_grid), 1L))
  expect_error(
    project(
      fit2, newdata = proj_grid, nsim = 1, sample_fe = FALSE,
      sample_historical_re = FALSE,
      sims_var = "not_a_report_element", silent = TRUE
    ),
    "not found"
  )
})

test_that("project() works with time-varying effects", {
  skip_on_cran()
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)

  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    time_varying = ~1,
    time_varying_type = "ar1",
    extra_time = all_years, #< note that all years here
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = tweedie()
  )
  p <- predict(fit, newdata = proj_grid)
  fit2 <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    time_varying = ~1,
    time_varying_type = "ar1",
    extra_time = historical_years, #< does *not* include projection years
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = tweedie()
  )
  set.seed(1)
  out <- project(
    fit2, newdata = proj_grid, nsim = 25, sample_fe = FALSE,
    sample_historical_re = FALSE
  )
  expect_identical(names(out), c("est"))
  i <- p$year == 2023
  expect_equal(mean(p$est[i]), 5.983172, tolerance = 1e-3)

  dogfish$depth_scaled <- as.numeric(scale(dogfish$depth))
  proj_grid$depth_scaled <- (proj_grid$depth - mean(dogfish$depth)) / sd(dogfish$depth)
  fit3 <- sdmTMB(
    catch_weight ~ 0,
    time = "year",
    time_varying = ~ 0 + depth_scaled + I(depth_scaled^2),
    time_varying_type = "ar1",
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    family = tweedie()
  )
  raw <- project(
    fit3, newdata = proj_grid, nsim = 2, sample_fe = FALSE,
    sample_historical_re = FALSE,
    return_tmb_report = TRUE, silent = TRUE
  )
  expect_equal(dim(raw[[1]]$b_rw_t), c(length(all_years), 2L, 1L))
  expect_equal(
    raw[[1]]$b_rw_t[seq_along(historical_years), , , drop = FALSE],
    get_pars(fit3)$b_rw_t,
    tolerance = 1e-3
  )
  raw_mean <- project(
    fit3, newdata = proj_grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, sample_future_re = FALSE,
    return_tmb_report = TRUE, silent = TRUE
  )[[1L]]
  rho_time <- 2 * plogis(get_pars(fit3)$rho_time_unscaled[1, 1]) - 1
  n_historical <- nrow(fit3$time_lu)
  expect_equal(
    raw_mean$b_rw_t[n_historical + 1L, 1L, 1L],
    rho_time * raw_mean$b_rw_t[n_historical, 1L, 1L]
  )

  raw_zero <- project(
    fit3, newdata = proj_grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, future_re = "zero",
    return_tmb_report = TRUE, silent = TRUE
  )[[1L]]
  expect_equal(raw_zero$b_rw_t[n_historical + 1L, , 1L], c(0, 0))

  raw_fix <- project(
    fit3, newdata = proj_grid, nsim = 1, sample_fe = FALSE,
    sample_historical_re = FALSE, future_re = "fix",
    return_tmb_report = TRUE, silent = TRUE
  )[[1L]]
  expect_equal(
    raw_fix$b_rw_t[n_historical + 1L, , 1L],
    raw_fix$b_rw_t[n_historical, , 1L]
  )
})


test_that("project() works/fails as expected in some less obvious situations", {
  skip_on_cran()

  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 35)
  historical_years <- 2004:2022
  to_project <- 1
  future_years <- seq(max(historical_years) + 1, max(historical_years) + to_project)
  all_years <- c(historical_years, future_years)
  proj_grid <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)

  # no time model:
  fit <- sdmTMB(
    catch_weight ~ 1,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    mesh = mesh,
    family = delta_gamma()
  )

  set.seed(1)
  expect_message(
    out <- project(
      fit, newdata = proj_grid, nsim = 2, sample_fe = FALSE,
      sample_historical_re = FALSE
    ),
    "structures"
  )
  expect_equal(unique(as.numeric(out$est1)), 0.8032636, tolerance = 1e-3)
  expect_message(out <- project(fit, newdata = proj_grid, nsim = 1, silent = TRUE), "structures")
  expect_message(
    out <- project(
      fit, newdata = proj_grid, nsim = 1,
      sample_fe = FALSE, sample_historical_re = TRUE, silent = TRUE
    )
  )

  # no time argument specified errors:
  fit <- sdmTMB(
    catch_weight ~ 1,
    # time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    data = dogfish,
    # mesh = mesh,
    family = delta_gamma()
  )

  expect_error(out2 <- project(fit, newdata = proj_grid, nsim = 1), regexp = "time")

  # no future time errors:
  fit <- sdmTMB(
    catch_weight ~ 0,
    time = "year",
    offset = log(dogfish$area_swept),
    extra_time = historical_years,
    spatial = "off",
    spatiotemporal = "off",
    time_varying = ~ 1,
    data = dogfish,
    family = tweedie()
  )
  expect_error(out <- project(fit, newdata = subset(proj_grid, year %in% historical_years), nsim = 1), "new")

  # newdata is missing a time step; make sure that's fine and matches not missing the time step:
  all_years <- c(historical_years, 2023:2025)
  proj_grid <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)
  set.seed(1)
  out <- project(fit, newdata = proj_grid, nsim = 1)

  all_years <- c(historical_years, c(2023, 2025))
  proj_grid2 <- replicate_df(wcvi_grid[seq(1, nrow(wcvi_grid), by = 20), ], "year", all_years)
  set.seed(1)
  out2 <- project(fit, newdata = proj_grid2, nsim = 1)

  i <- proj_grid$year %in% c(2023, 2025)
  i2 <- proj_grid2$year %in% c(2023, 2025)
  expect_equal(out$est[i,], out2$est[i2,])

  expect_error(project(fit, proj_grid, nsim = 1.5), "whole number")
  expect_error(project(fit, proj_grid, nsim = Inf), "finite")
  expect_error(project(fit, proj_grid, nsim = .Machine$integer.max + 1), "integer.max")
  expect_error(project(fit, proj_grid, sample_fe = NA), "sample_fe")
  expect_error(project(fit, proj_grid, sample_historical_re = c(TRUE, FALSE)), "sample_historical_re")
  expect_error(project(fit, proj_grid, sample_future_re = NA), "sample_future_re")
  expect_error(project(fit, proj_grid, future_re = "unknown"), "arg")
  expect_error(project(fit, proj_grid, silent = NA), "silent")
  expect_error(project(fit, proj_grid, return_tmb_report = c(TRUE, FALSE)), "return_tmb_report")
  expect_error(project(fit, proj_grid, uncertainty = "none"), "uncertainty")
  expect_error(project(fit, proj_grid, sim_re = c(FALSE, TRUE)), "sim_re")
  expect_error(project(fit, proj_grid, simulate_re = c(FALSE, TRUE)), "no longer supported")
  expect_error(project(fit, proj_grid, sims_var = character()), "sims_var")
  expect_error(project(fit, transform(proj_grid, year = NULL)), "missing")
})
