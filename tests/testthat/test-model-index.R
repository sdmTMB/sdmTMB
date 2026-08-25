test_that("model_index parser removes one direct special", {
  p <- .parse_model_index_svc(~ 0 + depth + model_index(factor(year)))
  expect_true(p$active)
  expect_identical(p$term, "factor(year)")
  expect_identical(p$label, "model_index(factor(year))")
  expect_identical(attr(terms(p$ordinary_formula), "term.labels"), "depth")
  expect_identical(attr(terms(p$ordinary_formula), "intercept"), 0L)

  expect_error(model_index(year), "only be used inside")
  expect_error(.parse_model_index_svc(~ model_index()), "exactly one")
  expect_error(.parse_model_index_svc(~ model_index(year, depth)), "exactly one")
  expect_error(
    .parse_model_index_svc(~ model_index(year) + model_index(depth)),
    "At most one"
  )
  expect_error(.parse_model_index_svc(~ depth:model_index(year)), "one direct term")
})

test_that("model_index creates one ordered SVC placeholder and metadata", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ depth_scaled + year_factor,
    spatial_varying = ~ 0 + depth_scaled + model_index(year_factor),
    data = d,
    mesh = mesh,
    time = "year",
    family = tweedie(link = "log"),
    do_fit = FALSE
  )

  expect_identical(
    colnames(fit$tmb_data$z_i),
    c("depth_scaled", "model_index(year_factor)")
  )
  expect_true(all(fit$tmb_data$z_i[, 2L] == 0))
  expect_identical(fit$tmb_data$model_index_z, 1L)
  expect_identical(fit$model_index$z_index, 2L)

  X <- fit$tmb_data$X_ij[[1L]]
  selected <- fit$model_index$coefficient_indices[[1L]]
  expect_identical(selected, which(attr(X, "assign") %in% c(0L, 2L)))
  expect_equal(nrow(fit$tmb_data$model_index_X_t[[1L]]), length(unique(d$year)))
  expect_equal(ncol(fit$tmb_data$model_index_X_t[[1L]]), ncol(X))
  expect_true(all(fit$tmb_data$model_index_X_t[[1L]][, -selected, drop = FALSE] == 0))
})

test_that("model_index time design handles cell means, row order, and delta formulas", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  args <- list(
    formula = density ~ 0 + year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    mesh = mesh,
    time = "year",
    family = delta_gamma(),
    do_fit = FALSE
  )
  fit1 <- do.call(sdmTMB, c(args, list(data = d)))
  set.seed(1)
  fit2 <- do.call(sdmTMB, c(args, list(data = d[sample(nrow(d)), ])))

  expect_length(fit1$tmb_data$model_index_X_t, 2L)
  expect_identical(ncol(fit1$tmb_data$z_i), 1L)
  expect_false(any(attr(fit1$tmb_data$X_ij[[1L]], "assign")[
    fit1$model_index$coefficient_indices[[1L]]
  ] == 0L))
  expect_equal(fit1$tmb_data$model_index_X_t, fit2$tmb_data$model_index_X_t)

  formulas <- list(
    density ~ year_factor + depth_scaled,
    density ~ depth_scaled + year_factor
  )
  fit3 <- sdmTMB(
    formulas,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d,
    mesh = mesh,
    time = "year",
    family = delta_gamma(),
    do_fit = FALSE
  )
  expect_length(fit3$model_index$coefficient_indices, 2L)
  expect_true(all(vapply(fit3$model_index$coefficient_indices, length, integer(1L)) > 1L))
})

test_that("model_index time design is padded after fixed-effect columns are appended", {
  d <- data.frame(year = rep(1:3, each = 2), f = factor(rep(1:3, each = 2)))
  X <- model.matrix(~ f, d)
  parsed <- .parse_model_index_svc(~ 0 + model_index(f))
  idx <- .model_index_coefficient_indices(parsed, list(terms(~ f)), list(X))
  X_padded <- .append_nonlocal_coef_columns(X, "diffusion(x)")
  X_t <- .build_model_index_X_t(list(X_padded), idx, rep(0:2, each = 2), 3L, "f")

  expect_identical(dim(X_t[[1L]]), c(3L, ncol(X_padded)))
  expect_true(all(X_t[[1L]][, ncol(X_padded)] == 0))
})

test_that("model_index validates version-one constraints", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  base <- list(
    formula = density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d,
    mesh = mesh,
    family = tweedie(link = "log"),
    do_fit = FALSE
  )
  expect_error(do.call(sdmTMB, base), "`time` must be supplied")
  expect_error(do.call(sdmTMB, c(base, list(time = "year", extra_time = 2020))), "`extra_time`")
  expect_error(
    sdmTMB(
      density ~ year_factor,
      spatial_varying = ~ 0 + model_index(year_factor),
      data = d, mesh = mesh, time = "year",
      family = gaussian(), do_fit = FALSE
    ),
    "not guaranteed to be positive"
  )

  bad <- d
  bad$within_time <- seq_len(nrow(bad))
  expect_error(
    sdmTMB(
      density ~ within_time,
      spatial_varying = ~ 0 + model_index(within_time),
      data = bad, mesh = mesh, time = "year",
      family = tweedie(link = "log"), do_fit = FALSE
    ),
    "constant within"
  )

  one_time <- d[d$year == d$year[[1L]], ]
  one_time$year_factor <- droplevels(one_time$year_factor)
  one_time_mesh <- make_mesh(one_time, c("X", "Y"), cutoff = 20)
  expect_error(
    sdmTMB(
      density ~ year_factor,
      spatial_varying = ~ 0 + model_index(year_factor),
      data = one_time, mesh = one_time_mesh, time = "year",
      family = tweedie(link = "log"), do_fit = FALSE
    ),
    "at least two fitted time values"
  )

  expect_error(
    sdmTMB(
      density ~ depth_scaled,
      spatial_varying = ~ 0 + model_index(year_factor),
      data = d, mesh = mesh, time = "year",
      family = tweedie(link = "log"), do_fit = FALSE
    ),
    "must occur as a complete term"
  )

  expect_error(
    sdmTMB(
      list(density ~ year_factor, density ~ depth_scaled),
      spatial_varying = ~ 0 + model_index(year_factor),
      data = d, mesh = mesh, time = "year",
      family = delta_gamma(), do_fit = FALSE
    ),
    "every model formula"
  )
})

test_that("one-part TMB model index uses fixed effects and mean epsilon_st", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = tweedie(link = "log"), do_fit = FALSE
  )

  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 4
  epsilon_mean <- seq(-0.4, 0.4, length.out = fit$tmb_data$n_t)
  vertex_deviation <- seq(-0.2, 0.2, length.out = dim(parameters$epsilon_st)[1L])
  vertex_deviation <- vertex_deviation - mean(vertex_deviation)
  parameters$epsilon_st[, , 1L] <- outer(vertex_deviation, epsilon_mean, `+`)
  parameters$zeta_s[] <- 0

  obj <- TMB::MakeADFun(
    data = fit$tmb_data, parameters = parameters, map = fit$tmb_map,
    random = NULL, DLL = "sdmTMB", silent = TRUE
  )
  report <- obj$report(obj$par)
  uncentered <- drop(fit$tmb_data$model_index_X_t[[1L]] %*% parameters$b_j) +
    apply(parameters$epsilon_st[, , 1L], 2L, mean)
  expected <- uncentered - mean(uncentered)

  expect_equal(report$model_index_t, expected, tolerance = 1e-10)
  expect_equal(sum(report$model_index_t), 0, tolerance = 1e-10)
})

test_that("one-part TMB SVC multiplies the model-derived index", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    present ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = binomial(link = "logit"), do_fit = FALSE
  )

  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 5
  parameters$epsilon_st[] <- 0
  parameters$zeta_s[] <- 0
  make_report <- function(parameters) {
    obj <- TMB::MakeADFun(
      data = fit$tmb_data, parameters = parameters, map = fit$tmb_map,
      random = NULL, DLL = "sdmTMB", silent = TRUE
    )
    obj$report(obj$par)
  }
  report_zero <- make_report(parameters)
  parameters$zeta_s[] <- 0.7
  report_svc <- make_report(parameters)

  component_eta <- drop(fit$tmb_data$model_index_X_t[[1L]] %*% parameters$b_j)
  q <- log(plogis(component_eta))
  expected_index <- q - mean(q)
  expected_difference <- 0.7 * expected_index[fit$tmb_data$year_i + 1L]

  expect_equal(report_svc$model_index_t, expected_index, tolerance = 1e-10)
  expect_equal(
    drop(report_svc$eta_i - report_zero$eta_i),
    expected_difference,
    tolerance = 1e-9
  )
})

test_that("standard delta model index combines components before centering", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = delta_gamma(), do_fit = FALSE
  )
  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 5 - 1
  parameters$b_j2[] <- rev(seq_along(parameters$b_j2)) / 4
  parameters$epsilon_st[] <- 0
  parameters$zeta_s[] <- 0
  make_report <- function(parameters) {
    obj <- TMB::MakeADFun(
      data = fit$tmb_data, parameters = parameters, map = fit$tmb_map,
      random = NULL, DLL = "sdmTMB", silent = TRUE
    )
    obj$report(obj$par)
  }
  report_zero <- make_report(parameters)

  a1 <- drop(fit$tmb_data$model_index_X_t[[1L]] %*% parameters$b_j)
  a2 <- drop(fit$tmb_data$model_index_X_t[[2L]] %*% parameters$b_j2)
  q <- log(plogis(a1) * exp(a2))
  expected <- q - mean(q)
  separately_centered <- log(plogis(a1 - mean(a1))) + a2 - mean(a2)
  separately_centered <- separately_centered - mean(separately_centered)

  expect_equal(report_zero$model_index_t, expected, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(expected, separately_centered, tolerance = 1e-6)))

  parameters$zeta_s[, , 1L] <- 0.4
  parameters$zeta_s[, , 2L] <- -0.2
  report_svc <- make_report(parameters)
  expected_by_observation <- expected[fit$tmb_data$year_i + 1L]
  expect_equal(
    report_svc$eta_i[, 1L] - report_zero$eta_i[, 1L],
    0.4 * expected_by_observation,
    tolerance = 1e-9
  )
  expect_equal(
    report_svc$eta_i[, 2L] - report_zero$eta_i[, 2L],
    -0.2 * expected_by_observation,
    tolerance = 1e-9
  )
})

test_that("Poisson-link delta model index is the centered component sum", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = delta_gamma(type = "poisson-link"), do_fit = FALSE
  )
  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 6
  parameters$b_j2[] <- rev(seq_along(parameters$b_j2)) / 7
  parameters$epsilon_st[, , 1L] <- 0.15
  parameters$epsilon_st[, , 2L] <- -0.05

  obj <- TMB::MakeADFun(
    data = fit$tmb_data, parameters = parameters, map = fit$tmb_map,
    random = NULL, DLL = "sdmTMB", silent = TRUE
  )
  report <- obj$report(obj$par)
  a1 <- drop(fit$tmb_data$model_index_X_t[[1L]] %*% parameters$b_j) + 0.15
  a2 <- drop(fit$tmb_data$model_index_X_t[[2L]] %*% parameters$b_j2) - 0.05
  expected <- a1 + a2
  expected <- expected - mean(expected)

  expect_equal(report$model_index_t, expected, tolerance = 1e-10)
})

test_that("truncated delta model index uses the conditional positive mean", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = delta_truncated_nbinom2(), do_fit = FALSE
  )
  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 8 - 0.5
  parameters$b_j2[] <- rev(seq_along(parameters$b_j2)) / 9
  parameters$epsilon_st[] <- 0
  parameters$ln_phi[2L] <- log(2.3)

  obj <- TMB::MakeADFun(
    data = fit$tmb_data, parameters = parameters, map = fit$tmb_map,
    random = NULL, DLL = "sdmTMB", silent = TRUE
  )
  report <- obj$report(obj$par)
  a1 <- drop(fit$tmb_data$model_index_X_t[[1L]] %*% parameters$b_j)
  a2 <- drop(fit$tmb_data$model_index_X_t[[2L]] %*% parameters$b_j2)
  mu <- exp(a2)
  phi <- exp(parameters$ln_phi[2L])
  nonzero_probability <- 1 - (phi / (phi + mu))^phi
  q <- log(plogis(a1) * mu / nonzero_probability)
  expected <- q - mean(q)

  expect_equal(report$model_index_t, expected, tolerance = 1e-10)
})

test_that("prediction uses the stored model index at fitted times", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ depth_scaled + year_factor,
    spatial_varying = ~ 0 + depth_scaled + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = tweedie(link = "log"), do_fit = FALSE
  )
  rows <- unique(round(seq(nrow(d), 1L, length.out = 40L)))
  nd <- d[rows, ]
  prediction_data <- predict(fit, newdata = nd, return_tmb_data = TRUE)

  expect_identical(
    colnames(prediction_data$proj_z_i),
    c("depth_scaled", "model_index(year_factor)")
  )
  expect_true(all(prediction_data$proj_z_i[, 2L] == 0))

  parameters <- fit$tmb_params
  parameters$b_j[] <- seq_along(parameters$b_j) / 5
  parameters$epsilon_st[] <- 0
  parameters$zeta_s[] <- 0
  make_report <- function(parameters) {
    obj <- TMB::MakeADFun(
      data = prediction_data, parameters = parameters, map = fit$tmb_map,
      random = NULL, DLL = "sdmTMB", silent = TRUE
    )
    obj$report(obj$par)
  }
  report_zero <- make_report(parameters)
  parameters$zeta_s[, 2L, 1L] <- 0.6
  report_svc <- make_report(parameters)
  expected <- 0.6 * report_svc$model_index_t[prediction_data$proj_year + 1L]

  expect_equal(
    drop(report_svc$proj_eta - report_zero$proj_eta),
    expected,
    tolerance = 1e-9
  )

  nd_repeated <- nd[c(seq_len(nrow(nd)), 1L, 1L, 5L), ]
  repeated_data <- predict(fit, newdata = nd_repeated, return_tmb_data = TRUE)
  repeated_obj <- TMB::MakeADFun(
    data = repeated_data, parameters = parameters, map = fit$tmb_map,
    random = NULL, DLL = "sdmTMB", silent = TRUE
  )
  expect_equal(
    repeated_obj$report(repeated_obj$par)$model_index_t,
    report_svc$model_index_t,
    tolerance = 1e-10
  )
})

test_that("simulation recomputes the model index from simulated epsilon_st", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)

  one_part <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = tweedie(link = "log"), do_fit = FALSE
  )
  one_part_reports <- simulate(
    one_part, nsim = 2L, seed = 101, re_form = NA,
    return_tmb_report = TRUE, silent = TRUE
  )
  expect_length(one_part_reports, 2L)
  for (report in one_part_reports) {
    epsilon_mean <- apply(report$epsilon_st[, , 1L], 2L, mean)
    expected <- epsilon_mean - mean(epsilon_mean)
    expect_equal(report$model_index_t, expected, tolerance = 1e-8)
  }
  expect_false(isTRUE(all.equal(
    one_part_reports[[1L]]$model_index_t,
    one_part_reports[[2L]]$model_index_t
  )))

  delta <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = delta_gamma(), do_fit = FALSE
  )
  delta_sim <- simulate(
    delta, nsim = 2L, seed = 102, re_form = NA,
    newdata = d[seq_len(30L), ], silent = TRUE
  )
  expect_identical(dim(delta_sim), c(30L, 2L))
  expect_true(all(is.finite(delta_sim)))
})

test_that("forward project rejects model-index models", {
  d <- pcod
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ year_factor,
    spatial_varying = ~ 0 + model_index(year_factor),
    data = d, mesh = mesh, time = "year",
    family = tweedie(link = "log"), do_fit = FALSE
  )
  expect_error(
    project(fit, newdata = d[1L, ], nsim = 1L),
    "does not currently support.*model_index"
  )
  future <- d[1L, ]
  future$year <- max(d$year) + 1
  expect_error(
    predict(fit, newdata = future, return_tmb_data = TRUE),
    "Some new time values"
  )
})

test_that("optimized model-index fits converge and predict fitted rows consistently", {
  # kind of slow!
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("dplyr")

  d <- dogfish
  d$year_factor <- factor(d$year)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)

  fit_model_index <- function(family) {
    sdmTMB(
      catch_weight ~ 0 + year_factor,
      offset = log(d$area_swept),
      spatial_varying = ~ 0 + model_index(year_factor),
      data = d, mesh = mesh, time = "year", family = family,
      spatial = "off", spatiotemporal = "iid", silent = TRUE,
      control = sdmTMBcontrol(multiphase = FALSE))
  }
  suppressWarnings({
  fits <- list(
    tweedie = fit_model_index(tweedie(link = "log")),
    delta = fit_model_index(delta_lognormal(type = "poisson-link"))
  )})

  for (fit in fits) {
    expect_identical(fit$model$convergence, 0L)
    expect_true(fit$sd_report$pdHess)
    expect_true(all(is.finite(fit$gradients)))
    expect_lt(max(abs(fit$gradients)), 0.01)

    p_fitted <- predict(fit)
    p_original <- predict(fit, newdata = d)
    expect_equal(
      p_fitted[, prediction_columns, drop = FALSE],
      p_original[, prediction_columns, drop = FALSE],
      tolerance = 1e-5
    )

    reordered <- rev(seq_len(nrow(d)))
    p_reordered <- predict(fit, newdata = d[reordered, ])
    expect_equal(
      p_reordered[, prediction_columns, drop = FALSE],
      p_fitted[reordered, prediction_columns, drop = FALSE],
      tolerance = 1e-5
    )

    years <- sort(unique(d$year))
    keep <- d$year %in% years[c(2L, 4L)]
    p_subset <- predict(fit, newdata = d[keep, ])
    expect_equal(
      p_subset[, prediction_columns, drop = FALSE],
      p_fitted[keep, prediction_columns, drop = FALSE],
      tolerance = 1e-5
    )
  }

  # two-part comparison:
  # index:
  fit1 <- sdmTMB(
    catch_weight ~ 0 + year_factor,
    offset = log(d$area_swept),
    data = d, mesh = mesh, time = "year", family = tweedie(),
    spatial = "off", spatiotemporal = "iid"
  )
  p <- get_pars(fit1)
  df <- data.frame(year = sort(unique(d$year)), B = coef(fit0), eps_mean = apply(p$epsilon_st, 2, mean))
  df$index <- df$B + df$eps_mean
  df$index <- df$index - mean(df$index)

  d2 <- dplyr::left_join(d, df, by = join_by(year))

  # log density SVC
  fit2 <- sdmTMB(
    catch_weight ~ 0 + year_factor,
    spatial_varying = ~ 0 + index,
    offset = log(d$area_swept),
    data = d2, mesh = mesh, time = "year", family = tweedie(),
    spatial = "off", spatiotemporal = "iid"
  )

  pred_fancy <- predict(fits$tweedie)
  pred2part <- predict(fit2)

  z1 <- pred_fancy$`zeta_s_model_index(year_factor)`
  z2 <- pred2part$zeta_s_index
  plot(z1, z2)
  abline(0, 1)
  expect_gt(cor(z2, z1), 0.98)

  tidy(fit2, "ran_pars")
  tidy(fits$tweedie, "ran_pars")
})
