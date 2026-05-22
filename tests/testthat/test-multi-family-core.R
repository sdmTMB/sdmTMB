test_that("family_spec normalizes ordinary single-family fits", {
  dat <- data.frame(y = c(1.2, 2.4, 3.6))

  spec <- .build_family_spec(gaussian(), data = dat)

  expect_identical(spec$n_f, 1L)
  expect_identical(spec$n_m, 1L)
  expect_identical(spec$distribution_column, NULL)
  expect_equal(spec$family_id_i, rep(1L, nrow(dat)))
  expect_true(spec$active[1, 1])
  expect_equal(spec$family_code[1, 1], unname(.valid_family["gaussian"]))
  expect_equal(spec$link_code[1, 1], unname(.valid_link["identity"]))
  expect_identical(unname(spec$combine_kind), "single")
  expect_equal(spec$param_slot$ln_phi, 1L)
  expect_true(is.na(spec$param_slot$thetaf))

  y_out <- .family_spec_build_response(dat$y, spec)
  expect_identical(ncol(y_out), 1L)
  expect_equal(y_out[, 1], dat$y)
})

test_that("family_spec normalizes mixed family metadata", {
  fam <- list(
    gauss = gaussian(),
    delta = delta_gamma(type = "poisson-link"),
    stud = student()
  )
  dat <- data.frame(
    y = c(1.2, 0, 2.4, 3.1, 0.8),
    dist = c("gauss", "delta", "delta", "stud", "gauss")
  )

  spec <- .build_family_spec(
    family = fam,
    data = dat,
    distribution_column = "dist"
  )

  expect_identical(spec$n_f, 3L)
  expect_identical(spec$n_m, 2L)
  expect_equal(spec$family_id_i, c(1L, 2L, 2L, 3L, 1L))
  expect_equal(unname(spec$combine_kind), c("single", "poisson_link_delta", "single"))
  expect_equal(spec$family_code[1, 1], unname(.valid_family["gaussian"]))
  expect_equal(spec$family_code[2, 1], unname(.valid_family["binomial"]))
  expect_equal(spec$family_code[2, 2], unname(.valid_family["Gamma"]))
  expect_equal(spec$link_code[2, 1], unname(.valid_link["log"]))
  expect_equal(spec$link_code[2, 2], unname(.valid_link["log"]))
  expect_true(is.na(spec$family_code[1, 2]))
  expect_equal(spec$param_slot$ln_phi, c(1L, 2L, 3L))
  expect_equal(spec$param_slot$ln_student_df, c(NA_integer_, NA_integer_, 1L))

  y_out <- .family_spec_build_response(dat$y, spec)
  expect_equal(y_out[1, ], c(1.2, NA_real_))
  expect_equal(y_out[2, ], c(0, NA_real_))
  expect_equal(y_out[3, ], c(1, 2.4))
  expect_equal(y_out[4, ], c(3.1, NA_real_))
})

test_that("family_spec processes rowwise binomial-like responses regardless of family order", {
  dat <- data.frame(
    y = c(0.2, 1, 0, 1, 0.4, 2, 3.5),
    dist = c("binom", "binom", "betabinom", "betabinom", "betabinom", "betabinom", "gauss")
  )
  weights <- c(10, NA, 5, 5, 12, 7, 1)

  run_case <- function(fam) {
    spec <- .build_family_spec(
      family = fam,
      data = dat,
      distribution_column = "dist"
    )
    .family_spec_process_response(
      y_i = dat$y,
      size = rep(1, nrow(dat)),
      weights = weights,
      family_spec = spec
    )
  }

  res <- run_case(list(
    binom = binomial(),
    betabinom = betabinomial(),
    gauss = gaussian()
  ))
  res_reversed <- run_case(list(
    gauss = gaussian(),
    betabinom = betabinomial(),
    binom = binomial()
  ))

  expect_equal(res, res_reversed)

  expect_equal(res$y_i[1], 2)
  expect_equal(res$size[1], 10)
  expect_equal(res$weights[1], 1)

  expect_equal(res$y_i[2], 1)
  expect_equal(res$size[2], 1)
  expect_equal(res$weights[2], 1)

  expect_equal(res$y_i[3:6], c(0, 1, 4.8, 2))
  expect_equal(res$size[3:6], c(5, 5, 12, 7))
  expect_equal(res$weights[3:6], rep(1, 4))
})

test_that("multi-family beta-binomial weighted integer responses reach TMB as counts", {
  dat <- data.frame(
    y = c(0, 1, 0.25, 2, 3.5),
    x = seq(-1, 1, length.out = 5),
    dist = c("betabinom", "betabinom", "betabinom", "betabinom", "gauss")
  )
  weights <- c(4, 5, 8, 7, 1)

  fit_a <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian(), betabinom = betabinomial()),
    distribution_column = "dist",
    weights = weights,
    do_fit = FALSE
  )
  fit_b <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(betabinom = betabinomial(), gauss = gaussian()),
    distribution_column = "dist",
    weights = weights,
    do_fit = FALSE
  )

  bb_rows <- dat$dist == "betabinom"
  expect_equal(fit_a$tmb_data$size[bb_rows], c(4, 5, 8, 7))
  expect_equal(fit_b$tmb_data$size[bb_rows], c(4, 5, 8, 7))
  expect_equal(as.numeric(fit_a$tmb_data$y_i[bb_rows, 1]), c(0, 1, 2, 2))
  expect_equal(as.numeric(fit_b$tmb_data$y_i[bb_rows, 1]), c(0, 1, 2, 2))
})

test_that("object_family_spec can rebuild stored single-family fields", {
  dat <- data.frame(y = c(1.2, 0, 2.4, 3.1))
  obj <- list(
    family = list(gauss = gaussian()),
    data = dat
  )

  spec <- .object_family_spec(obj)

  expect_identical(spec$n_f, 1L)
  expect_identical(spec$n_m, 1L)
  expect_identical(spec$distribution_column, NULL)
  expect_equal(spec$family_id_i, rep(1L, nrow(dat)))
})

test_that("object_family_spec requires canonical metadata for old multi-family objects", {
  dat <- data.frame(
    y = c(1.2, 0, 2.4, 3.1),
    dist = c("gauss", "delta", "delta", "gauss")
  )
  obj <- list(
    family = list(gauss = gaussian(), delta = delta_gamma()),
    data = dat,
    distribution_column = "dist"
  )

  expect_error(
    .object_family_spec(obj, caller = "`predict()`"),
    regexp = "Refit this model with the current version of sdmTMB"
  )
})

test_that("multi-family fits build family-based TMB data", {
  dat <- data.frame(
    y = c(1.2, 0, 2.4, 3.1),
    x = c(-1, -0.3, 0.4, 1),
    dist = c("gauss", "delta", "delta", "gauss")
  )

  fit <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian(), delta = delta_gamma()),
    distribution_column = "dist",
    do_fit = FALSE
  )

  expect_equal(fit$tmb_data$obs_family_id, c(0L, 1L, 1L, 0L))
  expect_equal(fit$tmb_data$component_active, matrix(c(1L, 0L, 1L, 1L), nrow = 2L, byrow = TRUE))
  expect_equal(fit$tmb_data$combine_kind, c(0L, 1L))
  expect_equal(fit$tmb_data$ln_phi_slot, c(0L, 1L))
  expect_equal(fit$tmb_data$thetaf_slot, c(-1L, -1L))
  expect_equal(fit$tmb_data$y_i[1, ], c(1.2, NA_real_))
  expect_equal(fit$tmb_data$y_i[2, ], c(0, NA_real_))
  expect_equal(fit$tmb_data$y_i[3, ], c(1, 2.4))
  expect_equal(fit$tmb_data$y_i[4, ], c(3.1, NA_real_))
  expect_length(fit$tmb_params$ln_phi, 2L)
})

test_that("NA-dropped fits keep row-aligned family and offset data", {
  dat_single <- data.frame(
    y = c(1.2, 2.4, 3.6, 4.8, 6.0),
    x = c(-1, NA, 0, 0.5, 1)
  )

  fit_single <- sdmTMB(
    y ~ x,
    data = dat_single,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    offset = c(10, 20, 30, 40, 50),
    do_fit = FALSE
  )

  expect_equal(fit_single$tmb_data$offset_i, c(10, 30, 40, 50))
  expect_equal(fit_single$tmb_data$obs_family_id, rep(0L, 4L))
  expect_equal(fit_single$family_spec$family_id_i, rep(1L, 5L))
  expect_equal(fit_single$tmb_data$y_i[, 1], c(1.2, 3.6, 4.8, 6.0))

  dat_multi <- data.frame(
    y = c(1.2, 0.0, 2.1, 0.0, 3.5),
    x = c(-1, NA, -0.1, 0.4, 1),
    dist = c("gauss", "delta", "delta", "gauss", "delta")
  )

  fit_multi <- sdmTMB(
    y ~ x,
    data = dat_multi,
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian(), delta = delta_gamma()),
    distribution_column = "dist",
    offset = c(5, 6, 7, 8, 9),
    do_fit = FALSE
  )

  expect_equal(fit_multi$tmb_data$offset_i, c(5, 7, 8, 9))
  expect_equal(fit_multi$tmb_data$obs_family_id, c(0L, 1L, 0L, 1L))
  expect_equal(fit_multi$family_spec$family_id_i, c(1L, 2L, 2L, 1L, 2L))
  expect_equal(fit_multi$tmb_data$y_i[1, ], c(1.2, NA_real_))
  expect_equal(fit_multi$tmb_data$y_i[2, ], c(1.0, 2.1))
  expect_equal(fit_multi$tmb_data$y_i[3, ], c(0.0, NA_real_))
  expect_equal(fit_multi$tmb_data$y_i[4, ], c(1.0, 3.5))
})

test_that("mixed gaussian plus delta fits reach the unified TMB path", {
  set.seed(13)
  x <- seq(-1, 1, length.out = 90)
  dist <- rep(c("gauss", "delta", "delta"), length.out = length(x))
  gauss_rows <- dist == "gauss"
  delta_rows <- !gauss_rows

  y <- numeric(length(x))
  y[gauss_rows] <- 1 + 0.4 * x[gauss_rows] + stats::rnorm(sum(gauss_rows), sd = 0.15)
  present <- stats::rbinom(sum(delta_rows), size = 1, prob = stats::plogis(-0.3 + 0.7 * x[delta_rows]))
  y_pos <- stats::rgamma(sum(delta_rows), shape = 6, scale = exp(0.2 + 0.3 * x[delta_rows]) / 6)
  y[delta_rows] <- ifelse(present == 1, y_pos, 0)

  dat <- data.frame(y = y, x = x, dist = dist)
  ctrl <- sdmTMBcontrol(newton_loops = 1, getsd = FALSE)

  fit <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian(), delta = delta_gamma()),
    distribution_column = "dist",
    control = ctrl
  )

  expect_s3_class(fit, "sdmTMB")
  expect_identical(fit$family_spec$n_f, 2L)
  expect_identical(fit$family_spec$n_m, 2L)
  expect_equal(fit$model$convergence, 0L)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_error(residuals(fit), regexp = "simulation-based residual workflow")
})

.mixed_family_step4_fixture <- function() {
  set.seed(16)
  x <- seq(-1, 1, length.out = 90)
  dist <- rep(c("gauss", "delta", "delta"), length.out = length(x))
  gauss_rows <- dist == "gauss"
  delta_rows <- !gauss_rows

  y <- numeric(length(x))
  y[gauss_rows] <- 1 + 0.4 * x[gauss_rows] + stats::rnorm(sum(gauss_rows), sd = 0.15)
  present <- stats::rbinom(sum(delta_rows), size = 1, prob = stats::plogis(-0.3 + 0.7 * x[delta_rows]))
  y_pos <- stats::rgamma(sum(delta_rows), shape = 6, scale = exp(0.2 + 0.3 * x[delta_rows]) / 6)
  y[delta_rows] <- ifelse(present == 1, y_pos, 0)

  dat <- data.frame(y = y, x = x, dist = dist, year = 1L)
  fit <- sdmTMB(
    y ~ x,
    data = dat,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian(), delta = delta_gamma()),
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = TRUE)
  )

  nd <- data.frame(
    x = c(-0.8, -0.3, 0.2, 0.6, 0.9),
    dist = c("gauss", "delta", "gauss", "delta", "delta"),
    year = 1L
  )

  list(fit = fit, newdata = nd)
}

.mixed_family_standard_delta_fixture <- function() {
  set.seed(103)
  n_each <- 120L
  n <- n_each * 3L
  x <- runif(n, -1, 1)
  area <- rep(1, n)
  dist <- factor(
    rep(c("binomial", "poisson", "delta"), each = n_each),
    levels = c("binomial", "poisson", "delta")
  )

  y <- numeric(length(dist))
  i_bin <- dist == "binomial"
  i_poi <- dist == "poisson"
  i_del <- dist == "delta"

  y[i_bin] <- stats::rbinom(sum(i_bin), 1, stats::plogis(-0.4 + 0.6 * x[i_bin]))
  y[i_poi] <- stats::rpois(sum(i_poi), exp(0.3 + 0.5 * x[i_poi]))

  present <- stats::rbinom(sum(i_del), 1, stats::plogis(-0.2 + 0.3 * x[i_del]))
  y_pos <- stats::rgamma(sum(i_del), shape = 4, scale = exp(0.4 + 0.5 * x[i_del]) / 4)
  y[i_del] <- ifelse(present == 1, y_pos, 0)

  dat <- data.frame(
    y = y,
    x = x,
    y_coord = 0,
    year = 1L,
    dist = dist,
    area_swept = area
  )

  ctrl <- sdmTMBcontrol(multiphase = FALSE, newton_loops = 0, getsd = FALSE)

  fit_mix <- sdmTMB(
    y ~ 0 + dist + dist:x,
    data = dat,
    family = list(
      binomial = binomial(link = "logit"),
      poisson = poisson(link = "log"),
      delta = delta_gamma()
    ),
    distribution_column = "dist",
    spatial = "off",
    spatiotemporal = "off",
    weights = rep(1, nrow(dat)),
    offset = log(dat$area_swept),
    control = ctrl
  )

  fit_bin <- sdmTMB(
    y ~ x,
    data = subset(dat, dist == "binomial"),
    family = binomial(link = "logit"),
    spatial = "off",
    spatiotemporal = "off",
    weights = rep(1, sum(dat$dist == "binomial")),
    control = ctrl
  )

  fit_poi <- sdmTMB(
    y ~ x,
    data = subset(dat, dist == "poisson"),
    family = poisson(link = "log"),
    spatial = "off",
    spatiotemporal = "off",
    offset = log(subset(dat, dist == "poisson")$area_swept),
    control = ctrl
  )

  fit_del <- sdmTMB(
    y ~ x,
    data = subset(dat, dist == "delta"),
    family = delta_gamma(),
    spatial = "off",
    spatiotemporal = "off",
    offset = log(subset(dat, dist == "delta")$area_swept),
    control = ctrl
  )

  make_newdata <- function(dist_label, area_swept = rep(1, 9L)) {
    data.frame(
      x = seq(-0.9, 0.9, length.out = 9),
      y_coord = 0,
      year = 1L,
      dist = factor(dist_label, levels = levels(dat$dist)),
      area_swept = area_swept
    )
  }

  list(
    fit_mix = fit_mix,
    fit_bin = fit_bin,
    fit_poi = fit_poi,
    fit_del = fit_del,
    nd_bin = make_newdata("binomial"),
    nd_poi = make_newdata("poisson"),
    nd_del = make_newdata("delta")
  )
}

test_that("mixed-family predict returns rowwise combined and masked component outputs", {
  fixture <- .mixed_family_step4_fixture()
  gauss_rows <- fixture$newdata$dist == "gauss"
  delta_rows <- !gauss_rows

  pred_link <- predict(fixture$fit, newdata = fixture$newdata, type = "link")
  expect_equal(pred_link$est[gauss_rows], pred_link$est1[gauss_rows], tolerance = 1e-6)
  expect_true(all(is.na(pred_link$est2[gauss_rows])))
  expect_true(all(is.finite(pred_link$est2[delta_rows])))
  expect_equal(
    pred_link$est[delta_rows],
    log(stats::plogis(pred_link$est1[delta_rows]) * exp(pred_link$est2[delta_rows])),
    tolerance = 1e-6
  )

  pred_lp1 <- predict(fixture$fit, newdata = fixture$newdata, type = "link", model = 1)
  pred_lp2 <- predict(fixture$fit, newdata = fixture$newdata, type = "link", model = 2)
  expect_equal(pred_lp1$est, pred_lp1$est1, tolerance = 1e-6)
  expect_equal(pred_lp2$est, pred_lp2$est2, tolerance = 1e-6)
  expect_true(all(is.na(pred_lp2$est[gauss_rows])))

  pred_resp <- predict(fixture$fit, newdata = fixture$newdata, type = "response")
  expect_equal(pred_resp$est[gauss_rows], pred_resp$est1[gauss_rows], tolerance = 1e-6)
  expect_equal(
    pred_resp$est[delta_rows],
    pred_resp$est1[delta_rows] * pred_resp$est2[delta_rows],
    tolerance = 1e-6
  )

  pred_se <- predict(fixture$fit, newdata = fixture$newdata, type = "link", se_fit = TRUE)
  expect_true(all(c("est", "est1", "est2", "est_se") %in% names(pred_se)))
  expect_equal(pred_se$est[gauss_rows], pred_se$est1[gauss_rows], tolerance = 1e-6)

  expect_error(
    predict(fixture$fit, newdata = fixture$newdata, type = "response", se_fit = TRUE),
    regexp = "not yet supported"
  )
})

test_that("mixed-family standard delta rows match standalone sdmTMB references", {
  fixture <- .mixed_family_standard_delta_fixture()
  nd_single <- subset(fixture$nd_del, select = c(x, y_coord, year, area_swept))

  pred_mix_bin <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_bin,
    type = "response",
    offset = rep(0, nrow(fixture$nd_bin))
  )
  pred_bin <- predict(fixture$fit_bin, newdata = nd_single, type = "response")

  pred_mix_poi <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_poi,
    type = "response",
    offset = rep(0, nrow(fixture$nd_poi))
  )
  pred_poi <- predict(
    fixture$fit_poi,
    newdata = nd_single,
    type = "response",
    offset = rep(0, nrow(nd_single))
  )

  pred_mix_del <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_del,
    type = "response",
    offset = rep(0, nrow(fixture$nd_del))
  )
  pred_del <- predict(
    fixture$fit_del,
    newdata = nd_single,
    type = "response",
    offset = rep(0, nrow(nd_single))
  )
  pred_mix_del1 <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_del,
    type = "response",
    model = 1,
    offset = rep(0, nrow(fixture$nd_del))
  )
  pred_del1 <- predict(
    fixture$fit_del,
    newdata = nd_single,
    type = "response",
    model = 1,
    offset = rep(0, nrow(nd_single))
  )
  pred_mix_del2 <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_del,
    type = "response",
    model = 2,
    offset = rep(0, nrow(fixture$nd_del))
  )
  pred_del2 <- predict(
    fixture$fit_del,
    newdata = nd_single,
    type = "response",
    model = 2,
    offset = rep(0, nrow(nd_single))
  )

  expect_equal(pred_mix_bin$est, pred_bin$est, tolerance = 1e-5)
  expect_equal(pred_mix_poi$est, pred_poi$est, tolerance = 1e-5)
  expect_equal(pred_mix_del$est, pred_del$est, tolerance = 1e-5)
  expect_equal(pred_mix_del1$est, pred_del1$est, tolerance = 1e-5)
  expect_equal(pred_mix_del2$est, pred_del2$est, tolerance = 1e-5)
})

test_that("mixed-family standard delta offsets only affect the positive component", {
  fixture <- .mixed_family_standard_delta_fixture()

  nd_off <- fixture$nd_del
  nd_off$area_swept <- exp(seq(log(0.7), log(1.3), length.out = nrow(nd_off)))
  off <- log(nd_off$area_swept)
  nd_single_off <- subset(nd_off, select = c(x, y_coord, year, area_swept))

  pred_mix_del1 <- predict(
    fixture$fit_mix,
    newdata = fixture$nd_del,
    type = "response",
    model = 1,
    offset = rep(0, nrow(fixture$nd_del))
  )
  pred_mix_del1_off <- predict(
    fixture$fit_mix,
    newdata = nd_off,
    type = "response",
    model = 1,
    offset = off
  )
  pred_mix_del2_off <- predict(
    fixture$fit_mix,
    newdata = nd_off,
    type = "response",
    model = 2,
    offset = off
  )
  pred_mix_del_off <- predict(
    fixture$fit_mix,
    newdata = nd_off,
    type = "response",
    offset = off
  )

  pred_del1_off <- predict(
    fixture$fit_del,
    newdata = nd_single_off,
    type = "response",
    model = 1,
    offset = off
  )
  pred_del2_off <- predict(
    fixture$fit_del,
    newdata = nd_single_off,
    type = "response",
    model = 2,
    offset = off
  )
  pred_del_off <- predict(
    fixture$fit_del,
    newdata = nd_single_off,
    type = "response",
    offset = off
  )

  expect_equal(pred_mix_del1_off$est, pred_mix_del1$est, tolerance = 1e-8)
  expect_equal(pred_mix_del1_off$est, pred_del1_off$est, tolerance = 1e-5)
  expect_equal(pred_mix_del2_off$est, pred_del2_off$est, tolerance = 1e-5)
  expect_equal(pred_mix_del_off$est, pred_del_off$est, tolerance = 1e-5)
})

test_that("mixed-family simulate combines rows the same way as predict", {
  fixture <- .mixed_family_step4_fixture()
  gauss_rows <- fixture$newdata$dist == "gauss"
  delta_rows <- !gauss_rows

  sim_combined <- simulate(
    fixture$fit,
    nsim = 2,
    seed = 101,
    newdata = fixture$newdata,
    model = NA,
    type = "mle-eb",
    silent = TRUE
  )
  sim_lp1 <- simulate(
    fixture$fit,
    nsim = 2,
    seed = 101,
    newdata = fixture$newdata,
    model = 1,
    type = "mle-eb",
    silent = TRUE
  )
  sim_lp2 <- simulate(
    fixture$fit,
    nsim = 2,
    seed = 101,
    newdata = fixture$newdata,
    model = 2,
    type = "mle-eb",
    silent = TRUE
  )

  expect_equal(sim_combined[gauss_rows, , drop = FALSE], sim_lp1[gauss_rows, , drop = FALSE], tolerance = 1e-10)
  expect_true(all(is.na(sim_lp2[gauss_rows, , drop = FALSE])))
  expect_equal(
    sim_combined[delta_rows, , drop = FALSE],
    sim_lp1[delta_rows, , drop = FALSE] * sim_lp2[delta_rows, , drop = FALSE],
    tolerance = 1e-10
  )

  expect_error(project(fixture$fit, newdata = fixture$newdata, nsim = 1), regexp = "not yet supported")
})

test_that("mixed-family weighted-average and EAO helpers use rowwise families", {
  fixture <- .mixed_family_step4_fixture()
  area <- c(1, 2, 1, 3, 2)

  pred_resp <- predict(fixture$fit, newdata = fixture$newdata, type = "response")
  pred_obj <- predict(fixture$fit, newdata = fixture$newdata, return_tmb_object = TRUE)

  wa <- get_weighted_average(
    pred_obj,
    vector = fixture$newdata$x,
    area = area,
    bias_correct = FALSE
  )
  manual_total <- sum(pred_resp$est * area)
  manual_wa <- sum(pred_resp$est * fixture$newdata$x * area) / manual_total
  expect_equal(wa$est, manual_wa, tolerance = 1e-6)

  eao <- get_eao(
    pred_obj,
    area = area,
    bias_correct = FALSE
  )
  manual_mean_dens <- sum(pred_resp$est * pred_resp$est) / sum(pred_resp$est)
  manual_eao <- manual_total / manual_mean_dens
  expect_equal(eao$est, manual_eao, tolerance = 1e-6)
})

test_that("mixed-family tidy and print report component and family summaries", {
  fixture <- .mixed_family_step4_fixture()

  tidy_lp2 <- tidy(fixture$fit, model = 2)
  expect_true(all(c("term", "estimate", "std.error") %in% names(tidy_lp2)))

  tidy_ran <- tidy(fixture$fit, effects = "ran_pars", model = 1)
  expect_true(any(tidy_ran$group_name == "gauss" & tidy_ran$term == "phi"))
  expect_true(any(tidy_ran$group_name == "delta" & tidy_ran$term == "phi"))

  printed <- paste(capture.output(print(fixture$fit)), collapse = "\n")
  expect_match(printed, "Linear predictor 1", fixed = TRUE)
  expect_match(printed, "Families:", fixed = TRUE)
  expect_match(printed, "distribution column = 'dist'", fixed = TRUE)
})

test_that("mixed-family methods use family-spec routing and guard unsupported summaries", {
  fixture <- .mixed_family_step4_fixture()
  fit <- fixture$fit

  pred_resp <- predict(fit, newdata = fit$data, type = "response")
  expect_equal(fitted(fit), pred_resp$est, tolerance = 1e-6)

  fam <- family(fit)
  expect_true(is.list(fam))
  expect_identical(names(fam), c("gauss", "delta"))

  vc_lp2 <- vcov(fit, model = 2)
  expect_identical(colnames(vc_lp2), rownames(vc_lp2))
  expect_equal(colnames(vc_lp2), tidy(fit, model = 2)$term)

  expect_error(sigma(fit), regexp = "not yet supported")
  expect_error(deviance(fit), regexp = "not yet supported")
  expect_error(spread_sims(fit), regexp = "not yet supported")
  expect_error(gather_sims(fit), regexp = "not yet supported")
})

test_that("multi-family do_index guard fails before fit setup continues", {
  dat <- data.frame(
    y = c(1.2, 0, 2.4, 3.1),
    x = c(-1, -0.3, 0.4, 1),
    dist = c("gauss", "delta", "delta", "gauss")
  )

  expect_error(
    sdmTMB(
      y ~ x,
      data = dat,
      spatial = "off",
      spatiotemporal = "off",
      family = list(gauss = gaussian(), delta = delta_gamma()),
      distribution_column = "dist",
      do_index = TRUE,
      do_fit = FALSE,
      predict_args = list(newdata = dat),
      index_args = list(area = 1)
    ),
    regexp = "not yet supported"
  )
})

test_that("named single-family list gaussian matches ordinary gaussian fit", {
  set.seed(11)
  x <- seq(-1, 1, length.out = 40)
  dat <- data.frame(
    y = 0.4 + 0.7 * x + stats::rnorm(length(x), sd = 0.15),
    x = x
  )
  ctrl <- sdmTMBcontrol(newton_loops = 0, getsd = FALSE)

  fit_regular <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    control = ctrl
  )
  fit_list <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(gauss = gaussian()),
    control = ctrl
  )

  expect_identical(fit_regular$family_spec$n_f, 1L)
  expect_identical(fit_list$family_spec$n_f, 1L)
  expect_equal(as.numeric(logLik(fit_list)), as.numeric(logLik(fit_regular)), tolerance = 1e-6)

  nd <- data.frame(x = seq(-0.8, 0.8, length.out = 12))
  expect_equal(
    predict(fit_list, newdata = nd, type = "link")$est,
    predict(fit_regular, newdata = nd, type = "link")$est,
    tolerance = 1e-6
  )
  expect_equal(
    predict(fit_list, newdata = nd, type = "response")$est,
    predict(fit_regular, newdata = nd, type = "response")$est,
    tolerance = 1e-6
  )
})

test_that("named single-family list delta matches ordinary delta predictions", {
  set.seed(12)
  x <- seq(-1, 1, length.out = 80)
  eta1 <- -0.2 + 0.6 * x
  eta2 <- 0.3 + 0.4 * x
  present <- stats::rbinom(length(x), size = 1, prob = stats::plogis(eta1))
  y_pos <- stats::rgamma(length(x), shape = 8, scale = exp(eta2) / 8)
  dat <- data.frame(
    y = ifelse(present == 1, y_pos, 0),
    x = x
  )
  ctrl <- sdmTMBcontrol(newton_loops = 0, getsd = FALSE)

  fit_regular <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    control = ctrl
  )
  fit_list <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = list(delta = delta_gamma()),
    control = ctrl
  )

  expect_identical(fit_list$family_spec$n_m, 2L)

  nd <- data.frame(x = seq(-0.8, 0.8, length.out = 12))
  for (type in c("link", "response")) {
    pred_regular <- predict(fit_regular, newdata = nd, type = type, model = NA)
    pred_list <- predict(fit_list, newdata = nd, type = type, model = NA)
    expect_equal(pred_list$est, pred_regular$est, tolerance = 1e-5)
    expect_equal(pred_list$est1, pred_regular$est1, tolerance = 1e-5)
    expect_equal(pred_list$est2, pred_regular$est2, tolerance = 1e-5)
  }
})

test_that("delta sdreport exposes generic combined prediction names", {
  set.seed(14)
  x <- seq(-1, 1, length.out = 60)
  eta1 <- -0.1 + 0.5 * x
  eta2 <- 0.2 + 0.3 * x
  present <- stats::rbinom(length(x), size = 1, prob = stats::plogis(eta1))
  y_pos <- stats::rgamma(length(x), shape = 7, scale = exp(eta2) / 7)
  dat <- data.frame(
    y = ifelse(present == 1, y_pos, 0),
    x = x
  )
  ctrl <- sdmTMBcontrol(newton_loops = 0, getsd = FALSE)

  fit <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    control = ctrl
  )

  nd <- data.frame(x = seq(-0.8, 0.8, length.out = 10))

  pred_obj <- predict(fit, newdata = nd, se_fit = TRUE, return_tmb_object = TRUE)
  pred_sr <- TMB::sdreport(pred_obj$obj, bias.correct = FALSE)
  pred_rep <- as.list(pred_sr, "Estimate", report = TRUE)

  expect_true("proj_eta_combined" %in% names(pred_rep))
  expect_false("proj_eta_delta" %in% names(pred_rep))

  pop_pred_obj <- predict(
    fit,
    newdata = nd,
    se_fit = TRUE,
    re_form = NA,
    return_tmb_object = TRUE
  )
  pop_pred_sr <- TMB::sdreport(pop_pred_obj$obj, bias.correct = FALSE)
  pop_pred_rep <- as.list(pop_pred_sr, "Estimate", report = TRUE)

  expect_true("proj_fe_combined" %in% names(pop_pred_rep))
  expect_false("proj_rf_delta" %in% names(pop_pred_rep))
})

.poisson_link_delta_step3_fixture <- function() {
  set.seed(15)
  x <- seq(-1, 1, length.out = 120)
  eta1 <- -0.7 + 0.6 * x
  eta2 <- 0.2 + 0.3 * x
  p <- 1 - exp(-exp(eta1))
  mu_pos <- exp(eta1 + eta2) / p
  present <- stats::rbinom(length(x), size = 1, prob = p)
  y_pos <- stats::rgamma(length(x), shape = 6, scale = mu_pos / 6)

  dat <- data.frame(
    y = ifelse(present == 1, y_pos, 0),
    x = x,
    year = 1L
  )
  ctrl <- sdmTMBcontrol(newton_loops = 1, getsd = FALSE)

  fit <- sdmTMB(
    y ~ x,
    data = dat,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(type = "poisson-link"),
    control = ctrl
  )

  list(fit = fit, data = dat)
}

test_that("poisson-link delta ordinary se_fit uses generic combined report", {
  fixture <- .poisson_link_delta_step3_fixture()
  nd <- data.frame(
    x = seq(-0.8, 0.8, length.out = 12),
    year = 1L
  )

  pred <- predict(
    fixture$fit,
    newdata = nd,
    type = "link",
    se_fit = TRUE
  )
  expect_equal(as.numeric(pred$est), as.numeric(pred$est1 + pred$est2), tolerance = 1e-6)

  pred_obj <- predict(
    fixture$fit,
    newdata = nd,
    type = "link",
    se_fit = TRUE,
    return_tmb_object = TRUE
  )
  pred_sr <- TMB::sdreport(pred_obj$obj, bias.correct = FALSE)
  pred_rep <- as.list(pred_sr, "Estimate", report = TRUE)

  expect_true("proj_eta_combined" %in% names(pred_rep))
  expect_false("proj_eta_delta" %in% names(pred_rep))
})

test_that("poisson-link delta population se_fit uses generic combined report", {
  fixture <- .poisson_link_delta_step3_fixture()
  nd <- data.frame(
    x = seq(-0.8, 0.8, length.out = 12),
    year = 1L
  )

  pred <- predict(
    fixture$fit,
    newdata = nd,
    type = "link",
    se_fit = TRUE,
    re_form = NA
  )
  expect_equal(as.numeric(pred$est), as.numeric(pred$est1 + pred$est2), tolerance = 1e-6)

  pred_obj <- predict(
    fixture$fit,
    newdata = nd,
    type = "link",
    se_fit = TRUE,
    re_form = NA,
    return_tmb_object = TRUE
  )
  pred_sr <- TMB::sdreport(pred_obj$obj, bias.correct = FALSE)
  pred_rep <- as.list(pred_sr, "Estimate", report = TRUE)

  expect_true("proj_fe_combined" %in% names(pred_rep))
  expect_false("proj_rf_delta" %in% names(pred_rep))
})

test_that("poisson-link delta get_index matches constant-grid response prediction", {
  fixture <- .poisson_link_delta_step3_fixture()
  nd <- data.frame(
    x = rep(0, 40),
    year = 1L
  )

  pred_obj <- predict(
    fixture$fit,
    newdata = nd,
    return_tmb_object = TRUE
  )
  ind <- get_index(
    pred_obj,
    area = rep(1 / nrow(nd), nrow(nd)),
    bias_correct = FALSE
  )
  pred_response <- predict(
    fixture$fit,
    newdata = nd[1, , drop = FALSE],
    type = "response"
  )

  expect_equal(ind$est, pred_response$est[1], tolerance = 1e-3)
})
