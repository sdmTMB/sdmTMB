test_that("family_spec normalizes ordinary single-family fits", {
  dat <- data.frame(y = c(1.2, 2.4, 3.6))

  spec <- sdmTMB:::.build_family_spec(gaussian(), data = dat)

  expect_identical(spec$n_f, 1L)
  expect_identical(spec$n_m, 1L)
  expect_identical(spec$distribution_column, NULL)
  expect_equal(spec$family_id_i, rep(1L, nrow(dat)))
  expect_true(spec$active[1, 1])
  expect_equal(spec$family_code[1, 1], unname(sdmTMB:::.valid_family["gaussian"]))
  expect_equal(spec$link_code[1, 1], unname(sdmTMB:::.valid_link["identity"]))
  expect_identical(unname(spec$combine_kind), "single")
  expect_equal(spec$param_slot$ln_phi, 1L)
  expect_true(is.na(spec$param_slot$thetaf))

  y_out <- sdmTMB:::.family_spec_build_response(dat$y, spec)
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

  spec <- sdmTMB:::.build_family_spec(
    family = fam,
    data = dat,
    distribution_column = "dist"
  )

  expect_identical(spec$n_f, 3L)
  expect_identical(spec$n_m, 2L)
  expect_equal(spec$family_id_i, c(1L, 2L, 2L, 3L, 1L))
  expect_equal(unname(spec$combine_kind), c("single", "poisson_link_delta", "single"))
  expect_equal(spec$family_code[1, 1], unname(sdmTMB:::.valid_family["gaussian"]))
  expect_equal(spec$family_code[2, 1], unname(sdmTMB:::.valid_family["binomial"]))
  expect_equal(spec$family_code[2, 2], unname(sdmTMB:::.valid_family["Gamma"]))
  expect_equal(spec$link_code[2, 1], unname(sdmTMB:::.valid_link["log"]))
  expect_equal(spec$link_code[2, 2], unname(sdmTMB:::.valid_link["log"]))
  expect_true(is.na(spec$family_code[1, 2]))
  expect_equal(spec$param_slot$ln_phi, c(1L, 2L, 3L))
  expect_equal(spec$param_slot$ln_student_df, c(NA_integer_, NA_integer_, 1L))

  y_out <- sdmTMB:::.family_spec_build_response(dat$y, spec)
  expect_equal(y_out[1, ], c(1.2, NA_real_))
  expect_equal(y_out[2, ], c(0, NA_real_))
  expect_equal(y_out[3, ], c(1, 2.4))
  expect_equal(y_out[4, ], c(3.1, NA_real_))
})

test_that("family_spec processes rowwise binomial-like responses", {
  fam <- list(
    binom = binomial(),
    betabinom = betabinomial(),
    gauss = gaussian()
  )
  dat <- data.frame(
    y = c(0.2, 1, 0.4, 2, 3.5),
    dist = c("binom", "binom", "betabinom", "betabinom", "gauss")
  )
  spec <- sdmTMB:::.build_family_spec(
    family = fam,
    data = dat,
    distribution_column = "dist"
  )

  res <- sdmTMB:::.family_spec_process_response(
    y_i = dat$y,
    size = rep(1, nrow(dat)),
    weights = c(10, 4, 12, 7, 1),
    family_spec = spec
  )

  expect_equal(res$y_i[1], 2)
  expect_equal(res$size[1], 10)
  expect_equal(res$weights[1], 1)

  expect_equal(res$y_i[2], 1)
  expect_equal(res$size[2], 1)
  expect_equal(res$weights[2], 4)

  expect_equal(res$y_i[3], 4.8)
  expect_equal(res$size[3], 12)
  expect_equal(res$weights[3], 1)

  expect_equal(res$y_i[4], 2)
  expect_equal(res$size[4], 7)
  expect_equal(res$weights[4], 1)
})

test_that("multi-family fits build family-based TMB payloads", {
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
  expect_error(predict(fit), regexp = "not yet supported")
  expect_error(residuals(fit), regexp = "not yet supported")
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
    if (is.null(pred_regular$est)) {
      expect_null(pred_list$est)
    } else {
      expect_equal(pred_list$est, pred_regular$est, tolerance = 1e-5)
    }
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
