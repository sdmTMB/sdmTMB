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

test_that("multi-family fits with more than one family still fail early on main", {
  dat <- data.frame(
    y = c(1.2, 0, 2.4, 3.1),
    dist = c("gauss", "pois", "gauss", "pois")
  )

  expect_error(
    sdmTMB(
      y ~ 1,
      data = dat,
      spatial = "off",
      spatiotemporal = "off",
      family = list(gauss = gaussian(), pois = poisson()),
      distribution_column = "dist",
      do_fit = FALSE
    ),
    regexp = "not wired into `main` yet"
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
    if (is.null(pred_regular$est)) {
      expect_null(pred_list$est)
    } else {
      expect_equal(pred_list$est, pred_regular$est, tolerance = 1e-5)
    }
    expect_equal(pred_list$est1, pred_regular$est1, tolerance = 1e-5)
    expect_equal(pred_list$est2, pred_regular$est2, tolerance = 1e-5)
  }
})
