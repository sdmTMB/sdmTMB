test_that("dispformula works for non-multi-family models", {
  d <- pcod[pcod$year %in% c(2003, 2004, 2005, 2006), ]

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = d,
    time = "year",
    family = tweedie(link = "log"),
    spatial = "off",
    spatiotemporal = "off",
    dispformula = ~ 0 + as.factor(year),
    control = sdmTMBcontrol(newton_loops = 0)
  )

  p <- get_pars(fit)
  expect_true("b_disp_k" %in% names(p))
  expect_gt(length(p$b_disp_k), 1L)
  expect_identical(fit$tmb_data$has_dispersion_model, 1L)
  expect_true(isTRUE(fit$has_dispformula))

  td_disp <- tidy(fit, effects = "dispersion")
  expect_true(nrow(td_disp) > 1L)
  expect_true(all(c("term", "estimate", "std.error") %in% names(td_disp)))
  expect_true(all(grepl("as.factor\\(year\\)", td_disp$term)))

  td_re <- tidy(fit, effects = "ran_pars")
  expect_false("phi" %in% td_re$term)

  expect_output(print(fit), "Dispersion model:")

  expect_error(
    sigma(fit),
    regexp = "not available when `dispformula` is used"
  )

  expect_error(
    residuals(fit, type = "mle-eb"),
    regexp = "dharma_residuals"
  )
  expect_equal(length(residuals(fit, type = "response")), nrow(fit$data))
})

test_that("dispformula is guarded for multi-family models", {
  d <- pcod_2011
  d$dist <- ifelse(d$year %% 2 == 0, "g", "tw")

  expect_error(
    sdmTMB(
      density ~ 1,
      data = d,
      family = list(
        g = gaussian(),
        tw = tweedie(link = "log")
      ),
      distribution_column = "dist",
      spatial = "off",
      spatiotemporal = "off",
      dispformula = ~ 0 + as.factor(year),
      do_fit = FALSE
    ),
    regexp = "not supported for multi-family"
  )
})

test_that("dispformula is ignored for families without phi", {
  d <- pcod_2011
  d$present <- as.integer(d$density > 0)

  expect_warning(
    fit <- sdmTMB(
      present ~ depth_scaled,
      data = d,
      family = binomial(),
      spatial = "off",
      spatiotemporal = "off",
      dispformula = ~ 0 + as.factor(year),
      do_fit = FALSE
    ),
    regexp = "ignored"
  )

  expect_false(isTRUE(fit$has_dispformula))
  expect_identical(fit$tmb_data$has_dispersion_model, 0L)
})

test_that("dispformula works with simulation and index prediction pathways", {
  skip_on_cran()
  d1 <- pcod[pcod$year == 2003, ]
  d2 <- pcod[pcod$year == 2004, ]
  d <- rbind(d1[seq_len(60), ], d2[seq_len(60), ])
  mesh <- make_mesh(d, c("X", "Y"), n_knots = 20)

  fit <- sdmTMB(
    density ~ 1,
    data = d,
    mesh = mesh,
    time = "year",
    family = tweedie(link = "log"),
    spatiotemporal = "off",
    dispformula = ~ 0 + as.factor(year),
    control = sdmTMBcontrol(newton_loops = 0)
  )

  nd <- d[seq_len(80), c("X", "Y", "year"), drop = FALSE]

  sims <- simulate(fit, nsim = 2L, newdata = nd, silent = TRUE)
  expect_equal(dim(sims), c(nrow(nd), 2L))

  pred <- predict(fit, newdata = nd, return_tmb_object = TRUE)
  idx <- get_index(pred, bias_correct = FALSE, area = 1)
  expect_true(all(c("year", "est", "lwr", "upr") %in% names(idx)))
  expect_true(nrow(idx) >= 1L)
})

test_that("dispformula works across supported single-family likelihoods", {
  skip_on_cran()
  set.seed(1)

  d <- pcod[pcod$year %in% c(2003, 2004, 2005, 2006), ]
  d <- d[seq_len(min(140L, nrow(d))), ]
  d$density_pos <- pmax(d$density, 1e-3)
  d$count <- as.integer(round(d$density * 5))
  d$trials <- sample(2:6, nrow(d), replace = TRUE)
  linpred <- -0.3 + 0.4 * as.numeric(scale(d$depth_scaled))
  d$success <- stats::rbinom(nrow(d), size = d$trials, prob = stats::plogis(linpred))
  d$prop <- d$success / d$trials

  specs <- list(
    list(name = "gaussian", formula = density ~ depth_scaled, family = gaussian(), weights = NULL),
    list(name = "Gamma", formula = density_pos ~ depth_scaled, family = Gamma(link = "log"), weights = NULL),
    list(name = "tweedie", formula = density ~ depth_scaled, family = tweedie(link = "log"), weights = NULL),
    list(name = "nbinom1", formula = count ~ depth_scaled, family = nbinom1(), weights = NULL),
    list(name = "nbinom2", formula = count ~ depth_scaled, family = nbinom2(), weights = NULL),
    list(name = "betabinomial", formula = prop ~ depth_scaled, family = betabinomial(), weights = d$trials)
  )

  for (spec in specs) {
    fit <- sdmTMB(
      formula = spec$formula,
      data = d,
      family = spec$family,
      weights = spec$weights,
      spatial = "off",
      spatiotemporal = "off",
      dispformula = ~ depth_scaled,
      control = sdmTMBcontrol(newton_loops = 0)
    )

    expect_true(isTRUE(fit$has_dispformula), info = spec$name)
    td <- tidy(fit, effects = "dispersion")
    expect_true(nrow(td) > 0L, info = spec$name)

    pred <- predict(fit, newdata = d[seq_len(30), , drop = FALSE], type = "response")
    expect_true(nrow(pred) == 30L, info = spec$name)
  }
})

test_that("dispformula works for regular delta models", {
  skip_on_cran()

  d <- pcod[pcod$year %in% c(2003, 2004, 2005, 2006), ]
  d <- d[seq_len(min(160L, nrow(d))), ]

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = d,
    family = delta_gamma(),
    spatial = "off",
    spatiotemporal = "off",
    dispformula = ~ depth_scaled,
    control = sdmTMBcontrol(newton_loops = 0)
  )

  expect_true(isTRUE(fit$has_dispformula))
  td1 <- tidy(fit, effects = "dispersion", model = 1)
  td2 <- tidy(fit, effects = "dispersion", model = 2)
  expect_equal(nrow(td1), 0L)
  expect_gt(nrow(td2), 0L)

  p <- predict(fit, newdata = d[seq_len(40), , drop = FALSE], type = "response")
  expect_true(all(c("est", "est1", "est2") %in% names(p)))
})
