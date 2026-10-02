test_that("Link/response type works", {

  # see also test-delta-population-predictions.R
  # https://github.com/sdmTMB/sdmTMB/issues/110
  skip_on_cran()

  fit <- sdmTMB(
    density ~ 1,
    family = tweedie(),
    data = pcod_2011,
    mesh = pcod_mesh_2011
  )

  p <- predict(fit, type = "link")
  expect_lt(mean(p$est), 5)

  p <- predict(fit, type = "link", re_form = NA)
  expect_lt(mean(p$est), 5)

  p <- predict(fit, type = "response")
  expect_gt(mean(p$est), 20)

  p <- predict(fit, type = "response", re_form = NA)
  expect_gt(mean(p$est), 20)

  fit_delt <- sdmTMB(
    density ~ 1,
    family = delta_gamma(),
    data = pcod_2011,
    mesh = pcod_mesh_2011
  )

  p <- predict(fit_delt, type = "link")
  # expect_false("est" %in% names(p))
  expect_lt(mean(p$est1), 1)
  expect_lt(mean(p$est2), 5)

  p <- predict(fit_delt, type = "link", re_form = NA)
  # expect_false("est" %in% names(p))
  expect_lt(mean(p$est1), 1)
  expect_lt(mean(p$est2), 5)

  p <- predict(fit_delt, type = "response")
  expect_gt(mean(p$est), 20)
  expect_gt(mean(p$est1), 0)
  expect_gt(mean(p$est2), 30)

  p <- predict(fit_delt, type = "response", re_form = NA)
  expect_gt(mean(p$est), 20)
  expect_gt(mean(p$est1), 0)
  expect_gt(mean(p$est2), 30)

  # with std. error
  expect_message(p <- predict(fit, type = "link", re_form = NULL, se_fit = TRUE), regexp = "slow")
  mean(p$est)

  p <- predict(fit_delt, type = "link", re_form = NA, se_fit = TRUE)
  mean(p$est)

  expect_error(
    predict(fit_delt, type = "response", re_form = NA, se_fit = TRUE), regexp = "link scale"
  )
})

test_that("Response prediction works as reported in https://github.com/sdmTMB/sdmTMB/issues/160#issuecomment-1380920333", {
  skip_on_cran()

  fit_dg <- sdmTMB(density ~ 1 + s(depth),
    data = pcod_2011,
    mesh = pcod_mesh_2011,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma()
  )

  nd <- data.frame(
    depth = seq(min(pcod$depth),
      max(pcod$depth),
      length.out = 20
    ),
    X = 0,
    Y = 0,
    year = 2015L # a chosen year
  )

  p1 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = TRUE,
    type = "link"
  )

  expect_true(all(c("est", "est1", "est2", "est_se") %in% names(p1)))
  expect_error(
    predict(fit_dg, newdata = nd, se_fit = TRUE, type = "response"),
    regexp = "link scale"
  )

  # without se_fit = TRUE
  p1 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = FALSE,
    type = "link"
  )

  p2 <- predict(fit_dg,
    newdata = nd,
    re_form = NULL,
    re_form_iid = NULL,
    se_fit = FALSE,
    type = "response"
  )

  expect_equal(plogis(p1$est1), p2$est1)
  expect_equal(exp(p1$est2), p2$est2)
  expect_equal(plogis(p1$est1) * exp(p1$est2), p2$est)
})

test_that("Predicting without newdata matches predicting on the fitted data", {
  skip_on_cran()

  d <- pcod_2011
  d$off <- log(seq(0.5, 2, length.out = nrow(d)))
  d$year_scaled <- as.numeric(scale(d$year))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 20)
  check_no_newdata <- function(fit) {
    p1 <- predict(fit)
    p2 <- predict(fit, newdata = fit$data, offset = fit$offset)
    expect_equal(names(p1), names(p2))
    expect_equal(p1, p2, tolerance = 1e-6)
  }

  fit <- sdmTMB(density ~ s(depth), offset = "off", data = d, mesh = mesh,
    family = tweedie(), time = "year")
  check_no_newdata(fit)

  fit <- sdmTMB(log(density + 1) ~ 1, data = d, mesh = mesh,
    spatial_varying = ~ 0 + year_scaled, time = "year",
    spatiotemporal = "off")
  check_no_newdata(fit)

  fit <- sdmTMB(density ~ 1, data = subset(d, density > 0), mesh = mesh,
    family = lognormal_mix(), spatial = "off")
  check_no_newdata(fit)
})

test_that("Population predictions order columns like full predictions", {
  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011,
    mesh = pcod_mesh_2011, family = delta_gamma())
  nd <- pcod_2011[1:5, c("depth_scaled", "X", "Y")]
  p <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE)
  expect_identical(names(p), c(names(nd), "est", "est1", "est2", "est_se"))
})

test_that("Non-spatial predictions keep the user's coordinate columns", {
  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, spatial = "off",
    family = tweedie())
  nd <- pcod_2011[1:5, c("depth_scaled", "X", "Y")]
  p <- predict(fit, newdata = nd)
  expect_true(all(c("X", "Y") %in% names(p)))
  expect_equal(p$X, nd$X)
  p <- predict(fit, newdata = nd, type = "response")
  expect_true(all(c("X", "Y") %in% names(p)))
})
