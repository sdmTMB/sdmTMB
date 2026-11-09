test_that("RTMB Gaussian fit matches lm() through fit, prediction, simulation, and reload", {
  set.seed(2026)
  d <- data.frame(x = seq(-1, 1, length.out = 24L))
  d$y <- 0.4 + 0.7 * d$x + rnorm(nrow(d), sd = 0.6)
  rtmb <- sdmTMB(y ~ x, data = d, spatial = "off")
  m <- stats::lm(y ~ x, data = d)

  expect_true(rtmb$sd_report$pdHess)
  expect_equal(tidy(rtmb)$estimate, unname(coef(m)), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(rtmb)), as.numeric(logLik(m)),
    tolerance = 1e-6)
  actual <- predict(rtmb, newdata = d[1:4, ], se_fit = TRUE)
  expected <- predict(m, newdata = d[1:4, ], se.fit = TRUE)
  expect_equal(actual$est, unname(expected$fit), tolerance = 1e-5)
  # ML and OLS standard errors differ by sqrt((n - p) / n).
  expect_equal(actual$est_se,
    unname(expected$se.fit) * sqrt((nrow(d) - 2) / nrow(d)), tolerance = 1e-4)
  sims <- simulate(rtmb, nsim = 2L)
  expect_identical(dim(sims), c(nrow(d), 2L))
  expect_true(all(is.finite(sims)))
  expect_equal(as.numeric(simulate(rtmb, nsim = 1L, observation_error = FALSE)),
    as.numeric(predict(rtmb)$est))

  saved <- tempfile(fileext = ".rds")
  saveRDS(rtmb, saved)
  rebuilt <- reload_model(readRDS(saved))
  expect_equal(rebuilt$tmb_obj$fn(rebuilt$model$par), rtmb$model$objective)
})

test_that("RTMB rejects unsupported MakeADFun options and random parameters", {
  set.seed(1)
  d <- data.frame(y = rnorm(20), X = runif(20), Y = runif(20), year = 1L)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.2)
  fit <- sdmTMB(y ~ 1, data = d, mesh = mesh, spatial = "off",
    spatiotemporal = "off", time = "year")
  expect_error(make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
    fit$tmb_random, hessian = TRUE),
    "Additional MakeADFun options")
  expect_error(make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
    random = "ln_phi"), "can't integrate over")
})

test_that("RTMB fits survive update, extra optimization, and CV", {
  set.seed(2027)
  d <- data.frame(x = runif(60L), y = runif(60L), z = rnorm(60L))
  d$response <- 0.3 + 0.5 * d$z + sin(4 * d$x) + cos(4 * d$y) +
    rnorm(60L, sd = 0.2)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 20L, type = "kmeans")
  control <- sdmTMBcontrol(multiphase = FALSE)
  fit <- sdmTMB(response ~ z, data = d, mesh = mesh, control = control)
  objective <- fit$model$objective

  updated <- update(fit, response ~ 1)
  expect_false(isTRUE(all.equal(updated$model$objective, objective)))

  extra <- run_extra_optimization(fit, newton_loops = 1L)
  expect_equal(extra$model$objective, objective, tolerance = 1e-6)

  predict(fit, newdata = d[1:5, ])
  expect_equal(fit$tmb_obj$fn(fit$model$par), objective, ignore_attr = TRUE)

  d$fold <- rep(1:3, length.out = nrow(d))
  cv <- sdmTMB_cv(response ~ z, data = d, mesh = mesh, fold_ids = d$fold,
    control = control)
  expect_true(is.finite(cv$sum_loglik))
})

test_that("RTMB fits random effects with the default multiphase start", {
  set.seed(2028)
  d <- data.frame(z = rnorm(120L), g = factor(rep(1:12, 10L)))
  d$y <- 1 + d$z + rnorm(12L)[d$g] + rnorm(120L, sd = 0.5)
  multiphase <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off")
  single <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE))
  expect_true(multiphase$sd_report$pdHess)
  expect_equal(multiphase$model$par, single$model$par, tolerance = 1e-5)
})

test_that("update() keeps the fitted controls", {
  set.seed(3)
  d <- data.frame(x = seq(-1, 1, length.out = 20L))
  d$y <- 0.2 + 0.5 * d$x + rnorm(nrow(d), sd = 0.5)

  # a named control that is later removed; its settings are kept
  ctl <- sdmTMBcontrol(multiphase = FALSE, profile = TRUE)
  fit <- sdmTMB(y ~ x, data = d, spatial = "off", control = ctl)
  rm(ctl)
  new <- update(fit, y ~ 1)
  expect_false(new$control$multiphase)
  expect_identical(new$control$profile, fit$control$profile)
  call <- update(fit, y ~ 1, evaluate = FALSE)
  expect_false(eval(call)$control$multiphase)

  # an explicit control replaces the fitted one
  new <- update(fit, control = sdmTMBcontrol(multiphase = TRUE))
  expect_true(new$control$multiphase)
})
