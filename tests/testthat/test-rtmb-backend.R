test_that("RTMB Gaussian fit agrees with TMB through fit, prediction, simulation, and reload", {
  set.seed(2026)
  d <- data.frame(x = seq(-1, 1, length.out = 24L))
  d$y <- 0.4 + 0.7 * d$x + rnorm(nrow(d), sd = 0.6)
  tmb <- sdmTMB(y ~ x, data = d, spatial = "off",
    control = sdmTMBcontrol(multiphase = FALSE, backend = "tmb"))
  rtmb <- sdmTMB(y ~ x, data = d, spatial = "off",
    control = sdmTMBcontrol(backend = "rtmb"))

  expect_equal(rtmb$model$par, tmb$model$par, tolerance = 1e-6)
  expect_equal(rtmb$model$objective, tmb$model$objective, tolerance = 1e-6)
  expect_true(rtmb$sd_report$pdHess)
  expect_sdreport_matches(rtmb, tmb)
  expect_equal(rtmb$tmb_obj$gr(rtmb$model$par),
    tmb$tmb_obj$gr(tmb$model$par), tolerance = 1e-5)
  expected <- predict(tmb, newdata = d[1:4, ], se_fit = TRUE)
  actual <- predict(rtmb, newdata = d[1:4, ], se_fit = TRUE)
  expect_equal(actual$est, expected$est, tolerance = 1e-5)
  expect_equal(actual$est_se, expected$est_se, tolerance = 1e-5)
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
    spatiotemporal = "off", time = "year",
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_error(make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
    fit$tmb_random, backend = "rtmb", ADreport = TRUE),
    "Additional MakeADFun options")
  expect_error(make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
    random = "ln_phi", backend = "rtmb"), "can't integrate over")
})

test_that("RTMB fits survive update, extra optimization, and CV", {
  set.seed(2027)
  d <- data.frame(x = runif(60L), y = runif(60L), z = rnorm(60L))
  d$response <- 0.3 + 0.5 * d$z + sin(4 * d$x) + cos(4 * d$y) +
    rnorm(60L, sd = 0.2)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 12L, type = "kmeans")
  control <- sdmTMBcontrol(backend = "rtmb", multiphase = FALSE)
  fit <- sdmTMB(response ~ z, data = d, mesh = mesh, control = control)
  objective <- fit$model$objective

  updated <- update(fit, response ~ 1)
  expect_identical(updated$backend, "rtmb")
  expect_identical(attr(updated$tmb_obj, "sdmTMB_backend"), "rtmb")
  expect_false(isTRUE(all.equal(updated$model$objective, objective)))

  extra <- run_extra_optimization(fit, newton_loops = 1L)
  expect_identical(extra$backend, "rtmb")
  expect_equal(extra$model$objective, objective, tolerance = 1e-6)

  predict(fit, newdata = d[1:5, ])
  expect_equal(fit$tmb_obj$fn(fit$model$par), objective, ignore_attr = TRUE)

  d$fold <- rep(1:3, length.out = nrow(d))
  cv <- sdmTMB_cv(response ~ z, data = d, mesh = mesh, fold_ids = d$fold,
    control = control)
  expect_true(all(vapply(cv$models, function(x) x$backend, "") == "rtmb"))
  expect_true(is.finite(cv$sum_loglik))
})

test_that("RTMB fits random effects with the default multiphase start", {
  set.seed(2028)
  d <- data.frame(z = rnorm(120L), g = factor(rep(1:12, 10L)))
  d$y <- 1 + d$z + rnorm(12L)[d$g] + rnorm(120L, sd = 0.5)
  tmb <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off",
    control = sdmTMBcontrol(backend = "tmb"))
  rtmb <- sdmTMB(y ~ z + (1 | g), data = d, spatial = "off",
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_equal(rtmb$model$par, tmb$model$par, tolerance = 1e-5)
})

test_that("update() keeps the fitted backend and controls", {
  set.seed(3)
  d <- data.frame(x = seq(-1, 1, length.out = 20L))
  d$y <- 0.2 + 0.5 * d$x + rnorm(nrow(d), sd = 0.5)
  old_options <- options(sdmTMB.backend = NULL)
  on.exit(options(old_options), add = TRUE)
  for (backend in c("rtmb", "tmb")) {
    other <- setdiff(c("tmb", "rtmb"), backend)
    options(sdmTMB.backend = backend)
    fit <- sdmTMB(y ~ x, data = d, spatial = "off")
    expect_identical(fit$backend, backend)
    options(sdmTMB.backend = other)
    expect_identical(update(fit, y ~ 1)$backend, backend)
    call <- update(fit, y ~ 1, evaluate = FALSE)
    expect_identical(eval(call)$backend, backend)
    # an explicit control replaces the fitted one, including its backend
    new <- update(fit, control = sdmTMBcontrol(backend = other))
    expect_identical(new$backend, other)
  }

  # a named control that is later removed; other settings are kept
  ctl <- sdmTMBcontrol(multiphase = FALSE, profile = TRUE, backend = "rtmb")
  fit <- sdmTMB(y ~ x, data = d, spatial = "off", control = ctl)
  rm(ctl)
  options(sdmTMB.backend = "tmb")
  new <- update(fit, y ~ 1)
  expect_identical(new$backend, "rtmb")
  expect_false(new$control$multiphase)
  expect_identical(new$control$profile, fit$control$profile)

  # a supplied legacy control list without a backend keeps the fitted one
  legacy <- fit$control
  legacy$backend <- NULL
  expect_identical(update(fit, control = legacy)$backend, "rtmb")

  # old objects without backend or control metadata fall back to TMB
  options(sdmTMB.backend = "rtmb")
  old <- fit
  old$backend <- NULL
  old$control <- NULL
  expect_identical(update(old, evaluate = FALSE)$control$backend, "tmb")
  old$control <- fit$control
  old$control$backend <- NULL
  expect_identical(update(old, y ~ 1)$backend, "tmb")
})
