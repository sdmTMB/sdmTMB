test_that("cAIC and EDF work", {
  skip_on_cran()
  skip_on_ci()
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 15)
  suppressMessages(
    fit <- sdmTMB(catch_weight ~ s(log(depth)),
      time_varying = ~1,
      time_varying_type = "ar1",
      time = "year",
      spatiotemporal = "off",
      mesh = mesh,
      family = tweedie(),
      data = dogfish,
      offset = log(dogfish$area_swept)
    )
  )
  expect_equal(AIC(fit), 12192.9613, tolerance = 1e-4)
  expect_equal(cAIC(fit), 12071.4289, tolerance = 1e-4)
  edf <- cAIC(fit, what = "EDF")
  expect_equal(sum(edf), 54.3870623, tolerance = 1e-4)
})

test_that("cAIC sparse EDF matches dense calculation in a delta model", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  fit <- sdmTMB(density ~ s(depth),
    data = pcod, mesh = mesh, family = delta_gamma(),
    time = "year", spatiotemporal = "iid"
  )
  obj <- fit$tmb_obj
  tmb_data <- fit$tmb_data
  tmb_data$weights_i[] <- 0
  obj_new <- make_sdmTMB_adfun(tmb_data, fit$parlist, fit$tmb_map,
    fit$tmb_random, profile = NULL)
  par <- fit$last.par.best
  H <- as.matrix(obj$env$spHess(par, random = TRUE))
  H_new <- as.matrix(obj_new$env$spHess(par, random = TRUE))
  negEDF <- diag(solve(H, H_new))
  group <- names(par[obj$env$random])
  edf <- cAIC(fit, what = "EDF")
  expect_equal(sum(edf), length(negEDF) - sum(negEDF), tolerance = 1e-8)
  expect_equal(as.numeric(edf[c("epsilon_st", "omega_s")]),
    c(
      sum(group == "epsilon_st") - sum(negEDF[group == "epsilon_st"]),
      sum(group == "omega_s") - sum(negEDF[group == "omega_s"])
    ), tolerance = 1e-8)
  q <- length(negEDF)
  p <- length(fit$model$par)
  cnll <- obj$env$f(par) - obj_new$env$f(par)
  expect_equal(cAIC(fit), 2 * cnll + 2 * (p + q) - 2 * sum(negEDF), tolerance = 1e-8)
})

test_that("cAIC compares multiple models like AIC", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  fit1 <- sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = tweedie())
  fit2 <- sdmTMB(density ~ depth_scaled, data = pcod, mesh = mesh,
    family = tweedie())
  fit_glm <- sdmTMB(density ~ depth_scaled, data = pcod, spatial = "off",
    family = tweedie())

  out <- cAIC(fit1, fit2, fit_glm)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("df", "cAIC"))
  expect_identical(row.names(out), c("fit1", "fit2", "fit_glm"))
  expect_equal(out$cAIC[1:2], c(cAIC(fit1), cAIC(fit2)))
  expect_equal(out$df[1:2], c(
    length(fit1$model$par) + sum(cAIC(fit1, what = "EDF")),
    length(fit2$model$par) + sum(cAIC(fit2, what = "EDF"))
  ))
  expect_equal(out$cAIC[3], AIC(fit_glm))
  expect_equal(out$df[3], length(fit_glm$model$par))

  expect_error(cAIC(fit1, fit2, what = "EDF"), "single model")
  expect_error(cAIC(fit1, "EDF"), "must now be named")
  expect_error(cAIC(fit1, list()), "fitted")
  fit_sub <- update(fit1, data = pcod[-1, ], mesh = make_mesh(pcod[-1, ], c("X", "Y"), cutoff = 25))
  expect_warning(cAIC(fit1, fit_sub), "same number of observations")
})

test_that("cAIC uses the mode and counts profiled fixed effects", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(density ~ depth_scaled + depth_scaled2, data = pcod_2011,
    mesh = mesh, family = tweedie())
  fit_prof <- update(fit, control = sdmTMBcontrol(profile = TRUE))
  out <- cAIC(fit, fit_prof)
  # same model, so equal up to optimizer noise (leaving out the 3 profiled
  # coefficients would lower df by 3 and cAIC by 6)
  expect_lt(abs(out$df[1] - out$df[2]), 0.5)
  expect_lt(abs(out$cAIC[1] - out$cAIC[2]), 1)

  # evaluating the objective elsewhere moves `last.par`, not the mode
  caic <- cAIC(fit)
  fit$tmb_obj$fn(fit$model$par + 0.1)
  expect_equal(cAIC(fit), caic)
})

test_that("cAIC counts the cutpoint of profile_breakpt() output", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(x = runif(300, -2, 2), g = factor(rep(1:15, each = 20)))
  d$y <- 1 + 1.5 * pmin(d$x, 0.4) + rnorm(15, sd = 0.5)[d$g] +
    rnorm(300, sd = 0.5)
  p <- profile_breakpt(y ~ breakpt(x) + (1 | g), data = d, spatial = "off")
  fit_lin <- sdmTMB(y ~ x + (1 | g), data = d, spatial = "off")

  expect_equal(cAIC(p), cAIC(p$fit) + 2)
  expect_equal(cAIC(p, what = "EDF"), cAIC(p$fit, what = "EDF"))
  out <- cAIC(p, fit_lin)
  expect_identical(row.names(out), c("p", "fit_lin"))
  expect_equal(out$cAIC, c(cAIC(p), cAIC(fit_lin)))
  expect_equal(out$df[1], length(p$fit$model$par) + 1 + sum(cAIC(p$fit, what = "EDF")))
  # either order dispatches to the same method
  expect_equal(cAIC(fit_lin, p)$cAIC, rev(out$cAIC))
  expect_equal(nobs(p), nobs(p$fit))

  p_nofit <- profile_breakpt(y ~ breakpt(x) + (1 | g), data = d,
    spatial = "off", final = FALSE)
  expect_error(cAIC(p_nofit), "final = TRUE")
})
