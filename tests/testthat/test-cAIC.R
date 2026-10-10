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
    fit$tmb_random, backend_sdmTMB(fit), profile = NULL)
  par <- obj$env$last.par
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
