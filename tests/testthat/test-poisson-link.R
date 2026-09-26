test_that("Poisson-link prediction issues from GitHub issues #389,", {
  skip_on_cran()
  dogfish$log_depth <- log(dogfish$depth)
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 30)
  fit_dpg <- sdmTMB(catch_weight ~ 0 + as.factor(year) + s(log_depth),
    family = delta_gamma(type = "poisson-link"),
    spatial = "on",
    mesh = mesh,
    data = dogfish,
    offset = log(dogfish$area_swept)
  )
  # Predicting on the fitted data uses the fitted offset. These values were
  # recorded at offset 0, so add the offset.
  offset <- log(dogfish$area_swept[1:5])
  p <- predict(fit_dpg, re_form = NA)
  p$est_est1est2 <- p$est1 + p$est2
  expect_equal(p$est, p$est_est1est2)
  expect_equal(p$est_est1est2[1:5],
    c(7.01028301430728, 2.34143881755487, 6.96979232578834, 6.99973559970208,
      7.03187981132451) + offset, tolerance = 1e-3)

  pp <- predict(fit_dpg, type = "response")
  expect_equal(pp$est[1:5],
    c(693.502617933497, 107.803298324919, 8414.54507288536, 5770.52404422525,
    6545.06096568627) * exp(offset), tolerance = 1e-3)
})

test_that("Poisson-link delta predictions include the offset", {
  skip_on_cran()
  mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 30)
  nd <- dogfish[1:20, ]
  for (backend in c("tmb", "rtmb")) {
    fit <- sdmTMB(catch_weight ~ 1, family = delta_gamma(type = "poisson-link"),
      mesh = mesh, data = dogfish, offset = log(dogfish$area_swept),
      control = sdmTMBcontrol(backend = backend))
    p0 <- predict(fit, newdata = nd, offset = rep(0, nrow(nd)))
    p1 <- predict(fit, newdata = nd, offset = rep(0.3, nrow(nd)))
    # The offset enters component 2, as in the expected catch
    # exp(offset + eta1 + eta2).
    expect_equal(p1$est1, p0$est1)
    expect_equal(p1$est2 - p0$est2, rep(0.3, nrow(nd)))
    expect_equal(p1$est - p0$est, rep(0.3, nrow(nd)))
    r0 <- predict(fit, newdata = nd, offset = rep(0, nrow(nd)), type = "response")
    r1 <- predict(fit, newdata = nd, offset = rep(0.3, nrow(nd)), type = "response")
    expect_equal(log(r1$est) - log(r0$est), rep(0.3, nrow(nd)))
    # Index with an offset scales by exp(offset).
    i0 <- get_index(fit, newdata = nd, offset = rep(0, nrow(nd)), area = 1)
    i1 <- get_index(fit, newdata = nd, offset = rep(0.3, nrow(nd)), area = 1)
    expect_equal(i1$est / i0$est, exp(0.3), tolerance = 1e-6)
  }
})
