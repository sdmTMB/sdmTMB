test_that("share_range mapping works with delta models", {
  skip_on_cran()

  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("off", "rw"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(NA, NA, 1, 1)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(TRUE, TRUE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(FALSE, TRUE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = list(FALSE, FALSE),
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = FALSE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("iid", "iid"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    spatiotemporal = list("off", "off"), time = "year",
    do_fit = FALSE, family = delta_gamma(),
    share_range = TRUE,
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(NA, NA, NA, NA)))

  # spatial models:
  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    do_fit = TRUE, family = tweedie()
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1)))

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    do_fit = TRUE, family = delta_gamma()
  )
  expect_identical(fit$tmb_map$ln_kappa, as.factor(c(1, 1, 2, 2)))
})

test_that("spatial field mapping/specification works with delta models", {
  skip_on_cran()

  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  fit1 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "on",
    family = delta_gamma()
  )
  s1 <- as.list(fit1$sd_report, "Estimate")
  s1$ln_tau_O

  fit2 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = list("on", "on"),
    family = delta_gamma()
  )
  s2 <- as.list(fit2$sd_report, "Estimate")
  s2$ln_tau_O
  expect_equal(s1, s2)

  fit3 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = delta_gamma()
  )
  s3 <- as.list(fit3$sd_report, "Estimate")
  s3$ln_tau_O

  fit4 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = list("off", "on"),
    family = delta_gamma()
  )
  s4 <- as.list(fit4$sd_report, "Estimate")
  s4$ln_tau_O
  expect_equal(s4$ln_tau_O[1], 0)
  expect_equal(s4$ln_tau_O[2], s2$ln_tau_O[2], tolerance = 0.01)


  pcod_spde2 <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  fit5 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde2, spatial = list("on", "on"),
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "iid")
  )

  s <- as.list(fit5$sd_report, "Estimate")
  expect_gt(abs(s$ln_tau_O[1]), 0)
  expect_gt(abs(s$ln_tau_O[2]), 0)
  expect_equal(s$ln_tau_E[1], 0)
  expect_gt(abs(s$ln_tau_E[2]), 0)
  expect_output(print(fit5), regexp = "Spatiotemporal model")
  expect_output(print(fit5), regexp = "Spatiotemporal IID SD")

  fit6 <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = list("off", "on"),
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "off")
  )
  s6 <- as.list(fit6$sd_report, "Estimate")
  expect_equal(s6$ln_tau_O[1], 0)
  expect_gt(abs(s6$ln_tau_O[2]), 0)
  expect_equal(s6$ln_tau_E[1], 0)
  expect_equal(s6$ln_tau_E[2], 0)
  expect_output(print(fit6), regexp = "Spatial model")
})

test_that("spatiotemporal field mapping/specification works with delta models", {
  skip_on_cran()

  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma()
  )
  s1 <- as.list(fit$sd_report, "Estimate")
  s1$ln_tau_E

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("iid", "iid")
  )
  s2 <- as.list(fit$sd_report, "Estimate")
  s2$ln_tau_E
  expect_equal(s1, s2)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("iid", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(abs(s$ln_tau_E[1]), 0)
  expect_equal(s$ln_tau_E[2], 0)
  expect_output(print(fit), regexp = "Spatiotemporal model")
  expect_output(print(fit), regexp = "Spatiotemporal IID SD")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "iid")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(abs(s$ln_tau_E[2]), 0)
  expect_equal(s$ln_tau_E[1], 0)
  print(fit)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_equal(s$ln_tau_E[1], 0)
  expect_equal(s$ln_tau_E[2], 0)
  print(fit)
  expect_output(print(fit), regexp = "Model fit")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("ar1", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(abs(s$ar1_phi[1]), 0)
  expect_equal(s$ar1_phi[2], 0)
  expect_identical(fit$tmb_map$ar1_phi, as.factor(c(1, NA)))
  tidy(fit, "ran_pars", model = 1)
  tidy(fit, "ran_pars", model = 2)
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("off", "ar1")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(abs(s$ar1_phi[2]), 0)
  expect_equal(s$ar1_phi[1], 0)
  expect_identical(fit$tmb_map$ar1_phi, as.factor(c(NA, 1)))
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("rw", "off")
  )
  s <- as.list(fit$sd_report, "Estimate")
  expect_gt(abs(s$ln_tau_E[1]), 0)
  expect_equal(s$ln_tau_E[2], 0)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = list("ar1", "ar1")
  )
  s1 <- as.list(fit$sd_report, "Estimate")
  print(fit)
  expect_output(print(fit), regexp = "rho")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    time = "year", family = delta_gamma(),
    spatiotemporal = "ar1"
  )
  s2 <- as.list(fit$sd_report, "Estimate")
  expect_equal(s1, s2, tolerance = 1e-4)
})

test_that("delta models work with different main effects", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)

  fit <- sdmTMB(
    formula = list(
      density ~ depth_scaled,
      density ~ depth_scaled + depth_scaled2
    ),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  s <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")
  expect_true("b_j2" %in% names(s))
  expect_true(all(!is.na(se$b_j2)))
  expect_true(all(!is.na(se$b_j)))
  tidy(fit)
  fit
  p <- predict(fit)

  # still works?
  fit <- sdmTMB(
    formula = density ~ depth_scaled + depth_scaled2,
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  fit
  p <- predict(fit)

  # one smoother works
  fit <- sdmTMB(
    formula = density ~ s(depth_scaled),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )

  # should be same:
  fit2 <- sdmTMB(
    formula = list(density ~ s(depth_scaled), density ~ s(depth_scaled)),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
  expect_equal(fit$sd_report, fit2$sd_report)

  # diff smoothers throw an error for now
  expect_error(
    {
      fit <- sdmTMB(
        formula = list(density ~ s(depth_scaled), density ~ s(year, k = 3)),
        data = pcod_2011,
        mesh = mesh,
        spatial = "off",
        family = delta_gamma()
      )
    },
    regexp = "smooth"
  )

  # diff random intercepts throw an error for now
  # expect_error(
  #   {
  #     pcod_2011$fyear <- as.factor(pcod_2011$year)
  #     fit <- sdmTMB(
  #       formula = list(density ~ 1 + (1 | fyear), density ~ 1),
  #       data = pcod_2011,
  #       mesh = mesh,
  #       spatial = "off",
  #       family = delta_gamma()
  #     )
  #   },
  #   regexp = "random intercepts"
  # )

  # OK:
  pcod_2011$fyear <- as.factor(pcod_2011$year)
  fit <- sdmTMB(
    formula = list(density ~ 1 + (1 | fyear), density ~ 1 + (1 | fyear)),
    data = pcod_2011,
    mesh = mesh,
    spatial = "off",
    family = delta_gamma()
  )
})

test_that("Offset works with delta models", {
  skip_on_cran()

  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  pcod_pos <- subset(pcod, density > 0)

  fit1 <- sdmTMB(present ~ 1,
    data = pcod, spatial = "off",
    family = binomial()
  )

  fit2 <- sdmTMB(density ~ 1,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  fit2_off <- sdmTMB(density ~ 1,
    offset = pcod_pos$offset,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  # error thrown if offset doesn't match data length:
  expect_error(
    {
      fit2_off_wrong <- sdmTMB(density ~ 1,
        offset = pcod$offset,
        data = pcod_pos, spatial = "off",
        family = Gamma(link = "log")
      )
    },
    regexp = "length"
  )

  fit_dg <- sdmTMB(density ~ 1,
    data = pcod, spatial = "off",
    family = delta_gamma()
  )

  fit_dg_off <- sdmTMB(density ~ 1,
    offset = pcod$offset,
    data = pcod, spatial = "off",
    family = delta_gamma()
  )

  pcod$offset2 <- log(1)
  fit_dg_off0 <- sdmTMB(density ~ 1,
    offset = pcod$offset2,
    data = pcod, spatial = "off",
    family = delta_gamma()
  )

  # intercept only models so order not an issue
  b1 <- tidy(fit1)$estimate[1]
  b_dg1 <- tidy(fit_dg)$estimate[1]
  b_dg1_offset <- tidy(fit_dg_off)$estimate[1]

  b2 <- tidy(fit2)$estimate[1]
  b2_offset <- tidy(fit2_off)$estimate[1]
  dg2 <- tidy(fit_dg, model = 2)$estimate[1]
  dg2_offset <- tidy(fit_dg_off, model = 2)$estimate[1]
  dg2_offset0 <- tidy(fit_dg_off0, model = 2)$estimate[1]

  # the offset is doing something for pos part of delta model
  expect_false(((dg2_offset - dg2) == 0))
  # binomial and delta model 1 without offset are same
  expect_equal(b_dg1, b1, tolerance = 1e-5)
  # offset doesn't affect binomial part of delta-Gamma
  expect_equal(b_dg1, b_dg1_offset, tolerance = 1e-5)
  # gamma on pos only and delta model 2 without offset are same
  expect_equal(dg2, b2, tolerance = 1e-5)
  # offset in Gamma part same in delta gamma as separate model:
  expect_equal(dg2_offset, b2_offset, tolerance = 1e-5)

  # the offset of 0 is same as no offset
  expect_equal(dg2_offset0, dg2, tolerance = 1e-8)
})

test_that("test that delta beta model works", {
  skip_on_cran()

  set.seed(1)
  y01 <- stats::rbinom(1000, 1, 0.5)
  npos <- sum(y01 == 1)
  y <- y01
  y[y == 1] <- stats::rbeta(npos, 3, 4)
  ypos <- y[y01 == 1]
  dat <- data.frame(y = y)

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    family = delta_beta(),
    control = sdmTMBcontrol(newton_loops = 1L)
  )

  d01 <- data.frame(y = y01)
  dpos <- data.frame(y = ypos)

  m1 <- glmmTMB::glmmTMB(y ~ 1, data = d01, family = binomial())
  s1 <- tidy(fit, effects = c("fixed"))

  m2 <- glmmTMB::glmmTMB(y ~ 1, data = dpos, family = glmmTMB::beta_family())
  s2 <- tidy(fit, effects = c("fixed"), model = 2)

  expect_equal(s1$estimate, m1$fit$par[[1]], tolerance = 1e-4)
  expect_equal(s2$estimate, m2$fit$par[[1]], tolerance = 1e-4)

  p <- predict(fit)

  p1 <- predict(m1)
  p2 <- predict(m2)
  # p <- predict(fit, type = "response")
  # glmmTMB_est <- stats::plogis(p1)[1] * stats::plogis(p2)[1]
  # expect_equal(p$est[1], glmmTMB_est, tolerance = 1e-4)

  r <- residuals(fit)
  qqnorm(r)
  qqline(r)

  set.seed(1)
  s <- simulate(fit)
  expect_gte(min(s), 0)
  expect_lte(min(s), 1)
})

test_that("one spatial off in a delta model works", {
  skip_on_cran()

  mesh0 <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  m0 <- sdmTMB(
    density ~ 1,
    mesh = mesh0,
    data = pcod,
    spatial = list("off", "on"), #<
    spatiotemporal = list("off", "off"),
    silent = TRUE,
    time = "year",
    family = delta_gamma()
  )

  # m0$tmb_map$omega_s
  # m0$tmb_map$ln_tau_O
  # m0$tmb_map$ln_kappa
  # m0$tmb_data$include_spatial
  # m0$tmb_data$spatial_only
  # m0$tmb_map$ln_tau_E
  # m0$tmb_map$epsilon_st
  # m0$tmb_params$ln_tau_E
  # m0$tmb_params$epsilon_re
  # m0$tmb_params$ln_tau_O

  pos <- subset(pcod, density > 0)
  mesh2 <- sdmTMB::make_mesh(pos, xy_cols = c("X", "Y"), mesh = mesh0$mesh)
  m2 <- sdmTMB(
    density ~ 1,
    mesh = mesh2,
    data = pos,
    spatial = "on", # <-
    spatiotemporal = "off",
    silent = TRUE,
    time = "year",
    family = Gamma(link = "log")
  )

  # m0$tmb_obj$report()$sigma_O
  # s0 <- as.list(m0$sd_report, what = "Estimate", report = TRUE)
  # s0$sigma_E
  # s0$sigma_O
  # s0$range
  #
  # s2 <- as.list(m2$sd_report, what = "Estimate", report = TRUE)
  # s2$sigma_E
  # s2$range

  t0 <- tidy(m0, "ran_pars", model = 2)
  t2 <- tidy(m2, "ran_pars")

  expect_equal(t0[,c("term","estimate","std.error","conf.low","conf.high")],
               t2[,c("term","estimate","std.error","conf.low","conf.high")], tolerance = 0.01)

  # ---------------------
  # with sigma_E

  m0 <- sdmTMB(
    density ~ 1,
    mesh = mesh0,
    data = pcod,
    spatial = list("off", "off"), #<
    spatiotemporal = list("off", "iid"), #<
    share_range = FALSE,
    silent = TRUE,
    time = "year",
    control = sdmTMBcontrol(newton_loops = 0L),
    family = delta_gamma()
  )
  m2 <- sdmTMB(
    density ~ 1,
    mesh = mesh2,
    data = pos,
    spatial = "off", # <-
    spatiotemporal = "iid",
    share_range = FALSE,
    silent = TRUE,
    time = "year",
    control = sdmTMBcontrol(newton_loops = 0L),
    family = Gamma(link = "log")
  )
  t0 <- tidy(m0, "ran_pars", model = 2)
  t2 <- tidy(m2, "ran_pars")
  expect_equal(t0$estimate, t2$estimate, tolerance = 0.01)
})
