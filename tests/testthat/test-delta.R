pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
pcod_pos <- subset(pcod, density > 0)
pcod_spde_pos <- make_mesh(pcod_pos, c("X", "Y"), mesh = pcod_spde$mesh)

test_that("Delta-Gamma family fits", {
  skip_on_cran()

  fit_dg <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    time = "year", family = delta_gamma()
  )
  fit_dg$sd_report
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))

  expect_equal(
    round(tidy(fit_dg, "ran_pars", model = 1)$estimate, 3),
    c(39.334, 2.289, 0.808)
  )
  expect_equal(
    round(tidy(fit_dg, "ran_pars", model = 2)$estimate, 3),
    c(16.224, 0.992, 0.656, 1.426)
  )

  p <- predict(fit_dg, newdata = nd, return_tmb_object = TRUE)
  ind_dg <- get_index(p, bias_correct = FALSE)

  # check
  fit_bin <- sdmTMB(present ~ 1,
    data = pcod, mesh = pcod_spde,
    time = "year", family = binomial(),
    control = sdmTMBcontrol(newton_loops = 1)
  )
  fit_gamma <- sdmTMB(density ~ 1,
    data = pcod_pos, mesh = pcod_spde_pos,
    time = "year", family = Gamma(link = "log"),
    control = sdmTMBcontrol(newton_loops = 1)
  )
  sr_bin <- as.list(fit_bin$sd_report, "Estimate")
  sr_gamma <- as.list(fit_gamma$sd_report, "Estimate")
  sr_dg <- as.list(fit_dg$sd_report, "Estimate")
  expect_equal(sr_bin$b_j[1], sr_dg$b_j[1], tolerance = 1e-4)
  expect_equal(sr_gamma$b_j[1], sr_dg$b_j2[1], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_phi, sr_dg$ln_phi[1], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_tau_O, sr_dg$ln_tau_O[2], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_tau_E, sr_dg$ln_tau_E[2], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_kappa[1,1], sr_dg$ln_kappa[1,2], tolerance = 1e-4)
  expect_equal(sr_bin$ln_kappa[1,1], sr_dg$ln_kappa[1,1], tolerance = 1e-4)
})

test_that("Delta-lognormal family fits", {
  skip_on_cran()

  fit_dln <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = delta_lognormal()
  )
  s <- as.list(fit_dln$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
})

test_that("delta_gamma() Poisson-link family fits", {
  skip_on_cran()

  fit_plg <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = delta_gamma(type = "poisson-link")
  )
  fit_plg$sd_report
  s <- as.list(fit_plg$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
  expect_equal(fit_plg$family[[1]]$link, "log")

  # p <- predict(fit_plg, newdata = qcs_grid, type = "response")
  # p <- predict(fit_plg, newdata = pcod, type = "response")
  # expect_error(p <- predict(fit_plg, newdata = NULL, type = "response"))
})

test_that("delta_lognormal() Poisson-link family fits", {
  skip_on_cran()

  fit_plg <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = delta_lognormal(type = "poisson-link")
  )
  fit_plg$sd_report
  s <- as.list(fit_plg$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
  expect_equal(fit_plg$family[[1]]$link, "log")
})

test_that("delta_truncated_nbinom2 family fits", {
  skip_on_cran()

  pcod$count <- round(pcod$density)
  fit_dtnb2 <- sdmTMB(count ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    family = delta_truncated_nbinom2()
  )
  s <- as.list(fit_dtnb2$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
  fit_dtnb2$sd_report
})


test_that("delta_truncated_nbinom1 family fits", {
  skip_on_cran()

  pcod$count <- round(pcod$density)
  fit_dtnb1 <- sdmTMB(count ~ 1,
    data = pcod, mesh = pcod_spde, spatial = "off",
    family = delta_truncated_nbinom1()
  )
  s <- as.list(fit_dtnb1$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
  fit_dtnb1$sd_report
})

test_that("Anisotropy with delta model", {
  skip_on_cran()

  suppressWarnings({
    fit_dg <- sdmTMB(density ~ 1,
      data = pcod, mesh = pcod_spde,
      time = "year", family = delta_gamma(),
      anisotropy = TRUE,
      control = sdmTMBcontrol(
        newton_loops = 1,
        map = list(ln_H_input = factor(c(1L, 2L, 3L, 4L)))
      )
    )
  })

  fit_bin <- sdmTMB(present ~ 1,
    data = pcod, mesh = pcod_spde,
    time = "year", family = binomial(),
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1)
  )
  p1 <- plot_anisotropy2(fit_bin)
  p2 <- plot_anisotropy2(fit_dg)
  expect_equal(p2, p1, tolerance = 1e-6)

  fit_gamma <- sdmTMB(density ~ 1,
    data = pcod_pos, mesh = pcod_spde_pos,
    time = "year", family = Gamma(link = "log"),
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1)
  )

  p3 <- plot_anisotropy2(fit_gamma)
  p4 <- plot_anisotropy2(fit_dg, model = 2)

  ## not sure why this isn't working
  # expect_equal(p4, p3, tolerance = 1e-6)

  fit_bin$sd_report
  fit_gamma$sd_report
  fit_dg$sd_report

  sr_bin <- as.list(fit_bin$sd_report, "Estimate")
  sr_gamma <- as.list(fit_gamma$sd_report, "Estimate")
  sr_dg <- as.list(fit_dg$sd_report, "Estimate")

  expect_equal(sr_bin$ln_H_input[1,1], sr_dg$ln_H_input[1,1], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_H_input[1,1], sr_dg$ln_H_input[1,2], tolerance = 1e-4)

  # all estimates still match?
  expect_equal(sr_bin$b_j[1], sr_dg$b_j[1], tolerance = 1e-4)
  expect_equal(sr_bin$ln_kappa[1,1], sr_dg$ln_kappa[1,1], tolerance = 1e-4)
  expect_equal(sr_gamma$b_j[1], sr_dg$b_j2[1], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_phi, sr_dg$ln_phi[1], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_tau_O, sr_dg$ln_tau_O[2], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_tau_E, sr_dg$ln_tau_E[2], tolerance = 1e-4)
  expect_equal(sr_gamma$ln_kappa[1,1], sr_dg$ln_kappa[1,2], tolerance = 1e-4)
})

test_that("Delta-Gengamma family fits", {
  skip_on_cran()

  fit_dgg <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = 'off',
    family = delta_gengamma(link1 = 'logit', link2 = 'log', type = 'standard')
  )
  fit_dgg$sd_report

  expect_equal(
    round(tidy(fit_dgg, "fixed", model = 1)$estimate, 3),
    -0.152
  )
  expect_equal(
    round(tidy(fit_dgg, "fixed", model = 2)$estimate, 3),
    c(4.418)
  )

  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  p <- predict(fit_dgg, newdata = nd, return_tmb_object = TRUE)

  # check
  fit_bin <- sdmTMB(present ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = binomial()
  )
  fit_gengamma <- sdmTMB(density ~ 1,
    data = pcod_pos, mesh = pcod_spde_pos,
    spatial = "off",
    family = gengamma(link = "log")
  )
  sr_bin <- as.list(fit_bin$sd_report, "Estimate")
  sr_gengamma <- as.list(fit_gengamma$sd_report, "Estimate")
  sr_dgg <- as.list(fit_dgg$sd_report, "Estimate")
  expect_equal(sr_bin$b_j[1], sr_dgg$b_j[1], tolerance = 1e-4)
  expect_equal(sr_gengamma$b_j[1], sr_dgg$b_j2[1], tolerance = 1e-4)
  expect_equal(sr_gengamma$ln_phi, sr_dgg$ln_phi[1], tolerance = 1e-4)
  expect_equal(sr_gengamma$ln_tau_O, sr_dgg$ln_tau_O[2], tolerance = 1e-4)
  expect_equal(sr_gengamma$ln_tau_E, sr_dgg$ln_tau_E[2], tolerance = 1e-4)
  expect_equal(sr_gengamma$ln_kappa[1,1], sr_dgg$ln_kappa[1,2], tolerance = 1e-4)
  expect_equal(sr_bin$ln_kappa[1,1], sr_dgg$ln_kappa[1,1], tolerance = 1e-4)
})

test_that("delta_gengamma() Poisson-link family fits", {
  skip_on_cran()

  fit_pgg <- sdmTMB(density ~ 1,
    data = pcod, mesh = pcod_spde,
    spatial = "off",
    family = delta_gengamma(link1 = 'log', link2 = 'log', type = "poisson-link")
  )
  fit_pgg$sd_report
  s <- as.list(fit_pgg$sd_report, "Std. Error")
  expect_true(sum(is.na(s$b_j)) == 0L)
})
