test_that("visreg works", {
  skip_on_cran()
  skip_if_not_installed("visreg")

  # first test against glmmTMB outputs
  pcod_2011$fyear <- as.factor(pcod_2011$year)

  # for binomial
  fit_glm <- glmmTMB::glmmTMB(
    present ~ poly(depth_scaled, 2) + fyear,
    data = pcod_2011,
    family = binomial("logit")
  )
  v1 <- visreg::visreg(fit_glm, xvar = "depth_scaled", nn = 10, plot = FALSE)
  v1r <- visreg::visreg(
    fit_glm, xvar = "depth_scaled", nn = 10,
    scale = "response", plot = FALSE
  )

  fit_sdm <- sdmTMB(
    present ~ poly(depth_scaled, 2) + fyear,
    data = pcod_2011, spatial = "off",
    family = binomial("logit")
  )
  v2 <- visreg::visreg(fit_sdm, xvar = "depth_scaled", nn = 10, plot = FALSE)
  v2r <- visreg::visreg(
    fit_sdm, xvar = "depth_scaled", nn = 10,
    scale = "response", plot = FALSE
  )

  # check that these models are matching
  expect_equal(AIC(fit_glm), AIC(fit_sdm), tolerance = 1e-6)
  # check that output is visreg
  expect_identical(class(v2), "visreg")
  # check which elements match between glmm and sdmTMB versions
  expect_equal(v2$fit, v1$fit, tolerance = 1e-6)
  expect_equal(v2r$fit, v1r$fit, tolerance = 1e-6)
  # residuals are calculated for all original data
  expect_identical(nrow(v2r$res), nrow(v1r$res))
  # The residuals are calculated differently and are larger for sdmTMB.

  # for gamma
  pcod_pos <- dplyr::filter(pcod_2011, present > 0)

  fit_glm <- glmmTMB::glmmTMB(
    density ~ poly(depth_scaled, 2) + fyear,
    data = pcod_pos,
    family = Gamma("log")
  )
  v1 <- visreg::visreg(fit_glm, xvar = "depth_scaled", nn = 10, plot = FALSE)
  v1r <- visreg::visreg(
    fit_glm, xvar = "depth_scaled", nn = 10,
    scale = "response", plot = FALSE
  )
  ## Plot customization requires ggplot2.
  # (g1 <- visreg::visreg(fit_glm, xvar = "depth_scaled", nn = 10,
  #                       scale = "response", rug = FALSE))

  fit_sdm <- sdmTMB(
    density ~ poly(depth_scaled, 2) + fyear,
    data = pcod_pos, spatial = "off",
    family = Gamma("log")
  )
  v2 <- visreg::visreg(fit_sdm, xvar = "depth_scaled", nn = 10, plot = FALSE)
  v2r <- visreg::visreg(
    fit_sdm, xvar = "depth_scaled", nn = 10,
    scale = "response", plot = FALSE
  )
  # (g2 <- visreg::visreg(fit_sdm, xvar = "depth_scaled", nn = 10,
  #                       scale = "response", rug = FALSE))
  # check that these models are matching
  expect_equal(AIC(fit_glm), AIC(fit_sdm), tolerance = 1e-6)
  # check that output is visreg
  expect_identical(class(v2), "visreg")
  # check which elements match between glmm and sdmTMB versions
  expect_equal(v2r$fit, v1r$fit, tolerance = 1e-6)
  # residuals are calculated for all original data
  expect_identical(nrow(v2r$res), nrow(v1r$res))
  # The residuals are calculated differently and are much larger for glmmTMB.

  # with smoother, tweedie, and reml
  fit <- sdmTMB(
    density ~ s(depth_scaled) + fyear,
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    reml = TRUE,
    family = tweedie()
  )
  visreg::visreg(fit, xvar = "depth_scaled", nn = 10)
  visreg::visreg(fit, xvar = "fyear", nn = 10)
  visreg::visreg(fit, xvar = "depth_scaled", scale = "response", nn = 10)
  visreg::visreg2d(
    fit, xvar = "fyear", yvar = "depth_scaled", nn = 10, plot = FALSE
  )
  v <- visreg::visreg(fit, xvar = "depth_scaled", nn = 10, plot = FALSE)
  expect_identical(class(v), "visreg")


  # works with a spatiotemporal model with time that isn't also a fixed effect
  pcod$fyear <- as.factor(pcod$year)
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(
    density ~ 1 + poly(depth_scaled, 2) + (1|fyear),
    data = pcod, mesh = mesh,
    time = "year",
    spatial = "on",
    spatiotemporal = TRUE,
    family = tweedie()
  )

  v <- visreg::visreg(fit, xvar = "depth_scaled", nn = 10, plot = FALSE)
  expect_identical(class(v), "visreg")
  #visreg::visreg(fit, xvar = "fyear", nn = 10)
  visreg::visreg(fit, xvar = "depth_scaled", scale = "response", nn = 10)
  visreg::visreg2d(
    fit, xvar = "fyear", yvar = "depth_scaled", plot = FALSE
  )

  # works with RW models
  fit_rw <- sdmTMB(
    density ~ 1 + poly(depth_scaled, 2),
    data = pcod_2011, mesh = pcod_mesh_2011,
    time = "year",
    spatial = "off",
    spatiotemporal = "rw",
    family = tweedie()
  )
  visreg::visreg(fit_rw, xvar = "depth_scaled", nn = 10)

  # works with ar1 models
  fit_ar1 <- sdmTMB(
    density ~ 1 + poly(depth_scaled, 2),
    data = pcod_2011, mesh = pcod_mesh_2011,
    time = "year",
    spatial = "off",
    spatiotemporal = "ar1",
    family = tweedie()
  )

  visreg::visreg(fit_ar1, xvar = "depth_scaled", nn = 10)

  # Delta model example:
  fit_dg <- sdmTMB(
    density ~ s(depth_scaled, year, k = 8),
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    family = delta_gamma()
  )
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 1, nn = 10)
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 2, nn = 10)
  v3 <- visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 1,
    scale = "response", nn = 10, plot = FALSE
  )

  fit_bin <- sdmTMB(
    present ~ s(depth_scaled, year, k = 8),
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    family = binomial()
  )

  v4 <- visreg::visreg(fit_bin,
    xvar = "depth_scaled",
    scale = "response", nn = 10, plot = FALSE
  )

  # outputs differ because response is still called density even for model 1 of a delta model
  # expect_equal(v3$fit, v4$fit, tolerance = 1e-6)

  # check that fit values match between binomial and model 1 of the delta
  fit_col <- grep("^visreg_?[Ff]it$", names(v3$fit), value = TRUE)
  expect_length(fit_col, 1L)
  expect_equal(v3$fit[[fit_col]], v4$fit[[fit_col]], tolerance = 1e-6)

  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response", nn = 10
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response", nn = 10
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 1, scale = "response", nn = 10
  )
  visreg2d_delta(fit_dg, xvar = "depth_scaled", yvar = "year", model = 1, nn = 10)
  visreg2d_delta(fit_dg, xvar = "depth_scaled", yvar = "year", model = 2, nn = 10)

  v <- visreg_delta(
    fit_dg, xvar = "depth_scaled", model = 1, nn = 10, plot = FALSE
  )
  expect_identical(class(v), "visreg")

  v <- visreg2d_delta(
    fit_dg, xvar = "depth_scaled", yvar = "year",
    model = 2, nn = 10, plot = FALSE
  )
  expect_identical(class(v), "visreg2d")
})

test_that("visreg works with extra time #330", {
  skip_on_cran()
  pcod_2011$fyear <- as.factor(pcod_2011$year)
  fit <- sdmTMB(
    density ~ 0 + depth_scaled,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    time_varying = ~ 1,
    data = pcod_2011,
    family = tweedie(),
    extra_time = c(2012, 2014, 2016)
  )
  set.seed(1)
  v <- visreg::visreg(fit, "depth_scaled", plot = FALSE)
  expect_s3_class(v, 'visreg')
})
