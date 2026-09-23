test_that("A logistic threshold model fits", {
  skip_on_cran()
  d <- subset(pcod, year >= 2011) # subset for speed
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + logistic(depth_scaled), data = d,
    mesh = pcod_spde, family = tweedie(link = "log"),
    time = "year")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("depth_scaled-s50" %in% tidy(m)$term)
  expect_true("depth_scaled-s95" %in% tidy(m)$term)
  expect_true("depth_scaled-smax" %in% tidy(m)$term)
  expect_equal(tidy(m)[,"estimate",drop=TRUE], c(1.555 , 1.655 , 1.718 , 1.138, -0.979, -0.937 , 1.760), tolerance = 1e-3)
})

test_that("A linear threshold model fits", {
  set.seed(41)
  d <- data.frame(x = seq(-2, 2, length.out = 300))
  d$y <- 2 + 1.5 * pmin(d$x, 0.35) + rnorm(nrow(d), sd = 0.25)
  m <- sdmTMB(y ~ 1 + breakpt(x), data = d,
    family = gaussian(), spatial = "off")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))

  expect_true("x-slope" %in% tidy(m)$term)
  expect_true("x-breakpt" %in% tidy(m)$term)
  expect_equal(tidy(m)[, "estimate", drop = TRUE], c(2, 1.5, 0.35),
    tolerance = 0.1)
})

test_that("A linear threshold *delta* model fits", {
  skip_on_cran()

  set.seed(1)
  predictor_dat <- data.frame(
    X = runif(1000), Y = runif(1000),
    a1 = rnorm(1000)
  )
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  s1 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = binomial(),
    range = 0.5,
    phi = 0.001,
    # sigma_O = 0.1 puts the binomial breakpt likelihood on a poorly-scaled
    # ridge near the threshold cutpoint (max gradient ~0.04-0.26 even with
    # newton polishing, varying by seed); 0.05 is well within the basin and
    # converges to ~2e-5, comfortably below the 0.001 warning threshold
    sigma_O = 0.05,
    seed = 4,
    B = 0,
    threshold_coefs = c(0.5, 0.3)
  )
  s2 <- sdmTMB_simulate(
    formula = ~ 1 + breakpt(a1),
    data = predictor_dat,
    mesh = mesh,
    family = Gamma(link = "log"),
    range = 0.5,
    phi = 1000,
    sigma_O = 0.1,
    seed = 4,
    B = 0,
    threshold_coefs = c(0.3, 0.3)
  )

  plot(predictor_dat$a1, s1$observed)
  plot(predictor_dat$a1, s2$observed)

  s <- s1
  s$observed <- s1$observed * s2$observed
  s$a1 <- predictor_dat$a1
  s1$a1 <- predictor_dat$a1
  s2$a1 <- predictor_dat$a1

  ctrl <- sdmTMBcontrol(nlminb_loops = 3L, newton_loops = 3L)

  # binomial works:
  fit1 <- sdmTMB(observed ~ breakpt(a1),
    data = s1,
    family = binomial(),
    # mesh = mesh,
    spatial = "off",
    control = ctrl
  )
  print(fit1)

  s2_pos <- subset(s2, s1$observed > 0)
  # mesh2 <- make_mesh(s2_pos, xy_cols = c("X", "Y"), mesh = mesh$mesh)
  # Gamma works:
  fit2 <- sdmTMB(observed ~ breakpt(a1),
    data = s2_pos,
    family = Gamma(link = "log"),
    spatial = "off",
    # mesh = mesh2,
    control = ctrl
  )
  print(fit2)

  fit <- sdmTMB(
    observed ~ breakpt(a1),
    data = s,
    family = delta_gamma(),
    # mesh = mesh,
    spatial = "off",
    control = ctrl
  )
  print(fit)
  sanity(fit)

  t1 <- tidy(fit1)
  t2 <- tidy(fit2)
  td1 <- tidy(fit, model = 1)
  td2 <- tidy(fit, model = 2)

  # standalone and joint-delta fits land at very slightly different optima
  expect_equal(t1$estimate, td1$estimate, tolerance = 1e-4)
  expect_equal(t2$estimate, td2$estimate, tolerance = 1e-4)
  expect_equal(t1$std.error, td1$std.error, tolerance = 1e-3)
  expect_equal(t2$std.error, td2$std.error, tolerance = 1e-3)
})
