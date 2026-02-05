test_that("SVC are estimated correctly for binomial and delta models", {
  skip_on_cran()
  local_edition(2)
  d <- pcod
  d$year_scaled <- as.numeric(scale(d$year))
  mesh10 <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m1 <- sdmTMB(
    data = d,
    formula = present ~ 1 + year_scaled,
    spatial_varying = ~ 0 + year_scaled,
    mesh = mesh10,
    family = binomial()
  )

  p <- predict(m1)
  pnd <- predict(m1, newdata = d)
  expect_identical(names(p), names(pnd))
  expect_equal(p$est, pnd$est)
  expect_equal(p$zeta_s_year_scaled, pnd$zeta_s_year_scaled)

  # m1.1 <- sdmTMB(
  #   data = d,
  #   formula = present ~ 1 + year_scaled,
  #   spatial_varying = ~ 1 + year_scaled, #<
  #   spatial = "off", #<
  #   mesh = mesh10,
  #   family = binomial()
  # )
  # expect_equal(m1$model$objective, m1.1$model$objective)

  b1 <- tidy(m1, effects = "ran_pars", conf.int = TRUE)
  expect_equal(b1$estimate[3], 0.312, tolerance = 0.1)

  m1.2 <- sdmTMB(
    data = d,
    formula = present ~ 1 + year_scaled,
    spatial_varying = ~ 1 + year_scaled,
    spatial = "on",
    mesh = mesh10,
    family = binomial()
  )

  expect_equal(m1$model$objective, m1.2$model$objective, tolerance = 1e-5)

  # warn: probably don't want to do this!
  expect_message({
    m1.3 <- sdmTMB(
      data = d,
      formula = present ~ 1 + year_scaled,
      spatial_varying = ~ 0 + as.factor(year), #<
      spatial = "on", #<
      mesh = mesh10,
      family = binomial(),
      do_fit = FALSE
    )
  }, regexp = "intercept")

  # better:
  m1.4 <- sdmTMB(
    data = d,
    formula = present ~ 1 + year_scaled,
    spatial_varying = ~ 0 + as.factor(year),
    spatial = "off",
    mesh = mesh10,
    family = binomial()
  )
  m1.4

  m1.5 <- sdmTMB(
    data = d,
    formula = present ~ 1 + year_scaled,
    spatial_varying = ~ 1 + as.factor(year),
    spatial = "on",
    mesh = mesh10,
    family = binomial()
  )
  m1.5
  p <- predict(m1.5, newdata = d)

  # also check that binomial portion of delta model matches the above
  m2 <- sdmTMB(
    data = d,
    formula = density ~ 1 + year_scaled,
    spatial_varying = ~ 0 + year_scaled,
    mesh = mesh10,
    family = delta_gamma()
  )
  b2 <- tidy(m2, effects = "ran_pars", conf.int = TRUE)
  expect_equal(b2$estimate[3], b1$estimate[3], tolerance = 1e-3)
})

test_that("Delta model with spatially varying factor predictor and no spatiotemporal field works #237", {
  # https://github.com/sdmTMB/sdmTMB/issues/237
  skip_on_cran()
  pcod_q2 <- pcod_2011
  pcod_q1 <- pcod_2011
  pcod_q1$quarter <- as.factor(1)
  pcod_q2$quarter <- as.factor(2)
  set.seed(123)
  pcod_q2$density <- pcod_q2$density + rnorm(10, 20, n = nrow(pcod_2011)) # just adding some difference between quarters..
  pcod2 <- rbind(pcod_q1, pcod_q2)
  # Fit delta model with spatially varying quarter effect
  mesh <- make_mesh(pcod2, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + as.factor(year) + quarter,
    data = pcod2,
    mesh = mesh,
    family = delta_gamma(link1 = "logit", link2 = "log"),
    spatiotemporal = "off",
    spatial = "off", # since spatially varying predictor is a factor
    spatial_varying = ~0 + quarter,
    time = "year",
    control = sdmTMBcontrol(newton_loops = 1L)
  )
  expect_s3_class(m, "sdmTMB")
  expect_true(sum(is.na(m$sd_report$sd)) == 0L)
})

test_that("Factor handling for SVC models works #269", {
  skip_on_cran()
  set.seed(1)
  pcod_2011$vessel <- sample(c("A", "B"), size = nrow(pcod_2011), replace = TRUE)
  pcod_2011$vessel <- as.factor(pcod_2011$vessel)
  fit <- sdmTMB(present ~ vessel,
    spatial_varying = ~ vessel,
    spatial = "on",
    mesh = pcod_mesh_2011,
    data = pcod_2011
  )
  p1 <- predict(fit, pcod_2011)
  p2 <- predict(fit, newdata = pcod_2011)
  expect_equal(p1$est, p2$est)

  p3 <- predict(fit, newdata = pcod_2011[pcod_2011$vessel == "A", ])
  p4 <- p2[p2$vessel == "A", ]
  expect_equal(p3$est, p4$est)
})

test_that("SVC throws a warning if character class #269", {
  skip_on_cran()
  pcod_2011$vessel <- sample(c("A", "B"), size = nrow(pcod_2011), replace = TRUE)
  expect_warning({
    fit <- sdmTMB(present ~ vessel,
      spatial_varying = ~ vessel,
      spatial = "on",
      mesh = pcod_mesh_2011,
      data = pcod_2011
    )
  }, regexp = "character")
})

