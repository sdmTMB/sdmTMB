test_that("Ohio areal SAR/CAR estimates and logLik stay stable", {
  skip_on_cran()

  dat <- ohio_df
  domain <- make_areal_domain(ohio_sf, id_column = "county")

  fit_sar <- sdmTMB(
    cases ~ 0 + as.factor(year) + pct_male,
    data = dat,
    mesh = domain,
    spatial_model = "sar",
    time = "year",
    family = stats::poisson(link = "log"),
    spatial = "on",
    spatiotemporal = "iid",
    offset = log(dat$pop)
  )

  fit_car <- sdmTMB(
    cases ~ 0 + as.factor(year) + pct_male,
    data = dat,
    mesh = domain,
    spatial_model = "car",
    time = "year",
    family = stats::poisson(link = "log"),
    spatial = "on",
    spatiotemporal = "iid",
    offset = log(dat$pop)
  )

  expect_equal(as.numeric(logLik(fit_sar)), -825.4355, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit_car)), -825.1994, tolerance = 1e-4)

  beta_sar <- coef(fit_sar)
  beta_car <- coef(fit_car)

  expect_equal(beta_sar[["as.factor(year)1968"]], -7.727077, tolerance = 1e-4)
  expect_equal(beta_sar[["as.factor(year)1978"]], -7.405386, tolerance = 1e-4)
  expect_equal(beta_sar[["as.factor(year)1988"]], -7.195807, tolerance = 1e-4)
  expect_equal(beta_sar[["pct_male"]], 0.076565, tolerance = 1e-4)

  expect_equal(beta_car[["as.factor(year)1968"]], -7.722861, tolerance = 1e-4)
  expect_equal(beta_car[["as.factor(year)1978"]], -7.400958, tolerance = 1e-4)
  expect_equal(beta_car[["as.factor(year)1988"]], -7.191256, tolerance = 1e-4)
  expect_equal(beta_car[["pct_male"]], 0.055420, tolerance = 1e-4)

  ran_sar <- tidy(fit_sar, "ran_pars")
  ran_car <- tidy(fit_car, "ran_pars")

  expect_equal(ran_sar$estimate[ran_sar$term == "rho_sar"], 0.3128827, tolerance = 1e-3)
  expect_equal(ran_car$estimate[ran_car$term == "alpha_car"], 0.5966396, tolerance = 1e-3)
})
