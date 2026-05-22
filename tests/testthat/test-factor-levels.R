test_that("Test that droplevels matches lm()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1)
  df <- data.frame(
    y = rnorm(100),
    a_char = sample(c("a", "b", "c", "d", "e"), size = 100, replace = T)
  )
  df$a_fac <- as.factor(df$a_char)
  df$a_extra_fac <- factor(df$a_fac, levels = c("a", "b", "c", "d", "e", "f"))

  fit_lm <- lm(y ~ -1 + a_extra_fac, data = df)
  fit_sdmTMB <- sdmTMB(y ~ -1 + a_extra_fac, data = df, spatial = FALSE)
  expect_equal(as.numeric(coef(fit_lm)), tidy(fit_sdmTMB)$estimate)

  # prediction to new levels fails on both
  newdf <- data.frame(a_char = sample(c("a", "b", "c", "d", "e", "f"), size = 100, replace = TRUE))
  newdf$a_fac <- as.factor(newdf$a_char)
  fit_lm <- lm(y ~ -1 + a_fac, data = df)
  fit_sdmTMB <- sdmTMB(y ~ -1 + a_fac, data = df, spatial = FALSE)
  expect_error(predict(fit_lm, newdf), regexp = "new levels")
  expect_error(predict(fit_sdmTMB, newdf), regexp = "new levels")

  # prediction with missing factor levels behaves the same
  fit_lm <- lm(y ~ -1 + a_fac, data = df)
  fit_sdmTMB <- sdmTMB(y ~ -1 + a_fac, data = df, spatial = FALSE)
  newdf <- df
  newdf <- newdf[newdf$a_fac != "a", , drop = FALSE]
  p_lm <- as.numeric(predict(fit_lm, newdata = newdf))
  p_sdmTMB <- predict(fit_sdmTMB, newdata = newdf)$est
  expect_equal(p_lm, p_sdmTMB)
})

test_that("Test that droplevels matches glmmTMB on (1 | factor)", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  d <- pcod
  d$fyear <- as.factor(d$year)
  fit_glmmTMB <- glmmTMB::glmmTMB(density ~ 1 + (1 | fyear), data = d, family = glmmTMB::tweedie())
  fit_sdmTMB <- sdmTMB(density ~ 1 + (1 | fyear), data = d, family = glmmTMB::tweedie(), spatial = FALSE)

  r1 <- ranef(fit_glmmTMB)$cond$fyear[, 1]
  r2 <- tidy(fit_sdmTMB, "ran_vals")$estimate
  expect_equal(r1, r2, tolerance = 1e-3)

  # extra level not included:
  d$fyear <- factor(d$fyear, levels = c(
    "2003", "2004", "2005", "2007", "2009", "2011", "2013", "2015",
    "2017", "9999"
  ))

  fit_glmmTMB <- glmmTMB::glmmTMB(density ~ 1 + (1 | fyear), data = d, family = glmmTMB::tweedie())
  fit_sdmTMB <- sdmTMB(density ~ 1 + (1 | fyear), data = d, family = glmmTMB::tweedie(), spatial = FALSE)
  r1 <- ranef(fit_glmmTMB)$cond$fyear[, 1]
  r2 <- tidy(fit_sdmTMB, "ran_vals")$estimate
  expect_equal(r1, r2, tolerance = 1e-3)

  # new level on predict
  nd <- d
  nd$fyear <- factor(nd$fyear, levels = c(
    "2003", "2004", "2005", "2007", "2009", "2011", "2013", "2015",
    "2017", "9999", "9998"
  ))
  p1 <- predict(fit_glmmTMB, newdata = nd, re.form = NULL)
  p2 <- predict(fit_sdmTMB, newdata = nd)$est
  expect_equal(p1, p2, tolerance = 1e-3)

  # drop level on predict
  nd <- d
  nd <- nd[nd$fyear != "2003", ]
  nd$fyear <- factor(nd$fyear, levels = c(
    "2003", "2004", "2005", "2007", "2009", "2011", "2013", "2015",
    "2017", "9999", "9998"
  ))

  p1 <- predict(fit_glmmTMB, newdata = nd, re.form = NULL)
  p2 <- predict(fit_sdmTMB, newdata = nd)$est
  expect_equal(p1, p2, tolerance = 1e-3)
})

test_that("re_form_iid is not specified but new levels in newdata doesn't blow up", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  sub <- pcod[pcod$year != 2017, ]
  sub$fyear <- as.factor(sub$year)
  fit <- sdmTMB(density ~ 1 + (1 | fyear),
    data = sub,
    family = tweedie(link = "log"),
    spatial = "off"
  )
  d <- pcod
  d$fyear <- as.factor(d$year)
  p <- predict(fit, newdata = d, re_form_iid = NA)
  expect_warning(
    p1 <- predict(fit, newdata = d),
    regexp = "Found new levels"
  )
  p2 <- predict(fit, newdata = d, allow_new_levels = TRUE)
  expect_equal(p1$est, p2$est)
  expect_equal(p1$est[d$fyear == "2017"], p$est[d$fyear == "2017"])

  # what about just with 1 level?
  fit_glmmTMB <- glmmTMB::glmmTMB(density ~ 1 + (1 | fyear),
    data = sub,
    family = glmmTMB::tweedie(link = "log")
  )
  nd <- sub[sub$year == 2009, ]
  p_glmmTMB <- predict(fit_glmmTMB, newdata = nd)
  p <- predict(fit, newdata = nd)$est
  expect_equal(p_glmmTMB, p, tolerance = 1e-4)
})

test_that("predict handles dropped factor levels in time_varying factors", {
  set.seed(1)

  dat <- data.frame(
    year = rep(2010:2014, each = 30),
    X = runif(150),
    Y = runif(150),
    substrate = factor(
      sample(c("mud", "sand", "gravel"), 150, replace = TRUE),
      levels = c("mud", "sand", "gravel")
    )
  )
  year_effect <- setNames(seq(-0.6, 0.6, length.out = 5), 2010:2014)
  substrate_effect <- c(mud = -0.4, sand = 0.3, gravel = 0.8)
  dat$y <- 1 +
    0.5 * dat$X -
    0.4 * dat$Y +
    year_effect[as.character(dat$year)] * substrate_effect[as.character(dat$substrate)] +
    rnorm(nrow(dat), sd = 0.1)

  for (tv_formula in list(~ substrate, ~ 0 + substrate)) {
    fit <- sdmTMB(
      y ~ X + Y,
      time_varying = tv_formula,
      data = dat,
      time = "year",
      spatial = "off",
      spatiotemporal = "off",
      family = gaussian(),
      control = sdmTMBcontrol(newton_loops = 0L, getsd = FALSE),
      priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5))
    )

    nd <- dat[dat$substrate == "sand", ]
    nd <- nd[1:10, ]
    nd$substrate <- droplevels(nd$substrate)

    expect_equal(levels(nd$substrate), "sand")
    expect_no_error(predict(fit, newdata = nd))

    nd_char <- nd
    nd_char$substrate <- as.character(nd_char$substrate)
    expect_no_error(predict(fit, newdata = nd_char))

    nd_new_level <- nd
    nd_new_level$substrate <- factor("shell", levels = "shell")
    expect_error(predict(fit, newdata = nd_new_level), regexp = "new level")
  }
})
