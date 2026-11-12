
# https://github.com/sdmTMB/sdmTMB/issues/60
test_that("smoothers with 'bs = re' error", {
  skip_on_cran()
  expect_error({
    m <- sdmTMB(
      density ~ s(depth_scaled, bs = "re"),
      data = pcod_2011,
      mesh = pcod_mesh_2011, spatial = "off"
    )
  }, regexp = "re")
})

test_that("A model with 2 s() splines works", {
  skip_on_cran()
  d <- subset(pcod, year >= 2000 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
    mesh = pcod_spde, spatial = "off", spatiotemporal = "off"
  )
  expect_equal(ncol(m$tmb_data$X_ij[[1]]), 1L)
  expect_equal(length(m$tmb_data$Zs), 2L)
  # head(m$tmb_data$Zs[[1]])
  # head(m$tmb_data$Zs[[2]])
  expect_equal(ncol(m$tmb_data$Xs), 2L)
  expect_equal(m$tmb_data$b_smooth_start, c(0L, 8L))
  # exp(as.list(m$sd_report, "Estimate")$ln_smooth_sigma)
  # as.list(m$sd_report, "Estimate")$b_smooth
  expect_equal(sum(is.na(as.list(m$sd_report, "Std. Error")$b_smooth)), 0L)
  p <- predict(m, newdata = NULL)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled) + s(year, k = 5), data = d)
  p_mgcv <- predict(m_mgcv, newdata = d)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.999)

  pnd <- predict(m, newdata = d[1:8,])
  expect_equal(p$est[1:8], pnd$est, tolerance = 0.001)

  set.seed(23402)
  .s <- sample(seq_len(nrow(d)), 200L)
  pnd_mgcv <- predict(m_mgcv, newdata = d[.s, ])
  pnd <- predict(m, newdata = d[.s, ])
  expect_gt(cor(pnd_mgcv, pnd$est), 0.999)
})

test_that("A model with t2() works", {
  skip_on_cran()
  set.seed(2938)
  dat <- mgcv::gamSim(1, n = 400, dist = "normal", scale = 2)

  dat$f <- NULL
  dat$f0 <- NULL
  dat$f1 <- NULL
  dat$f2 <- NULL
  dat$f3 <- NULL
  dat$x3 <- NULL
  dat$observed <- dat$y
  dat$y <- NULL
  dat$.x0 <- dat$x0
  dat$.x1 <- dat$x1
  dat$x0 <- NULL
  dat$x1 <- NULL
  dat$x2 <- NULL

  m_mgcv <- mgcv::gam(observed ~ t2(.x0, .x1, k = 4),
                      data = dat,
                      method = "REML"
  )
  p_mgcv <- predict(m_mgcv)
  m <- sdmTMB(observed ~ t2(.x0, .x1, k = 4),
              data = dat,
              spatial = 'off', reml = TRUE
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
  expect_equal(as.numeric(p$est), as.numeric(p_mgcv), tolerance = 0.001)

  # prediction on newdata
  pnd <- predict(m, newdata = dat)
  # plot(p$est, pnd$est)
  expect_equal(as.numeric(p$est), as.numeric(pnd$est), tolerance = 0.00001)
})

test_that("A model with dimensions specified in t2() works", {
  skip_on_cran()
  set.seed(2938)
  dat <- mgcv::gamSim(1, n = 400, dist = "normal", scale = 1)

  dat$f <- NULL
  dat$f0 <- NULL
  dat$f1 <- NULL
  dat$f2 <- NULL
  dat$f3 <- NULL
  #dat$x3 <- NULL
  dat$observed <- dat$y
  dat$y <- NULL
  dat$.x0 <- dat$x0
  dat$.x1 <- dat$x1
  dat$.x2 <- dat$x2
  #dat$x0 <- NULL
  #dat$x1 <- NULL
  #dat$x2 <- NULL

  m_mgcv <- mgcv::gam(observed ~ t2(x0, x1, x2, d = c(2,1), k = c(5,3)),
                      data = dat,
                      method = "REML"
  )
  p_mgcv <- predict(m_mgcv)
  m <- sdmTMB(observed ~ t2(x0, x1, x2, d = c(2,1), k = c(5,3)),
              data = dat,
              spatial = 'off'
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
  expect_equal(as.numeric(p$est), as.numeric(p_mgcv), tolerance = 0.001)
})


test_that("A model with by in spline (and s(x, y)) works", {
  skip_on_cran()
  set.seed(19203)
  # examples from ?mgcv::gam.models
  # continuous by example:
  dat <- mgcv::gamSim(3, n = 800)
  m_mgcv <- mgcv::gam(y ~ s(x2, by = x1), data = dat)
  p_mgcv <- predict(m_mgcv)
  dat$X <- runif(nrow(dat))
  dat$Y <- runif(nrow(dat))
  spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
  m <- sdmTMB(y ~ s(x2, by = x1),
              data = dat,
              mesh = spde,spatial = "off"
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
  pnd <- predict(m, newdata = dat)
  plot(p$est, pnd$est)
  expect_gt(cor(p$est, pnd$est), 0.9999)

  # s(x, y)
  m_mgcv <- mgcv::gam(y ~ s(x2, x1), data = dat)
  p_mgcv <- predict(m_mgcv)
  m <- sdmTMB(y ~ s(x2, x1),
              data = dat,
              mesh = spde, spatial = "off"
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.999)
  p2 <- predict(m, newdata = dat)
  plot(p2$est, p$est)
  expect_gt(cor(p2$est, p$est), 0.999)

  pnd <- predict(m, newdata = dat)
  plot(p$est, pnd$est)
  expect_gt(cor(p$est, pnd$est), 0.9999)


  # Factor `by' variable example (with a spurious covariate x0)
  set.seed(1)
  dat <- mgcv::gamSim(4)
  m_mgcv <- mgcv::gam(y ~ fac + s(x2, by = fac) + s(x0), data = dat)
  p_mgcv <- predict(m_mgcv)
  dat$X <- runif(nrow(dat))
  dat$Y <- runif(nrow(dat))
  spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
  m <- sdmTMB(y ~ fac + s(x2, by = fac) + s(x0),
              data = dat,
              mesh = spde, spatial = "off"
  )
  p <- predict(m, newdata = NULL)
  plot(p$est, p_mgcv)
  abline(a = 0, b = 1)
  expect_gt(cor(p$est, p_mgcv), 0.9999)
  pnd <- predict(m, newdata = dat)
  plot(p$est, pnd$est)
  expect_gt(cor(p$est, pnd$est), 0.9999)

  set.seed(291823)
  .s <- sample(seq_len(nrow(dat)), 200L)
  pnd <- predict(m, newdata = dat[.s,])
  expect_equal(p$est[.s], pnd$est, tolerance = 0.001)
})

test_that("Formula removal of s and t2 works", {
  expect_identical(remove_s_and_t2(y ~ x + s(z)), y ~ x)
  expect_identical(remove_s_and_t2(y ~ s(x) + s(z)), y ~ 1)
  expect_identical(remove_s_and_t2(y ~ s(x) + t2(z)), y ~ 1)
  expect_identical(remove_s_and_t2(y ~ ignore(x) + t2(z)), y ~ ignore(x))
})

test_that("Smooth plotting works", {
  skip_on_cran()
  d <- subset(pcod, year >= 2000 & density > 0)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
    mesh = pcod_spde, spatial = "off", spatiotemporal = "off"
  )

  plot_smooth(m)
  plot_smooth(m, level = 0.4)
  plot_smooth(m, select = 2)
  plot_smooth(m, ggplot = TRUE)
  plot_smooth(m, ggplot = TRUE, rug = FALSE)
  plot_smooth(m, rug = FALSE)
  plot_smooth(m, n = 3)
  out <- plot_smooth(m, return_data = TRUE)
  expect_equal(class(out), "data.frame")

  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled) + year,
    mesh = pcod_spde, spatial = "off", spatiotemporal = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )

  plot_smooth(m, select = 1)
  expect_error(plot_smooth(m, select = 2), regexp = "select")

  suppressMessages({
    m <- sdmTMB(
      data = d, time = "year",
      formula = log(density) ~ s(depth_scaled),
      mesh = pcod_spde, spatial = "on", spatiotemporal = "off"
    )
  })
  plot_smooth(m)

  # with a factor
  suppressMessages({
    m1 <- sdmTMB(
      data = d, time = "year",
      formula = log(density) ~ 0 + s(depth_scaled) + as.factor(year),
      mesh = pcod_spde, spatial = "on", spatiotemporal = "off"
    )
  })
  plot_smooth(m)
})

test_that("print works with s(X, Y)", {
  skip_on_cran()

  d <- subset(pcod_2011, density > 0)
  m <- sdmTMB(log(density) ~ s(X, Y, k = 5), data = d, spatial = "off")
  print(m)
  expect_output(summary(m), regexp = "sXY_1")
  expect_output(summary(m), regexp = "sXY_2")

  m2 <- sdmTMB(log(density) ~ s(X) + s(Y), data = d, spatial = "off")
  summary(m2)

  m3 <- sdmTMB(log(density) ~ s(X) + s(Y) + s(depth, year), data = d, spatial = "off")
  summary(m3)
  expect_output(summary(m3), regexp = "sX")
  expect_output(summary(m3), regexp = "sY")
  expect_output(summary(m3), regexp = "sdepthyear_1")
  expect_output(summary(m3), regexp = "sdepthyear_2")
  # m <- brms::brm(log(density) ~ s(X) + s(Y) + s(depth, year), data = d, chains = 1, iter = 800)
  # m
})

test_that("smoothers with 'bs = cc' work", {
  skip_on_cran()
  m <- sdmTMB(
    density ~ s(depth_scaled, bs = "cc"),
    data = pcod_2011,
    mesh = pcod_mesh_2011, spatial = "off"
  )
  p <- predict(m)
  print(m)

  m1 <- mgcv::gam(
    density ~ s(depth_scaled, bs = "cc"),
    data = pcod_2011
  )
  p1 <- predict(m1)

  plot(p$est, p1)
  abline(0, 1)
  expect_gt(cor(p$est, p1), 0.999)
})


test_that("smoothers with 'bs = cc' work with knots specified", {
  skip_on_cran()
  m <- sdmTMB(
    density ~ s(depth_scaled, bs = "cc", k = 5), knots = list(depth_scaled = c(-3, -1, 0, 1, 3)),
    data = pcod_2011, spatial = "off"
  )
  p <- predict(m, newdata = pcod_2011)
  print(m)

  m1 <- mgcv::gam(
    density ~ s(depth_scaled, bs = "cc", k = 5), knots = list(depth_scaled = c(-3, -1, 0, 1, 3)),
    data = pcod_2011
  )
  p1 <- predict(m1)

  plot(p$est, p1)
  abline(0, 1)
  expect_gt(cor(p$est, p1), 0.999)
})

test_that("smoothers with 'bs = cr' work", {
  skip_on_cran()
  m <- sdmTMB(
    density ~ s(depth_scaled, bs = "cr"),
    data = pcod_2011,
    mesh = pcod_mesh_2011, spatial = "off"
  )
  p <- predict(m)
  print(m)

  m1 <- mgcv::gam(
    density ~ s(depth_scaled, bs = "cr"),
    data = pcod_2011
  )
  p1 <- predict(m1)

  plot(p$est, p1)
  abline(0, 1)
  expect_gt(cor(p$est, p1), 0.999)
})

test_that("prediction with smoothers error helpfully if missing variable", {
  skip_on_cran()
  suppressWarnings({
    m <- sdmTMB(
      density ~ s(year, k = 3) + s(depth_scaled),
      data = pcod_2011,
      mesh = pcod_mesh_2011, spatial = "off"
    )
  })
  nd <- data.frame(year = 2007)
  expect_error({
    p <- predict(m, newdata = nd, re_form = NA, se_fit = TRUE)
  }, regexp = "depth_scaled")

  nd <- data.frame(aaa = 2007)
  expect_error({
    p <- predict(m, newdata = nd, re_form = NA, se_fit = TRUE)
  }, regexp = "year")
})

test_that("A model with s(x, bs = 'cs') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "cs"),
    spatial = "off"
  )
  print(m)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "cs"), data = d)
  p <- predict(m)
  p2 <- predict(m_mgcv)
  # plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.999)
})

test_that("A model with s(x, bs = 'cr') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "cr"),
    spatial = "off"
  )
  print(m)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "cr"), data = d)
  p <- predict(m)
  p2 <- predict(m_mgcv)
  # plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.999)
})

test_that("A model with s(x, bs = 'ds') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "ds"),
    spatial = "off"
  )
  print(m)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "ds"), data = d)
  p <- predict(m)
  p2 <- predict(m_mgcv)
  # plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.999)
})

test_that("A model with s(x, bs = 'ps') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "ps"),
    spatial = "off"
  )
  print(m)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "ps"), data = d, method = "REML")
  p <- predict(m)
  p2 <- predict(m_mgcv)
  plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.999)
})

test_that("A model with s(x, bs = 're') errors", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  expect_error(m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "re"),
    spatial = "off"
  ))
})

test_that("A model with s(x, bs = 'gp') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, bs = "gp"),
    spatial = "off"
  )
  print(m)
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, bs = "gp"), data = d, method = "REML")
  p <- predict(m)
  p2 <- predict(m_mgcv)
  plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.999)
})

test_that("A model with s(x, bs = 'fs') works", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  d$yearf <- as.factor(d$year)
  m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, by = year, bs = "fs"),
    spatial = "off", reml = TRUE
  )
  suppressWarnings(print(m))
  m_mgcv <- mgcv::gam(log(density) ~ s(depth_scaled, by = year, bs = "fs"), data = d, method = "REML")
  p <- predict(m)
  p2 <- predict(m_mgcv)
  plot(p$est, p2)
  expect_gt(stats::cor(p$est, p2), 0.9999)
})

test_that("An fx=TRUE smoother errors out", {
  skip_on_cran()
  d <- subset(pcod, density > 0)
  expect_error(m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, fx = TRUE),
    spatial = "off"
  ))
  expect_error(m <- sdmTMB(
    data = d,
    formula = log(density) ~ s(depth_scaled, fx = T),
    spatial = "off"
  ))
})

test_that("print works with t2()", {
  skip_on_cran()
  skip_on_ci()

  # Test basic t2() printing
  fit <- sdmTMB(density ~ t2(X, Y, k = 4),
    spatial = "off",
    data = pcod
  )

  # Check that coefficient names are indexed correctly using tidy
  coef_tidy <- tidy(fit, silent = TRUE)
  expect_true(any(grepl("^t2\\(X,Y\\)_1$", coef_tidy$term)))
  expect_true(any(grepl("^t2\\(X,Y\\)_2$", coef_tidy$term)))
  expect_true(any(grepl("^t2\\(X,Y\\)_3$", coef_tidy$term)))
  # Check smooth SD names via tidy ran_pars
  ran_pars <- tidy(fit, effects = "ran_pars", silent = TRUE)
  expect_true(any(grepl("^sd__t2\\(X,Y\\)_1$", ran_pars$term)))
  expect_true(any(grepl("^sd__t2\\(X,Y\\)_2$", ran_pars$term)))
  expect_true(any(grepl("^sd__t2\\(X,Y\\)_3$", ran_pars$term)))

  # add in an s() too:
  fit2 <- sdmTMB(density ~ t2(X, Y, k = 4) + s(year, k = 3),
    spatial = "off",
    data = pcod
  )
  print(fit2)
  ran_pars <- tidy(fit2, effects = "ran_pars", silent = TRUE)
  expect_true(any(grepl("^sd__s\\(year\\)$", ran_pars$term)))
  expect_true(any(grepl("^sd__t2\\(X,Y\\)_1$", ran_pars$term)))
})

test_that("smoother linear terms are named with by smooths and cyclic smooths", {
  d <- pcod_2011
  d$fyear <- factor(d$year)
  fit <- sdmTMB(
    density ~ s(depth_scaled, by = fyear, k = 4) + s(X, bs = "cc", k = 5) +
      s(Y, k = 4),
    data = d, spatial = "off", family = tweedie()
  )
  # one linear column per `by` level, none for the cyclic smooth:
  expect_identical(
    tidy(fit)$term[-1],
    c(paste0("sdepth_scaled):fyear", c(2011, 2013, 2015, 2017)), "sY")
  )
})
