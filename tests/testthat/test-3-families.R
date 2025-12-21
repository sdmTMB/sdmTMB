test_that("Families return a name to list with the correct names", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_true(all(.names %in% names(student(link = "identity"))))
  expect_true(all(.names %in% names(lognormal(link = "log"))))
  expect_true(all(.names %in% names(tweedie(link = "log"))))
  expect_true(all(.names %in% names(nbinom2(link = "log"))))
})

test_that("The supplementary families work with appropriate links", {
  expect_identical(class(tweedie(link = "log")), "family")
  expect_identical(class(tweedie(link = log)), "family")
  expect_error(class(tweedie(link = "banana")))
  expect_error(class(tweedie(link = banana)))

  expect_identical(class(lognormal(link = "log")), "family")
  expect_identical(class(lognormal(link = log)), "family")
  expect_error(class(lognormal(link = "banana")))
  expect_error(class(lognormal(link = banana)))

  expect_identical(class(nbinom2(link = "log")), "family")
  expect_identical(class(nbinom2(link = log)), "family")
  expect_error(class(nbinom2(link = "banana")))
  expect_error(class(nbinom2(link = inverse)))

  expect_identical(class(student(link = "identity")), "family")
  expect_identical(class(student(link = identity)), "family")
  expect_error(class(student(link = "banana")))
  expect_error(class(student(link = banana)))
})

set.seed(1)
x <- stats::runif(100, -1, 1)
y <- stats::runif(100, -1, 1)
loc <- data.frame(x = x, y = y)

spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

test_that("Student family fits", {
  skip_on_cran()
  set.seed(3)
  initial_betas <- 0.5
  range <- 0.5
  sigma_O <- 0.3
  phi <- 0.01
  s <- sdmTMB_simulate(~ 1, data = loc,
    B = initial_betas,
    phi = phi, range = range, sigma_O = sigma_O, sigma_E = 0,
    seed = 1, mesh = spde
  )
  m <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde,
    family = student(link = "identity", df = 7),
    spatial = "off", spatiotemporal = "off"
  )
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(s))
})

test_that("Student-t df parameter (fixed vs estimated)", {
  skip_on_cran()
  set.seed(1234)
  # Simulate data from Student-t with known df
  x <- stats::runif(300, -1, 1)
  y <- stats::runif(300, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  # Simulate response with Student-t errors
  true_df <- 5
  mu <- 2
  sigma <- 0.3
  loc$observed <- mu + sigma * rt(nrow(loc), df = true_df)

  # Test 1: Estimated df (new default)
  suppressMessages(
    m_est_default <- sdmTMB(data = loc, formula = observed ~ 1, mesh = spde,
      family = student(link = "identity"),  # df = NULL by default, will be estimated
      spatial = "off", spatiotemporal = "off"
    )
  )
  # Check that ln_student_df is not mapped (estimated)
  expect_null(m_est_default$tmb_map[["ln_student_df"]])
  # Check tidy output shows estimated df with std.error
  tidy_est_default <- tidy(m_est_default, "ran_pars")
  df_row_est_default <- tidy_est_default[tidy_est_default$term == "student_df", ]
  expect_true(df_row_est_default$estimate > 1)
  expect_false(is.na(df_row_est_default$std.error))

  # Test 2: Fixed df at 3
  suppressMessages(
    m_fixed_3 <- sdmTMB(data = loc, formula = observed ~ 1, mesh = spde,
      family = student(link = "identity", df = 3),
      spatial = "off", spatiotemporal = "off"
    )
  )
  expect_true(is.na(m_fixed_3$tmb_map$ln_student_df[1]))
  tidy_3 <- tidy(m_fixed_3, "ran_pars")
  df_row_3 <- tidy_3[tidy_3$term == "student_df", ]
  expect_equal(df_row_3$estimate, 3, tolerance = 1e-6)
  expect_true(is.na(df_row_3$std.error))

  # Test 3: Fixed df at custom value (5)
  suppressMessages(
    m_fixed_5 <- sdmTMB(data = loc, formula = observed ~ 1, mesh = spde,
      family = student(link = "identity", df = 5),
      spatial = "off", spatiotemporal = "off"
    )
  )
  expect_true(is.na(m_fixed_5$tmb_map$ln_student_df[1]))
  tidy_5 <- tidy(m_fixed_5, "ran_pars")
  df_row_5 <- tidy_5[tidy_5$term == "student_df", ]
  expect_equal(df_row_5$estimate, 5, tolerance = 1e-6)
  expect_true(is.na(df_row_5$std.error))

  # Test 4: Explicitly estimated df with df = NULL
  suppressMessages(
    m_est_explicit <- sdmTMB(data = loc, formula = observed ~ 1, mesh = spde,
      family = student(link = "identity", df = NULL),  # explicitly estimate df
      spatial = "off", spatiotemporal = "off"
    )
  )
  # Check that ln_student_df is not mapped (estimated)
  expect_null(m_est_explicit$tmb_map[["ln_student_df"]])
  # Check tidy output shows estimated df with std.error
  tidy_est_explicit <- tidy(m_est_explicit, "ran_pars")
  df_row_est_explicit <- tidy_est_explicit[tidy_est_explicit$term == "student_df", ]
  expect_true(df_row_est_explicit$estimate > 1)  # df should be > 1
  expect_false(is.na(df_row_est_explicit$std.error))  # should have std.error
  expect_false(is.na(df_row_est_explicit$conf.low))
  expect_false(is.na(df_row_est_explicit$conf.high))

  # Check print method includes df
  expect_output(print(m_est_explicit), "Student-t df")
})

test_that("Lognormal fits", {
  skip_on_cran()
  range <- 1
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 70, type = "kmeans")
  sigma_O <- 0.3
  sigma_E <- 0
  phi <- 0.2
  s <- sdmTMB_simulate(~ 1, loc, mesh = spde, family = lognormal(),
    B = 1,
    phi = phi, range = range, sigma_O = sigma_O, seed = 1
  )
  mlog <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde,
    family = lognormal(link = "log"))
  expect_equal(exp(mlog$model$par[["ln_phi"]]), phi, tolerance = 0.1)
})

test_that("NB2 fits", {
  skip_on_cran()
  set.seed(1)
  x <- stats::runif(300, -1, 1)
  y <- stats::runif(300, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 80, type = "kmeans")
  s <- sdmTMB_simulate(~ 1, loc, B = 0.4, phi = 1.5, range = 0.8,
    sigma_O = 0.4, seed = 1, mesh = spde, family = nbinom2())
  m <- sdmTMB(data = s, formula = observed ~ 1,
    mesh = spde, family = nbinom2(),
    control = sdmTMBcontrol(newton_loops = 1))
  expect_equal(round(tidy(m)[,"estimate", drop=TRUE], 6), 0.601897)
})

test_that("Truncated NB2, truncated NB1, and regular NB1 fit", {
  skip_on_cran()
  set.seed(1)
  x <- stats::runif(300, -1, 1)
  y <- stats::runif(300, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 80, type = "kmeans")
  s <- sdmTMB_simulate(~ 1, loc, B = 0.4, phi = 1.5, range = 0.8,
    sigma_O = 0.4, seed = 1, mesh = spde, family = nbinom2())

  m_sdmTMB <- sdmTMB(data = s, formula = observed ~ 1,
    mesh = spde, family = nbinom1(),
    spatial = "off")
  m_glmmTMB <- glmmTMB::glmmTMB(data = s, formula = observed ~ 1,
    family = glmmTMB::nbinom1())
  expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
  expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)

  s_trunc <- subset(s, observed > 0)
  spde <- make_mesh(s_trunc, c("x", "y"), n_knots = 80, type = "kmeans")
  m_sdmTMB <- sdmTMB(data = s_trunc, formula = observed ~ 1,
    mesh = spde, family = truncated_nbinom2(), spatial = "off")
  m_glmmTMB <- glmmTMB::glmmTMB(data = s_trunc, formula = observed ~ 1,
    family = glmmTMB::truncated_nbinom2())
  expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
  expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)

  m_sdmTMB <- sdmTMB(data = s_trunc, formula = observed ~ 1,
    mesh = spde, family = truncated_nbinom1(), spatial = "off")
  m_glmmTMB <- glmmTMB::glmmTMB(data = s_trunc, formula = observed ~ 1,
    family = glmmTMB::truncated_nbinom1())
  expect_equal(m_glmmTMB$fit$par[[1]], m_sdmTMB$model$par[[1]], tolerance = 0.00001)
  expect_equal(m_glmmTMB$fit$par[[2]], m_sdmTMB$model$par[[2]], tolerance = 0.00001)
})

test_that("Poisson fits", {
  skip_on_cran()
  d <- pcod
  spde <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  set.seed(3)
  d$density <- rpois(nrow(pcod), 3)
  m <- sdmTMB(data = d, formula = density ~ 1,
    mesh = spde, family = poisson(link = "log"),
    control = sdmTMBcontrol(newton_loops = 1)
  )
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(pcod))
})

test_that("Binomial fits", {
  skip_on_cran()
  d <- pcod[pcod$year == 2017, ]
  d$density <- round(d$density)
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  d$present <- ifelse(d$density > 0, 1, 0)
  m <- sdmTMB(data = d, formula = present ~ 1,
    mesh = spde, family = binomial(link = "logit"),
    control = sdmTMBcontrol(newton_loops = 1))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(d))
})

test_that("Gamma fits", {
  skip_on_cran()
  d <- pcod[pcod$year == 2017 & pcod$density > 0, ]
  spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(data = d, formula = density ~ 1,
    mesh = spde, family = Gamma(link = "log"),
    spatial = "off",
    control = sdmTMBcontrol(newton_loops = 1))
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(d))
  set.seed(123)
  d$test_gamma <- stats::rgamma(nrow(d), shape = 0.5, scale = 1 / 0.5)
  m <- sdmTMB(data = d, formula = test_gamma ~ 1,
    mesh = spde, family = Gamma(link = "inverse"), spatiotemporal = "off",
    control = sdmTMBcontrol(newton_loops = 1))
})

test_that("Beta fits", {
  skip_on_cran()
  set.seed(1)
  x <- stats::runif(400, -1, 1)
  y <- stats::runif(400, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 90, type = "kmeans")
  s <- sdmTMB_simulate(~ 1, loc, mesh = spde, sigma_O = 0.2,
    range = 0.8, family = Beta(), phi = 4, B = 1)
  m <- sdmTMB(data = s, formula = observed ~ 1,
    mesh = spde, family = Beta(link = "logit"),
    control = sdmTMBcontrol(newton_loops = 1),
    spatial = "off")
  expect_true(all(!is.na(summary(m$sd_report)[,"Std. Error"])))
  expect_length(residuals(m), nrow(s))

  m_glmmTMB<- glmmTMB::glmmTMB(data = s, formula = observed ~ 1,
    family = glmmTMB::beta_family(link = "logit"))

  expect_equal(m$model$par[[2]], m_glmmTMB$fit$par[[2]], tolerance = 1e-4)
  expect_equal(m$model$par[[1]], m_glmmTMB$fit$par[[1]], tolerance = 1e-4)
})

test_that("Censored Poisson fits", {
  skip_on_cran()
  set.seed(1)

  predictor_dat <- data.frame(X = runif(300), Y = runif(300))
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  sim_dat <- sdmTMB_simulate(
    formula = ~1,
    data = predictor_dat,
    mesh = mesh,
    family = poisson(),
    range = 0.5,
    sigma_O = 0.2,
    seed = 1,
    B = 2 # B0 = intercept
  )
  m_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = poisson(link = "log")
  )
  m_nocens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    mesh = mesh, family = censored_poisson(link = "log"),
    control = sdmTMBcontrol(censored_upper = sim_dat$observed)
  )
  expect_equal(m_nocens_pois$tmb_data$y_i[,1], m_nocens_pois$tmb_data$upr)
  expect_equal(names(m_nocens_pois$tmb_data$family), "censored_poisson")
  expect_equal(m_pois$model, m_nocens_pois$model)

  # # left-censored version
  # L_1 <- 5 # zeros and ones cannot be observed directly - observed as <= L1
  # y <- sim_dat$observed
  # lwr <- ifelse(y <= L_1, 0, y)
  # upr <- ifelse(y <= L_1, L_1, y)
  # m_left_cens_pois <- sdmTMB(
  #   data = sim_dat, formula = observed ~ 1,
  #   mesh = mesh, family = censored_poisson(link = "log"),
  #   control = sdmTMBcontrol(censored_lower = lwr, censored_upper = upr),
  #   spatial = "off"
  # )

  # right-censored version
  U_1 <- 8 # U_1 and above cannot be directly observed - instead we see >= U1
  y <- sim_dat$observed
  lwr <- ifelse(y >= U_1, U_1, y)
  upr <- ifelse(y >= U_1, NA, y)

  # old:
  expect_error(m_right_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    family = censored_poisson(link = "log"),
    experimental = list(lwr = lwr, upr = upr),
    spatial = "off"
  ), regexp = "upr")

  # new:
  m_right_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    family = censored_poisson(link = "log"),
    control = sdmTMBcontrol(censored_upper = upr),
    spatial = "off"
  )

  # interval-censored tough example
  # unique bounds per observation with upper limit 500 to test numerical underflow issues
  set.seed(123)
  U_2 <- sample(c(5:9), size = length(y), replace = TRUE)
  L_2 <- sample(c(1, 2, 3, 4), size = length(y), replace = TRUE)
  # lwr <- ifelse(y >= U_2, U_2, ifelse(y <= L_2, 0, y))
  upr <- ifelse(y >= U_2, 500, ifelse(y <= L_2, L_2, y))
  m_interval_cens_pois <- sdmTMB(
    data = sim_dat, formula = observed ~ 1,
    family = censored_poisson(link = "log"),
    control = sdmTMBcontrol(censored_upper = upr),
    spatial = "off"
  )
  expect_true(all(!is.na(summary(m_interval_cens_pois$sd_report)[, "Std. Error"])))

  # # reversed upr and lwr:
  # expect_error(
  #   m <- sdmTMB(
  #     data = sim_dat, formula = observed ~ 1,
  #     mesh = mesh, family = censored_poisson(link = "log"),
  #     experimental = list(lwr = upr, upr = lwr)
  #   ), regexp = "lwr")

  # wrong length lwr and upr
  expect_error(
    m <- sdmTMB(
      data = sim_dat, formula = observed ~ 1,
      mesh = mesh, family = censored_poisson(link = "log"),
      control = sdmTMBcontrol(censored_upper = c(4, 5, 6))
    ), regexp = "upr")

  # missing lwr/upr
  expect_error(
    m <- sdmTMB(
      data = sim_dat, formula = observed ~ 1,
      mesh = mesh, family = censored_poisson(link = "log"),
    ), regexp = "censored_upper")

})

test_that("Censored Poisson upper limit function works", {
  dat <- structure(
    list(
      n_catch = c(
        78L, 63L, 15L, 6L, 7L, 11L, 37L, 99L, 34L, 100L, 77L, 79L,
        98L, 30L, 49L, 33L, 6L, 28L, 99L, 33L
      ),
      prop_removed = c(
        0.61, 0.81, 0.96, 0.69, 0.99, 0.98, 0.25, 0.95, 0.89, 1, 0.95, 0.95,
        0.94, 1, 0.95, 1, 0.84, 0.3, 1, 0.99
      ), n_hooks = c(
        140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L,
        140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L
      )
    ),
    class = "data.frame", row.names = c(NA, -20L)
  )
  upr <- get_censored_upper(dat$prop_removed, dat$n_catch, dat$n_hooks, pstar = 0.9)
  expect_identical(upr,
    c(78, 63, 28, 6, 39, 34, 37, 109, 34, 140, 87, 89, 105, 70, 59, 73, 6, 28, 139, 65))
  x <- get_censored_upper(
    prop_removed = c(0.5, 0.3, 0.2),
    n_catch = c(3, 3, 3),
    n_hooks = c(5, 5, 5),
    pstar = 0.9
  )
  expect_equal(x, c(3, 3, 3))
  expect_error(
    get_censored_upper(
      prop_removed = c(0.5, 0.3, 0.2),
      n_catch = c(3, 3),
      n_hooks = c(5, 5, 5),
      pstar = 0.9
    ),
    regexp = "length"
  )
  expect_error(
    get_censored_upper(
      prop_removed = c(0.5, 0.3, 0.2),
      n_catch = c(3, 3, 3),
      n_hooks = c(5, 5),
      pstar = 0.9
    ),
    regexp = "length"
  )
  expect_error(
    get_censored_upper(
      prop_removed = 1.1,
      n_catch = 1,
      n_hooks = 1
    ),
    regexp = "1"
  )
  expect_error(
    get_censored_upper(
      prop_removed = -0.1,
      n_catch = 1,
      n_hooks = 1
    ),
    regexp = "0"
  )
  expect_error(
    get_censored_upper(
      prop_removed = 0.5,
      n_catch = -1,
      n_hooks = 1
    ),
    regexp = "0"
  )
  expect_error(
    get_censored_upper(
      prop_removed = 0.5,
      n_catch = 0,
      n_hooks = -1
    ),
    regexp = "0"
  )
  expect_error(
    get_censored_upper(
      prop_removed = 0.5,
      n_catch = NA_integer_,
      n_hooks = -1
    ),
    regexp = "missing"
  )
  expect_error(
    get_censored_upper(
      prop_removed = 0.5,
      n_catch = 1,
      n_hooks = NA_integer_
    ),
    regexp = "missing"
  )
})

test_that("Binomial simulation/residuals works with weights argument or cbind()", {
  set.seed(1)
  w <- sample(1:9, size = 300, replace = TRUE)
  dat <- data.frame(y = stats::rbinom(300, size = w, 0.5))
  dat$prop <- dat$y / w
  m <- sdmTMB(prop ~ 1, data = dat, weights = w, family = binomial(), spatial = "off")
  r <- residuals(m)
  expect_true(sum(is.infinite(r)) == 0L)
  set.seed(1)
  stats::qqnorm(r)
  stats::qqline(r)
  set.seed(1)
  s <- simulate(m, nsim = 500)
  expect_equal(mean(dat$y), mean(s), tolerance = 0.1)
  expect_equal(mean(apply(s, 1, mean) - dat$y), 0, tolerance = 0.01)

  # cbind() approach:
  dat$y0 <- w - dat$y
  m2 <- sdmTMB(cbind(y, y0) ~ 1, data = dat, family = binomial(), spatial = "off")

  expect_equal(m$model$par, m2$model$par)

  set.seed(1)
  r2 <- residuals(m2)
  expect_true(sum(is.infinite(r2)) == 0L)
  stats::qqnorm(r2)
  stats::qqline(r2)

  set.seed(1)
  s2 <- simulate(m2, nsim = 500)
  expect_equal(mean(dat$y), mean(s2), tolerance = 0.1)
  expect_equal(mean(apply(s2, 1, mean) - dat$y), 0, tolerance = 0.01)
})

test_that("Generalized gamma works", {
  skip_on_cran()
  d <- subset(pcod_2011, density > 0)
  fit1 <- sdmTMB(
    density ~ 1 + depth_scaled,
    data = d,
    spatial = "off",
    family = lognormal(link = "log")
  )
  fit2 <- sdmTMB(
    density ~ 1 + depth_scaled,
    data = d,
    spatial = "off",
    family = Gamma(link = "log")
  )
  fit3 <- sdmTMB(
    density ~ 1 + depth_scaled,
    data = d,
    spatial = "off",
    family = gengamma(link = "log")
  )
  expect_s3_class(fit2, "sdmTMB")
  logLik(fit1)
  logLik(fit2)
  logLik(fit3)
  get_df <- function(x) {
    L <- logLik(x)
    attr(L, "df")
  }
  df1 <- get_df(fit1)
  df2 <- get_df(fit2)
  df3 <- get_df(fit3)
  expect_identical(df1, 3L)
  expect_identical(df3, 4L)

  b <- as.list(fit3$sd_report, "Estimate")
  expect_equal(b$gengamma_Q, 0.04212623, tolerance = 0.001)
  AIC(fit1)
  AIC(fit2)
  AIC(fit3)

})


test_that("Generalized gamma matches Gamma when Q = sigma", {
  skip_on_cran()

  # Generate values drawn from generaliased gamma distribution given the mean of those values
  rgengamma <- function(n, mean, sigma, Q) {
    # Get mu from mean
    k <- Q^-2
    beta <- Q / sigma
    log_theta <- log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k )
    mu <- log_theta + log(k) / beta

    if (Q != 0) {
      w <- log(Q^2 * rgamma(n, 1 / Q^(2), 1)) / Q
      y <- exp(mu + (sigma * w))

    } else {
      y <- rlnorm(n, mu, sigma)
    }
    return(y)
  }

  sigma <- 0.5
  Q <- sigma
  mean <- 5
  n <- 10000

  # Regression coefficients (effects)
  intercept <- 1
  b1 <- 1.8
  # Generate covariate values
  set.seed(1)
  x <- runif(n, min = 0, max = 2)

  # Compute mu's
  coefs_true <- matrix(c(intercept, b1))
  X <- matrix(cbind(1, x), ncol = 2)
  y_mean <- exp(X %*% coefs_true)

  set.seed(10)
  y <- rgengamma(n = n, mean = y_mean, sigma = sigma, Q = Q)
  # Should get the same answers with flexsurv::rgengamma
  # set.seed(10)
  # y_flex <- flexsurv::rgengamma(n = n, mu = y_mu, sigma = sigma, Q = Q)

  d <- data.frame(x = x, y = y)

  fit1 <- sdmTMB(
    y ~ x,
    data = d,
    spatial = "off",
    family = gengamma(link = "log")
  )

  fit2 <- sdmTMB(
    y ~ x,
    data = d,
    spatial = "off",
    family = Gamma(link = "log")
  )

  b <- as.list(fit1$sd_report, "Estimate")
  expect_equal(b$gengamma_Q, 0.5, tolerance = 0.1)
  expect_equal(b$b_j[1], 1, tolerance = 0.01)
  expect_equal(b$b_j[2], 1.8, tolerance = 0.01)

})
