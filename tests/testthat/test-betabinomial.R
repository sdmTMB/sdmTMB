test_that("Beta-binomial family returns correct structure", {
  .names <- c("family", "link", "linkfun", "linkinv")
  expect_true(all(.names %in% names(betabinomial(link = "logit"))))
  expect_true(all(.names %in% names(betabinomial(link = "cloglog"))))
})

test_that("Beta-binomial family works with appropriate links", {
  expect_identical(class(betabinomial(link = "logit")), "family")
  expect_identical(class(betabinomial(link = logit)), "family")
  expect_identical(class(betabinomial(link = "cloglog")), "family")
  expect_identical(class(betabinomial(link = cloglog)), "family")
  expect_error(class(betabinomial(link = "probit")))
  expect_error(class(betabinomial(link = "banana")))
  expect_error(class(betabinomial(link = banana)))
})


test_that("Beta-binomial fits with logit link", {
  skip_on_cran()
  set.seed(1)

  # Simulate beta-binomial data
  n_trials <- 20
  prob <- 0.3
  phi <- 2
  n_obs <- 200

  # Create spatial coordinates
  x <- stats::runif(n_obs, -1, 1)
  y <- stats::runif(n_obs, -1, 1)
  dat <- data.frame(x = x, y = y)

  # Simulate data manually since we don't have betabinomial simulation yet
  # Use rbeta then rbinom to create beta-binomial
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  dat$y <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$n_trials <- n_trials

  # Create mesh
  spde <- make_mesh(dat, c("x", "y"), n_knots = 50, type = "kmeans")

  # Fit model using cbind() syntax
  m <- sdmTMB(
    cbind(y, n_trials - y) ~ 1,
    data = dat,
    mesh = spde,
    family = betabinomial(link = "logit"),
    spatial = "off",
    spatiotemporal = "off"
  )

  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  expect_length(residuals(m), nrow(dat))
  expect_true(m$model$convergence == 0)
})

test_that("Beta-binomial fits with cloglog link", {
  skip_on_cran()
  set.seed(42)

  # Simulate smaller dataset for cloglog test
  n_trials <- 10
  prob <- 0.2
  phi <- 1.5
  n_obs <- 100

  # Create spatial coordinates
  x <- stats::runif(n_obs, -1, 1)
  y <- stats::runif(n_obs, -1, 1)
  dat <- data.frame(x = x, y = y)

  # Simulate beta-binomial data
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  dat$y <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$n_trials <- n_trials

  # Fit model
  m <- sdmTMB(
    cbind(y, n_trials - y) ~ 1,
    data = dat,
    family = betabinomial(link = "cloglog"),
    spatial = "off"
  )

  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  expect_length(residuals(m), nrow(dat))
  expect_true(m$model$convergence == 0)
})

test_that("Beta-binomial with proportion and weights works", {
  skip_on_cran()
  set.seed(123)

  # Create test data with proportions
  n_obs <- 100
  n_trials <- sample(5:20, n_obs, replace = TRUE)
  prob <- 0.4
  phi <- 3

  x <- stats::runif(n_obs, -1, 1)
  y <- stats::runif(n_obs, -1, 1)
  dat <- data.frame(x = x, y = y)

  # Simulate beta-binomial data
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  successes <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$prop <- successes / n_trials

  # Fit model using proportion and weights
  m <- sdmTMB(
    prop ~ 1,
    data = dat,
    weights = n_trials,
    family = betabinomial(link = "logit"),
    spatial = "off",
  )

  expect_true(all(!is.na(summary(m$sd_report)[, "Std. Error"])))
  expect_length(residuals(m), nrow(dat))
  expect_true(m$model$convergence == 0)
})


test_that("Beta-binomial simulation works", {
  skip_on_cran()
  set.seed(789)

  n_trials <- 10
  n_obs <- 50

  x <- stats::runif(n_obs, -1, 1)
  y <- stats::runif(n_obs, -1, 1)
  dat <- data.frame(x = x, y = y)

  # Create minimal dataset for simulation test
  beta_vals <- stats::rbeta(n_obs, 1.5, 3.5)
  dat$y <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$n_trials <- n_trials

  spde <- make_mesh(dat, c("x", "y"), n_knots = 30, type = "kmeans")

  m <- sdmTMB(
    cbind(y, n_trials - y) ~ 1,
    data = dat,
    mesh = spde,
    family = betabinomial(),
    spatial = "off",
    spatiotemporal = "off"
  )

  # Test simulation
  s <- simulate(m, nsim = 10)
  expect_true(is.matrix(s))
  expect_equal(nrow(s), nrow(dat))
  expect_equal(ncol(s), 10)

  # Check that simulated values are within valid range
  expect_true(all(s >= 0 & s <= n_trials))
})

test_that("Beta-binomial matches glmmTMB when spatial='off'", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  set.seed(1234)

  # Create test data
  n_obs <- 100
  n_trials <- 20
  prob <- 0.4
  phi <- 2.5

  x <- stats::runif(n_obs, -1, 1)
  dat <- data.frame(
    x = x,
    y = 1:n_obs,  # dummy y for mesh
    covariate = x
  )

  # Simulate beta-binomial data
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  dat$successes <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$failures <- n_trials - dat$successes

  # Fit with glmmTMB
  m_glmmTMB <- glmmTMB::glmmTMB(
    cbind(successes, failures) ~ covariate,
    data = dat,
    family = glmmTMB::betabinomial()
  )

  # Fit with sdmTMB (no spatial effects)
  m_sdmTMB <- sdmTMB(
    cbind(successes, failures) ~ covariate,
    data = dat,
    family = betabinomial(),
    spatial = "off"
  )

  # Extract coefficients
  glmmTMB_coefs <- fixef(m_glmmTMB)$cond
  sdmTMB_coefs <- tidy(m_sdmTMB, effects = "fixed")$estimate

  # Compare fixed effects (should be very close)
  expect_equal(as.numeric(glmmTMB_coefs), sdmTMB_coefs, tolerance = 0.0001)

  # Extract dispersion parameters
  glmmTMB_phi <- exp(fixef(m_glmmTMB)$disp)
  sdmTMB_phi <- exp(m_sdmTMB$model$par[["ln_phi"]])

  # Compare dispersion parameters
  expect_equal(as.numeric(glmmTMB_phi), sdmTMB_phi, tolerance = 0.001)

  # Compare log-likelihoods (should be very close)
  expect_equal(logLik(m_glmmTMB)[1], logLik(m_sdmTMB)[1], tolerance = 0.001)
})

test_that("Beta-binomial handles NA in response", {
  skip_on_cran()
  set.seed(1)

  n_trials <- 10
  n_obs <- 100

  dat <- data.frame(id = seq_len(n_obs))

  prob <- 0.4
  phi <- 2
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  dat$successes <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$failures <- n_trials - dat$successes

  # Test 1: cbind syntax - both successes and failures NA (same rows)
  dat$successes[c(5, 10, 15)] <- NA
  dat$failures[c(5, 10, 15)] <- NA

  m <- sdmTMB(
    cbind(successes, failures) ~ 1,
    data = dat,
    family = betabinomial(),
    spatial = "off"
  )

  expect_true(m$model$convergence == 0)

  # Test 2: cbind syntax - only failures is NA (triggers the original bug)
  dat2 <- data.frame(id = seq_len(n_obs))
  dat2$successes <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat2$failures <- n_trials - dat2$successes
  dat2$failures[c(5, 10, 15)] <- NA

  m2 <- sdmTMB(
    cbind(successes, failures) ~ 1,
    data = dat2,
    family = betabinomial(),
    spatial = "off"
  )

  expect_true(m2$model$convergence == 0)

  # Test 3: proportions + weights syntax - NA in proportion
  dat3 <- data.frame(id = 1:n_obs)
  dat3$prop <- beta_vals
  dat3$size <- n_trials
  dat3$prop[c(5, 10, 15)] <- NA

  m3 <- sdmTMB(
    prop ~ 1,
    data = dat3,
    family = betabinomial(),
    weights = dat3$size,
    spatial = "off"
  )

  expect_true(m3$model$convergence == 0)

  # try same with pre-filtering to make sure it matches exactly
  dat3 <- dat3[!is.na(dat3$prop),,drop=FALSE]
  m3filtered <- sdmTMB(
    prop ~ 1,
    data = dat3,
    family = betabinomial(),
    weights = dat3$size,
    spatial = "off"
  )
  expect_equal(logLik(m3), logLik(m3filtered))

  # now with offset
  dat3 <- data.frame(id = 1:n_obs)
  dat3$prop <- beta_vals
  dat3$size <- n_trials
  dat3$offset <- rnorm(n_obs)
  dat3$prop[c(5, 10, 15)] <- NA

  expect_error({m3 <- sdmTMB(
    prop ~ 1,
    data = dat3,
    family = betabinomial(),
    offset = dat3$offset,
    spatial = "off"
  )}, regexp = "weights")

  m3 <- sdmTMB(
    prop ~ 1,
    data = dat3,
    family = betabinomial(),
    weights = dat3$size,
    offset = dat3$offset,
    spatial = "off"
  )
  dat3 <- dat3[!is.na(dat3$prop),,drop=FALSE]
  m3filtered <- sdmTMB(
    prop ~ 1,
    data = dat3,
    family = betabinomial(),
    weights = dat3$size,
    offset = dat3$offset,
    spatial = "off"
  )
  expect_equal(logLik(m3), logLik(m3filtered))

  # Test 4: proportions + weights syntax - NA in weights (size)
  dat4 <- data.frame(id = 1:n_obs)
  dat4$prop <- beta_vals
  dat4$size <- rep(n_trials, n_obs)
  dat4$size[c(5, 10, 15)] <- NA

  m4 <- sdmTMB(
    prop ~ 1,
    data = dat4,
    family = betabinomial(),
    weights = dat4$size,
    spatial = "off"
  )

  expect_true(m4$model$convergence == 0)
})

test_that("Betabinomial can be used with index standardization", {
  skip_on_cran()
  skip_on_ci()
  set.seed(1)
  n_trials <- 10
  n_obs <- 100
  dat <- data.frame(id = seq_len(n_obs))
  prob <- 0.4
  phi <- 2
  beta_vals <- stats::rbeta(n_obs, prob * phi, (1 - prob) * phi)
  dat$successes <- stats::rbinom(n_obs, size = n_trials, prob = beta_vals)
  dat$trials <- n_trials
  dat$prop <- dat$successes / dat$trials
  dat$time <- 1L
  m <- sdmTMB(
    prop ~ 1,
    data = dat,
    family = betabinomial(),
    weights = dat$trials,
    spatial = "off",
    spatiotemporal = "off",
    time = "time"
  )
  nd <- data.frame(time = 1L)
  p <- predict(m, newdata = nd, return_tmb_object = TRUE)
  i <- get_index(p, area = 10) # 10 'trials'
  expect_equal(i$est, stats::plogis(p$data$est) * 10)
  expect_equal(i$est, 3.452362, tolerance = 1e-5)
})
