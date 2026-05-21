test_that("ordbeta() family returns correct structure", {
  .names <- c("family", "link", "linkfun", "linkinv", "initialize")
  expect_true(all(.names %in% names(ordbeta(link = "logit"))))
  expect_identical(class(ordbeta()), "family")
  expect_identical(ordbeta()$family, "ordbeta")
  expect_identical(ordbeta()$link, "logit")
})

test_that("ordbeta() fit recovers parameters", {
  skip_on_cran()
  set.seed(1)
  n <- 500
  x <- rnorm(n)
  eta <- 0.3 + 0.8 * x
  mu <- plogis(eta)
  phi <- 6
  psi <- c(-1.2, 1.0)
  p0 <- plogis(psi[1] - eta)
  p1 <- plogis(eta - psi[2])
  u <- runif(n)
  y <- numeric(n)
  zero <- u < p0
  one <- u > 1 - p1
  mid <- !zero & !one
  y[zero] <- 0
  y[one] <- 1
  if (any(mid)) y[mid] <- rbeta(sum(mid), mu[mid] * phi, (1 - mu[mid]) * phi)
  dat <- data.frame(y = y, x = x)

  fit <- sdmTMB(y ~ x, data = dat, family = ordbeta(), spatial = "off")
  expect_s3_class(fit, "sdmTMB")

  b <- tidy(fit)
  expect_equal(b$estimate[b$term == "(Intercept)"], 0.3, tolerance = 0.15)
  expect_equal(b$estimate[b$term == "x"], 0.8, tolerance = 0.15)

  r <- tidy(fit, "ran_pars")
  expect_true("phi" %in% r$term)
  expect_true("ordbeta_cutpoint_lower" %in% r$term)
  expect_true("ordbeta_cutpoint_upper" %in% r$term)
  expect_equal(r$estimate[r$term == "ordbeta_cutpoint_lower"],
    plogis(psi[1]), tolerance = 0.25)
  expect_equal(r$estimate[r$term == "ordbeta_cutpoint_upper"],
    plogis(psi[2]), tolerance = 0.25)
  expect_equal(r$estimate[r$term == "phi"], phi, tolerance = 0.5)
})

test_that("ordbeta() validates [0, 1] response", {
  expect_error(
    sdmTMB(y ~ 1,
      data = data.frame(y = c(0.5, 1.1, 0.2)),
      family = ordbeta(), spatial = "off"),
    regexp = "\\[0, 1\\]"
  )
})

test_that("ordbeta() residuals are well calibrated", {
  skip_on_cran()
  set.seed(2)
  n <- 400
  x <- rnorm(n)
  eta <- 0.1 + 0.6 * x
  mu <- plogis(eta)
  phi <- 5
  psi <- c(-1.0, 1.2)
  p0 <- plogis(psi[1] - eta); p1 <- plogis(eta - psi[2])
  u <- runif(n); y <- numeric(n)
  zero <- u < p0; one <- u > 1 - p1; mid <- !zero & !one
  y[zero] <- 0; y[one] <- 1
  if (any(mid)) y[mid] <- rbeta(sum(mid), mu[mid] * phi, (1 - mu[mid]) * phi)
  fit <- sdmTMB(y ~ x, data = data.frame(y = y, x = x),
    family = ordbeta(), spatial = "off")
  r <- residuals(fit)
  expect_equal(mean(r), 0, tolerance = 0.15)
  expect_equal(sd(r), 1, tolerance = 0.15)
})

test_that("ordbeta() likelihood matches glmmTMB exactly", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(42)
  n <- 600
  x <- rnorm(n)
  eta <- 0.4 + 0.7 * x
  mu <- plogis(eta)
  phi <- 5
  psi <- c(-1.3, 0.9)
  p0 <- plogis(psi[1] - eta); p1 <- plogis(eta - psi[2])
  u <- runif(n); y <- numeric(n)
  zero <- u < p0; one <- u > 1 - p1; mid <- !zero & !one
  y[zero] <- 0; y[one] <- 1
  if (any(mid)) y[mid] <- rbeta(sum(mid), mu[mid] * phi, (1 - mu[mid]) * phi)
  dat <- data.frame(y = y, x = x)

  f_s <- sdmTMB(y ~ x, data = dat, family = ordbeta(), spatial = "off")
  f_g <- glmmTMB::glmmTMB(y ~ x, data = dat, family = glmmTMB::ordbeta())

  expect_equal(as.numeric(logLik(f_s)), as.numeric(logLik(f_g)),
    tolerance = 1e-5)

  b_s <- tidy(f_s)
  b_g <- summary(f_g)$coef$cond
  expect_equal(b_s$estimate[b_s$term == "(Intercept)"],
    unname(b_g["(Intercept)", "Estimate"]), tolerance = 1e-4)
  expect_equal(b_s$estimate[b_s$term == "x"],
    unname(b_g["x", "Estimate"]), tolerance = 1e-4)

  r_s <- tidy(f_s, "ran_pars")
  expect_equal(r_s$estimate[r_s$term == "phi"], glmmTMB::sigma(f_g),
    tolerance = 1e-3)

  psi_g <- f_g$fit$parfull
  psi_g <- unname(psi_g[names(psi_g) == "psi"])
  expect_equal(r_s$estimate[r_s$term == "ordbeta_cutpoint_lower"],
    plogis(psi_g[1]), tolerance = 1e-4)
  expect_equal(r_s$estimate[r_s$term == "ordbeta_cutpoint_upper"],
    plogis(psi_g[2]), tolerance = 1e-4)
})

test_that("ordbeta() is rejected in multi-family mode", {
  expect_error(
    sdmTMB(y ~ 1,
      data = data.frame(y = c(0, 0.5, 1), g = "a"),
      family = list(ordbeta = ordbeta(), normal = gaussian()),
      distribution_column = "g", spatial = "off"),
    regexp = "ordbeta"
  )
})
