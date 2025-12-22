test_that("Multi-family recovers shared linear predictors and dispersion params", {
  set.seed(123)

  n_each <- 3000
  beta0 <- 0.3
  beta1 <- 0.7
  sigma <- 0.2
  gamma_shape <- 8
  binom_size <- 40
  bb_size <- 30
  bb_phi <- 25
  nb_phi <- 12
  student_phi <- 0.4
  student_df <- 6

  families <- c(
    "gaussian", "poisson", "binomial", "betabinomial",
    "nbinom2", "student", "delta_gamma_pl"
  )
  dist <- rep(families, each = n_each)
  n <- length(dist)
  x <- stats::runif(n, -1, 1)
  eta <- beta0 + beta1 * x

  y <- numeric(n)
  idx <- dist == "gaussian"
  y[idx] <- eta[idx] + stats::rnorm(sum(idx), sd = sigma)
  idx <- dist == "poisson"
  y[idx] <- stats::rpois(sum(idx), lambda = exp(eta[idx]))
  idx <- dist == "binomial"
  p <- stats::plogis(eta[idx])
  y[idx] <- stats::rbinom(sum(idx), size = binom_size, prob = p) / binom_size
  idx <- dist == "betabinomial"
  mu <- stats::plogis(eta[idx])
  p <- stats::rbeta(sum(idx), shape1 = mu * bb_phi, shape2 = (1 - mu) * bb_phi)
  y[idx] <- stats::rbinom(sum(idx), size = bb_size, prob = p)
  idx <- dist == "nbinom2"
  y[idx] <- stats::rnbinom(sum(idx), size = nb_phi, mu = exp(eta[idx]))
  idx <- dist == "student"
  y[idx] <- eta[idx] + student_phi * stats::rt(sum(idx), df = student_df)
  idx <- dist == "delta_gamma_pl"
  n_val <- exp(eta[idx])
  p <- 1 - exp(-n_val)
  z <- stats::rbinom(sum(idx), size = 1, prob = p)
  w <- exp(eta[idx])
  r <- (n_val * w) / p
  y_pos <- stats::rgamma(sum(idx), shape = gamma_shape, scale = r / gamma_shape)
  y[idx] <- ifelse(z == 1, y_pos, 0)

  weights <- rep(1, n)
  weights[dist == "binomial"] <- binom_size
  weights[dist == "betabinomial"] <- bb_size

  dat <- data.frame(
    y = y,
    x = x,
    dist = dist,
    weights = weights
  )

  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial(),
    betabinomial = betabinomial(),
    nbinom2 = nbinom2(),
    student = student(),
    delta_gamma_pl = delta_gamma(type = "poisson-link")
  )

  fit <- sdmTMB(
    formula = list(y ~ 1 + x, y ~ 1 + x),
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    weights = dat$weights
  )

  truth <- c("(Intercept)" = beta0, "x" = beta1)
  tol <- 0.02

  fixed_1 <- tidy(fit, "fixed", model = 1)
  est_1 <- fixed_1$estimate[match(names(truth), fixed_1$term)]
  expect_equal(as.numeric(est_1), unname(truth), tolerance = tol)

  fixed_2 <- tidy(fit, "fixed", model = 2)
  est_2 <- fixed_2$estimate[match(names(truth), fixed_2$term)]
  expect_equal(as.numeric(est_2), unname(truth), tolerance = tol)

  ran_1 <- tidy(fit, "ran_pars", model = 1)
  phi_bb <- ran_1$estimate[match("phi_betabinomial", ran_1$term)]
  phi_nb <- ran_1$estimate[match("phi_nbinom2", ran_1$term)]
  phi_student <- ran_1$estimate[match("phi_student", ran_1$term)]
  df_student <- ran_1$estimate[match("student_df_student", ran_1$term)]
  expect_equal(as.numeric(phi_bb), bb_phi, tolerance = 0.1)
  expect_equal(as.numeric(phi_nb), nb_phi, tolerance = 0.1)
  expect_equal(as.numeric(phi_student), student_phi, tolerance = 0.1)
  expect_equal(as.numeric(df_student), student_df, tolerance = 0.1)

  ran_2 <- tidy(fit, "ran_pars", model = 2)
  phi_dg <- ran_2$estimate[match("phi_delta_gamma_pl", ran_2$term)]
  expect_equal(as.numeric(phi_dg), gamma_shape, tolerance = 0.1)
})
