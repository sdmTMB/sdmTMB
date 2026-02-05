test_that("Stage 1 multi-family builds TMB data", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    do_fit = FALSE
  )
  expect_identical(fit$tmb_data$multi_family, 1L)
  expect_equal(fit$tmb_data$e_i, c(0L, 1L, 2L, 0L, 1L, 2L))
  expect_equal(
    fit$tmb_data$family_e1,
    c(
      sdmTMB:::.enum_family("gaussian"),
      sdmTMB:::.enum_family("poisson"),
      sdmTMB:::.enum_family("binomial")
    )
  )
  expect_equal(
    fit$tmb_data$link_e1,
    c(
      sdmTMB:::.enum_link("identity"),
      sdmTMB:::.enum_link("log"),
      sdmTMB:::.enum_link("logit")
    )
  )
  expect_equal(fit$tmb_data$family_e2, rep(-1L, 3))
  expect_equal(fit$tmb_data$link_e2, rep(-1L, 3))
  expect_equal(fit$tmb_data$delta_family_e, rep(0L, 3))
})

test_that("Stage 1 multi-family can fit a simple model", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist"
  )
  expect_s3_class(fit, "sdmTMB")
  expect_true(is.list(fit$family))
})

test_that("Multi-family fit rejects NAs in modeling columns", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    x = c(0.1, 0.5, NA, 0.2, 0.3, 0.7),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  expect_error(
    sdmTMB(
      y ~ x,
      data = dat,
      spatial = "off",
      spatiotemporal = "off",
      family = fam,
      distribution_column = "dist",
      do_fit = FALSE
    ),
    regexp = "NAs are not allowed in multi-family modeling columns"
  )
})

test_that("Stage 1 multi-family binomial proportions use weights as size", {
  dat <- data.frame(
    y = c(0.2, 0.7, 1, 0, 0.5, 0.3),
    dist = c("binomial", "binomial", "gaussian", "gaussian", "binomial", "binomial"),
    w = c(10, 8, 2, 3, 12, 6)
  )
  fam <- list(
    gaussian = gaussian(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    weights = dat$w,
    do_fit = FALSE
  )
  prop_rows <- dat$dist == "binomial" & dat$y > 0 & dat$y < 1
  expect_equal(
    fit$tmb_data$size[prop_rows],
    dat$w[prop_rows]
  )
  expect_equal(
    fit$tmb_data$y_i[prop_rows, 1],
    dat$y[prop_rows] * dat$w[prop_rows]
  )
  expect_equal(
    fit$tmb_data$weights_i[prop_rows],
    rep(1, sum(prop_rows))
  )
  binary_rows <- dat$dist == "binomial" & dat$y %in% c(0, 1)
  expect_equal(
    fit$tmb_data$size[binary_rows],
    rep(1, sum(binary_rows))
  )
  expect_equal(
    fit$tmb_data$weights_i[binary_rows],
    dat$w[binary_rows]
  )
})

test_that("Stage 2 multi-family supports extra-parameter families", {
  dat <- data.frame(
    y = c(0, 2, 4, 0.5, 1.1, 2.2, 0, 3.5),
    dist = c(
      "nbinom2", "nbinom2", "nbinom2",
      "Gamma", "Gamma", "Gamma",
      "tweedie", "tweedie"
    )
  )
  fam <- list(
    nbinom2 = nbinom2(),
    Gamma = Gamma(link = "log"),
    tweedie = tweedie(link = "log")
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$ln_phi_len, c(1L, 1L, 1L))
  expect_equal(fit$tmb_data$ln_phi_start, c(0L, 1L, 2L))
  expect_equal(fit$tmb_data$thetaf_len, c(0L, 0L, 1L))
  expect_equal(fit$tmb_data$thetaf_start, c(-1L, -1L, 0L))
})

test_that("Stage 2 multi-family supports Beta and betabinomial", {
  dat <- data.frame(
    y = c(0.2, 0.5, 0.8, 0.1, 0.4, 0.7),
    dist = c(
      "Beta", "Beta", "Beta",
      "betabinomial", "betabinomial", "betabinomial"
    ),
    w = c(1, 1, 1, 10, 12, 8)
  )
  fam <- list(
    Beta = Beta(),
    betabinomial = betabinomial()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    weights = dat$w,
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$ln_phi_len, c(1L, 1L))
  expect_equal(fit$tmb_data$ln_phi_start, c(0L, 1L))
  expect_equal(
    fit$tmb_data$size[dat$dist == "betabinomial"],
    dat$w[dat$dist == "betabinomial"]
  )
  expect_true(all(fit$tmb_data$weights_i[dat$dist == "betabinomial"] == 1))
  expect_equal(
    fit$tmb_data$y_i[dat$dist == "betabinomial", 1],
    dat$y[dat$dist == "betabinomial"] * dat$w[dat$dist == "betabinomial"]
  )
})

test_that("Stage 3 multi-family supports delta families", {
  dat <- data.frame(
    y = c(0, 2, 0, 3, 1, 0, 5),
    dist = c(
      "delta_gamma", "delta_gamma", "delta_gamma", "delta_gamma",
      "poisson", "poisson", "poisson"
    )
  )
  fam <- list(
    delta_gamma = delta_gamma(),
    poisson = poisson()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_s3_class(fit, "sdmTMB")
  expect_equal(fit$tmb_data$delta_family_e, c(1L, 0L))
  expect_equal(
    fit$tmb_data$family_e1,
    c(sdmTMB:::.enum_family("binomial"), sdmTMB:::.enum_family("poisson"))
  )
  expect_equal(
    fit$tmb_data$family_e2,
    c(sdmTMB:::.enum_family("Gamma"), -1L)
  )
  expect_equal(ncol(fit$tmb_data$y_i), 2L)
  expect_true(all(is.na(fit$tmb_data$y_i[dat$dist == "poisson", 2])))
})

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
