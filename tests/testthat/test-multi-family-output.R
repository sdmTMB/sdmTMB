test_that("Multi-family residuals are blocked", {
  dat <- data.frame(
    y = c(1, 0, 2, 0),
    dist = c("poisson", "binomial", "poisson", "binomial")
  )
  fam <- list(
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
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  expect_error(
    residuals(fit),
    regexp = "multi-family"
  )
})

test_that("Multi-family simulation respects model = 1 vs model = 2", {
  skip_on_cran()

  dat <- data.frame(
    y = c(0, 2, 0, 3, 1, 0, 5),
    dist = "delta_gamma"
  )
  fam <- list(
    delta_gamma = delta_gamma()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist"
  )

  set.seed(101)
  sims1 <- predict(fit, newdata = dat, nsim = 5, model = 1, type = "link")
  set.seed(101)
  sims2 <- predict(fit, newdata = dat, nsim = 5, model = 2, type = "link")

  expect_true(is.matrix(sims1))
  expect_true(is.matrix(sims2))
  expect_equal(dim(sims1), dim(sims2))
  expect_true(any(abs(sims1 - sims2) > 1e-10))
})

test_that("Multi-family delta simulations return combined responses", {
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
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  sims <- simulate(fit, nsim = 3, silent = TRUE)
  expect_true(is.matrix(sims))
  expect_equal(nrow(sims), nrow(dat))
  expect_equal(ncol(sims), 3L)
  expect_true(rlang::is_integerish(sims[dat$dist == "poisson",1L]))
})

test_that("Multi-family predict rejects NAs in prediction columns", {
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    x = c(0.1, 0.5, 0.2, 0.2, 0.3, 0.7),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  fit <- sdmTMB(
    y ~ x,
    data = dat,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
  nd <- dat
  nd$x[2] <- NA_real_
  expect_error(
    predict(fit, newdata = nd),
    regexp = "NAs are not allowed in multi-family prediction columns"
  )
})

test_that("tidy works for multi-family models", {
  dat <- data.frame(
    y = c(1, 2.5, 1, 3.2),
    dist = c("poisson", "gaussian", "poisson", "gaussian")
  )
  fam <- list(
    poisson = poisson(),
    gaussian = gaussian()
  )
  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    spatial = "off",
    family = fam,
    distribution_column = "dist"
  )
  td <- tidy(fit, effects = "fixed")
  expect_true("(Intercept)" %in% td$term)

  td_re <- tidy(fit, effects = "ran_pars")
  expect_true("phi_gaussian" %in% td_re$term)
})

test_that("sdmTMB_cv works for multi-family models", {
  skip_on_cran()
  skip_on_ci()

  set.seed(110)
  n <- 60
  dist <- rep(c("gaussian", "poisson"), each = n / 2)
  x <- stats::runif(n, -1, 1)
  eta <- 0.2 + 0.5 * x
  y <- numeric(n)
  idx_gauss <- dist == "gaussian"
  y[idx_gauss] <- eta[idx_gauss] + stats::rnorm(sum(idx_gauss), sd = 0.2)
  y[!idx_gauss] <- stats::rpois(sum(!idx_gauss), lambda = exp(eta[!idx_gauss]))
  dat <- data.frame(
    y = y,
    x = x,
    dist = dist,
    X = stats::runif(n),
    Y = stats::runif(n)
  )
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 0.5)

  fam <- list(
    gaussian = gaussian(),
    poisson = poisson()
  )

  cv <- sdmTMB_cv(
    y ~ 1 + x,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = fam,
    distribution_column = "dist",
    k_folds = 2,
    parallel = FALSE
  )

  expect_s3_class(cv, "sdmTMB_cv")
  expect_true(all(is.finite(cv$data$cv_predicted)))
  expect_true(all(is.finite(cv$data$cv_loglik)))

  fold1 <- cv$data[cv$data$cv_fold == 1, , drop = FALSE]
  pred_link <- predict(cv$models[[1]], newdata = fold1, type = "link")
  manual <- pred_link$est
  manual[fold1$dist == "poisson"] <- exp(manual[fold1$dist == "poisson"])
  expect_equal(fold1$cv_predicted, manual, tolerance = 1e-6)
})

test_that("Multi-family se_fit rejects response scale", {
  skip_on_cran()

  dat <- data.frame(
    y = c(1, 2, 0, 3, 1, 0),
    dist = rep(c("poisson", "gaussian"), 3)
  )
  fam <- list(
    poisson = poisson(),
    gaussian = gaussian()
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
  expect_error(
    predict(fit, newdata = dat, type = "response", se_fit = TRUE),
    regexp = "se_fit"
  )
})

test_that("SE works for single-family multi-family model (delta_gamma)", {
  skip_on_cran()

  # Create data with single family in multi-family format
  dat <- pcod
  dat$data_type <- "delta_gamma"

  # Fit multi-family model with only delta_gamma
  family_list <- list(delta_gamma = delta_gamma())

  fit_mf <- sdmTMB(
    density ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off"
  )

  # Fit regular delta model for comparison
  fit_delta <- sdmTMB(
    density ~ depth_scaled,
    data = pcod,
    family = delta_gamma(),
    spatial = "off"
  )

  # Prediction data
  nd <- data.frame(
    depth_scaled = seq(min(pcod$depth_scaled), max(pcod$depth_scaled), length.out = 30),
    data_type = "delta_gamma"
  )

  # Test: multi-family predictions with se_fit should work
  p_mf <- predict(fit_mf, newdata = nd, re_form = NA, se_fit = TRUE, type = "link")

  expect_true("est" %in% names(p_mf))
  expect_true("est_se" %in% names(p_mf))
  expect_true(all(!is.na(p_mf$est_se)))
  expect_true(all(p_mf$est_se > 0))

  # Test: combined estimate should match manual calculation
  expect_equal(
    as.numeric(p_mf$est),
    as.numeric(log(plogis(p_mf$est1) * exp(p_mf$est2))),
    tolerance = 1e-6
  )

  # Test: multi-family should match regular delta model
  nd_delta <- nd
  nd_delta$data_type <- NULL
  p_delta <- predict(fit_delta, newdata = nd_delta, re_form = NA, se_fit = TRUE, type = "link")

  expect_equal(p_mf$est, p_delta$est, tolerance = 1e-4)
  expect_equal(p_mf$est_se, p_delta$est_se, tolerance = 1e-4)
})

test_that("SE works for multi-family combined predictions (non-pop)", {
  skip_on_cran()

  dat <- pcod
  dat$data_type <- "delta_gamma"

  family_list <- list(delta_gamma = delta_gamma())

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off"
  )

  nd <- data.frame(
    depth_scaled = seq(min(dat$depth_scaled), max(dat$depth_scaled), length.out = 30),
    data_type = "delta_gamma"
  )

  p <- predict(fit, newdata = nd, se_fit = TRUE, type = "link")

  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  expect_true(all(!is.na(p$est_se)))
  expect_true(all(p$est_se > 0))

  # Combined estimate should match manual calculation
  expect_equal(
    as.numeric(p$est),
    as.numeric(log(plogis(p$est1) * exp(p$est2))),
    tolerance = 1e-6
  )
  expect_true(any(abs(p$est - p$est1) > 1e-6))
})

test_that("SE works for mixed delta and non-delta families", {
  skip_on_cran()

  set.seed(456)

  # Create data with mixed families
  dat <- pcod
  n <- nrow(dat)
  dat$data_type <- sample(c("delta_gamma", "gaussian"), size = n, replace = TRUE, prob = c(0.6, 0.4))

  # For gaussian rows, use a transformed response to avoid negative values
  dat$response_mf <- ifelse(dat$data_type == "gaussian",
                             log(dat$density + 1),
                             dat$density)

  family_list <- list(
    delta_gamma = delta_gamma(),
    gaussian = gaussian()
  )

  fit <- sdmTMB(
    response_mf ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off"
  )

  # Prediction data with mixed families
  nd <- data.frame(
    depth_scaled = rep(seq(min(dat$depth_scaled), max(dat$depth_scaled), length.out = 20), 2),
    data_type = rep(c("delta_gamma", "gaussian"), each = 20)
  )

  # Test: se_fit should work without error
  p <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link")

  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  expect_true(all(!is.na(p$est_se)))
  expect_true(all(p$est_se > 0))

  # Test: delta rows should have est1 and est2
  delta_rows <- which(nd$data_type == "delta_gamma")
  if (length(delta_rows) > 0) {
    expect_true(all(!is.na(p$est1[delta_rows])))
    expect_true(all(!is.na(p$est2[delta_rows])))
    # Combined est should differ from components
    expect_true(any(abs(p$est[delta_rows] - p$est1[delta_rows]) > 1e-6))
  }

  # Test: gaussian rows should have est = est1 (no second component)
  gaussian_rows <- which(nd$data_type == "gaussian")
  if (length(gaussian_rows) > 0) {
    expect_equal(as.numeric(p$est[gaussian_rows]), as.numeric(p$est1[gaussian_rows]), tolerance = 1e-6)
    expect_true(all(is.na(p$est2[gaussian_rows])))
  }
})

test_that("SE works for population predictions (pop_pred = TRUE)", {
  skip_on_cran()

  set.seed(789)
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)

  dat <- pcod
  dat$data_type <- "delta_lognormal"

  family_list <- list(delta_lognormal = delta_lognormal())

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    mesh = mesh,
    spatial = "on"
  )

  nd <- data.frame(
    depth_scaled = seq(min(dat$depth_scaled), max(dat$depth_scaled), length.out = 25),
    data_type = "delta_lognormal"
  )
  nd$X <- mean(dat$X)
  nd$Y <- mean(dat$Y)

  # Test: population-level predictions with se_fit
  p <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link")

  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  expect_true(all(!is.na(p$est_se)))
  expect_true(all(p$est_se > 0))

  # Verify combined transformation
  expect_equal(
    as.numeric(p$est),
    as.numeric(log(plogis(p$est1) * exp(p$est2))),
    tolerance = 1e-6
  )
})

test_that("SE works for Poisson-link delta families", {
  skip_on_cran()

  dat <- pcod
  dat$data_type <- "delta_pl"

  family_list <- list(delta_pl = delta_gamma(type = "poisson-link"))

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off"
  )

  nd <- data.frame(
    depth_scaled = seq(min(dat$depth_scaled), max(dat$depth_scaled), length.out = 25),
    data_type = "delta_pl"
  )

  # Test: Poisson-link delta predictions with se_fit
  p <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link")

  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  expect_true(all(!is.na(p$est_se)))
  expect_true(all(p$est_se > 0))

  # Verify Poisson-link transformation (additive in log space)
  expect_equal(
    as.numeric(p$est),
    as.numeric(p$est1 + p$est2),
    tolerance = 1e-6
  )
})

test_that("Component-specific SEs work with multi-family", {
  skip_on_cran()

  set.seed(202)

  dat <- pcod
  dat$data_type <- "delta_gamma"

  family_list <- list(delta_gamma = delta_gamma())

  fit <- sdmTMB(
    density ~ depth_scaled,
    data = dat,
    family = family_list,
    distribution_column = "data_type",
    spatial = "off"
  )

  nd <- data.frame(
    depth_scaled = seq(min(dat$depth_scaled), max(dat$depth_scaled), length.out = 20),
    data_type = "delta_gamma"
  )

  # Test: model = 1 should give component 1 SE
  p1 <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link", model = 1)
  expect_true("est_se" %in% names(p1))
  expect_equal(p1$est, p1$est1, tolerance = 1e-6)

  # Test: model = 2 should give component 2 SE
  p2 <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link", model = 2)
  expect_true("est_se" %in% names(p2))
  expect_equal(p2$est, p2$est2, tolerance = 1e-6)

  # Component SEs should differ from combined SE
  p_combined <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, type = "link", model = NA)
  expect_true(any(abs(p1$est_se - p_combined$est_se) > 1e-6))
  expect_true(any(abs(p2$est_se - p_combined$est_se) > 1e-6))
})
