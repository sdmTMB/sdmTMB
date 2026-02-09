make_dl_predict_data <- function() {
  set.seed(101)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 3) + y / max(y)))
  eta <- 0.2 + 0.5 * x1 - 0.3 * x2
  data.frame(
    y = eta + rnorm(length(eta), sd = 0.1),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

make_dl_predict_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

make_dl_predict_delta_data <- function() {
  set.seed(202)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year / 2) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 3) + y / max(y)))
  positive <- exp(0.2 + 0.4 * x1)
  present <- as.integer(x2 > 0)
  data.frame(
    y = ifelse(present == 1L, positive, 0),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

test_that("distributed lag predict works for default and newdata pathways", {
  skip_on_cran()

  dat <- make_dl_predict_data()
  mesh <- make_dl_predict_mesh(dat)

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1) + time(x2),
    control = sdmTMBcontrol(newton_loops = 0)
  )

  p_fit <- predict(fit)
  p_new <- predict(fit, newdata = dat)
  expect_equal(nrow(p_fit), nrow(dat))
  expect_equal(p_fit$est, p_new$est, tolerance = 1e-6)
  expect_equal(p_fit$est_non_rf, p_new$est_non_rf, tolerance = 1e-6)
  expect_true(all(c("dl_cov_space_x1", "dl_cov_time_x2") %in% names(p_fit)))
  expect_equal(p_fit$dl_cov_space_x1, p_new$dl_cov_space_x1, tolerance = 1e-6)
  expect_equal(p_fit$dl_cov_time_x2, p_new$dl_cov_time_x2, tolerance = 1e-6)

  p_se <- predict(fit, newdata = dat, re_form = NA, se_fit = TRUE)
  expect_true("est_se" %in% names(p_se))
  expect_true(all(is.finite(p_se$est_se)))

  sims <- predict(fit, newdata = dat, nsim = 3)
  expect_equal(dim(sims), c(nrow(dat), 3L))
})

test_that("distributed lag newdata covariate changes prediction direction when lag beta is fixed positive", {
  skip_on_cran()

  dat <- make_dl_predict_data()
  mesh <- make_dl_predict_mesh(dat)

  proto <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ time(x1),
    do_fit = FALSE
  )

  lag_col <- proto$distributed_lags_data$term_coef_name
  lag_idx <- match(lag_col, colnames(proto$tmb_data$X_ij[[1]]))
  b_map <- seq_along(proto$tmb_params$b_j)
  b_map[lag_idx] <- NA_integer_
  b_start <- proto$tmb_params$b_j
  b_start[lag_idx] <- 1

  fit <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ time(x1),
    control = sdmTMBcontrol(
      start = list(b_j = b_start),
      map = list(b_j = factor(b_map)),
      newton_loops = 0,
      getsd = FALSE
    )
  )

  nd_low <- dat
  nd_high <- dat
  nd_high$x1 <- nd_high$x1 + 0.5

  p_low <- predict(fit, newdata = nd_low)
  p_high <- predict(fit, newdata = nd_high)

  expect_true("dl_cov_time_x1" %in% names(p_low))
  expect_gt(mean(p_high$dl_cov_time_x1 - p_low$dl_cov_time_x1), 0)
  expect_gt(mean(p_high$est - p_low$est), 0)
})

test_that("delta distributed lag in component 2 changes combined response predictions", {
  skip_on_cran()

  dat <- make_dl_predict_delta_data()
  mesh <- make_dl_predict_mesh(dat)

  proto <- sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    distributed_lags = ~ time(x1),
    do_fit = FALSE
  )

  lag_col <- proto$distributed_lags_data$term_coef_name
  lag_idx1 <- match(lag_col, colnames(proto$tmb_data$X_ij[[1]]))
  lag_idx2 <- match(lag_col, colnames(proto$tmb_data$X_ij[[2]]))

  b_map1 <- seq_along(proto$tmb_params$b_j)
  b_map1[lag_idx1] <- NA_integer_
  b_start1 <- proto$tmb_params$b_j
  b_start1[lag_idx1] <- 0

  b_map2 <- seq_along(proto$tmb_params$b_j2)
  b_map2[lag_idx2] <- NA_integer_
  b_start2 <- proto$tmb_params$b_j2
  b_start2[lag_idx2] <- 1

  fit <- suppressWarnings(sdmTMB(
    y ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = delta_gamma(),
    distributed_lags = ~ time(x1),
    control = sdmTMBcontrol(
      start = list(b_j = b_start1, b_j2 = b_start2),
      map = list(b_j = factor(b_map1), b_j2 = factor(b_map2)),
      newton_loops = 0,
      getsd = FALSE
    )
  ))

  nd_low <- dat
  nd_high <- dat
  nd_high$x1 <- nd_high$x1 + 0.5

  p_low <- predict(fit, newdata = nd_low, type = "response")
  p_high <- predict(fit, newdata = nd_high, type = "response")

  expect_gt(mean(p_high$est - p_low$est), 0)
})

test_that("temporal distributed lags require full modeled time coverage in newdata", {
  skip_on_cran()

  dat <- make_dl_predict_data()
  mesh <- make_dl_predict_mesh(dat)

  fit_time <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ time(x1),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  nd_subset <- dat[dat$year %in% sort(unique(dat$year))[1:3], , drop = FALSE]
  expect_error(
    predict(fit_time, newdata = nd_subset),
    regexp = "requires full time coverage"
  )

  fit_space <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  expect_silent({
    p_space <- predict(fit_space, newdata = nd_subset)
    expect_equal(nrow(p_space), nrow(nd_subset))
  })
})

test_that("space-only distributed lag predict works without modeled time", {
  skip_on_cran()

  dat <- make_dl_predict_data()
  dat$year <- NULL
  mesh <- make_dl_predict_mesh(dat)

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1),
    control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  p <- predict(fit)
  expect_equal(nrow(p), nrow(dat))
  expect_true("dl_cov_space_x1" %in% names(p))
})
