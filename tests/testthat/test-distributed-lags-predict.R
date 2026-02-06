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

  expect_gt(mean(p_high$est - p_low$est), 0)
})
