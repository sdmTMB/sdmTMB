# Data are simulated *through* the nonlocal diffusion/time-lag operator
# itself (via simulate_new(), following vignettes/articles/nonlocal-covariates.Rmd)
# so that kappaS_nl/rhoT are actually identifiable from the data; a plain
# regression-style toy dataset with no real diffusion signal leaves kappaS_nl
# on a flat likelihood ridge and produces NaN standard errors.
make_nl_output_data <- function(nonlocal_formula) {
  set.seed(1)
  n_t <- 6L
  n_sites <- 60L
  site <- data.frame(X = runif(n_sites), Y = runif(n_sites))
  dat <- data.frame(
    X = rep(site$X, times = n_t),
    Y = rep(site$Y, times = n_t),
    year = rep(seq_len(n_t), each = n_sites)
  )
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.15)
  grid <- expand.grid(vertex = seq_len(nrow(mesh$mesh$loc)), year = seq_len(n_t))
  grid$X <- mesh$mesh$loc[grid$vertex, 1]
  grid$Y <- mesh$mesh$loc[grid$vertex, 2]
  grid$x1 <- as.numeric(scale(
    sin(2 * pi * (grid$X + grid$year / 8)) + cos(2 * pi * (grid$Y - grid$year / 10)) +
      0.6 * sin(4 * pi * grid$X) * cos(grid$year / 3) + rnorm(nrow(grid), sd = 0.15)
  ))
  grid$x2 <- as.numeric(scale(
    cos(2 * pi * (grid$X - grid$year / 6)) + sin(2 * pi * (grid$Y + grid$year / 9)) +
      rnorm(nrow(grid), sd = 0.15)
  ))
  grid$vertex <- NULL
  dat$x1 <- NA_real_
  dat$x2 <- NA_real_
  for (tt in seq_len(n_t)) {
    rows <- dat$year == tt
    dat$x1[rows] <- as.numeric(mesh$A_st[rows, ] %*% grid$x1[grid$year == tt])
    dat$x2[rows] <- as.numeric(mesh$A_st[rows, ] %*% grid$x2[grid$year == tt])
  }
  dat$y <- 0 # placeholder; simulate_new() needs the response column to exist
  n_nl_covariates <- length(unique(all.vars(nonlocal_formula)))
  sim <- simulate_new(
    formula = y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    phi = 0.1,
    range = 0.3,
    sigma_O = 0,
    B = c(0.2, 0.4, -0.2, rep(0.5, n_nl_covariates)),
    nonlocal_formula = nonlocal_formula,
    nonlocal_data = grid,
    lags_kappaS = 4.4,
    lags_rhoT = 0.3,
    seed = 123
  )
  dat$y <- sim$observed
  list(dat = dat, mesh = mesh, grid = grid)
}

test_that("covariate diffusion fixed effects are named consistently in tidy/coef/vcov", {
  skip_on_cran()
  d <- make_nl_output_data(~ diffusion(x1) + time_lag(x2))

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = d$dat,
    mesh = d$mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
    nonlocal_data = d$grid,
    control = sdmTMBcontrol(newton_loops = 1)
  )

  lag_terms <- fit$nonlocal_parsed$term_coef_name
  td <- tidy(fit, effects = "fixed", silent = TRUE)
  expect_true(all(lag_terms %in% td$term))
  expect_equal(length(unique(td$term)), nrow(td))

  cf <- coef(fit)
  expect_true(all(lag_terms %in% names(cf)))

  vc <- vcov(fit)
  expect_true(all(lag_terms %in% rownames(vc)))
  expect_true(all(lag_terms %in% colnames(vc)))
})

test_that("covariate diffusion ran_pars include lag scales and derived diagnostics", {
  skip_on_cran()
  d <- make_nl_output_data(~ diffusion(x1) + time_lag(x1))

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = d$dat,
    mesh = d$mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
    nonlocal_data = d$grid,
    control = sdmTMBcontrol(newton_loops = 1)
  )

  td <- tidy(fit, effects = "ran_pars", silent = TRUE)
  expected_terms <- c(
    "kappaS_nl[x1]",
    "kappaT_nl[x1]",
    "rhoT[x1]",
    "MSDK[x1]",
    "RMSDK[x1]"
  )
  expect_true(all(expected_terms %in% td$term), info = paste(setdiff(expected_terms, td$term), collapse = ", "))

  rep_est <- as.list(fit$sd_report, "Estimate", report = TRUE)
  rep_se <- as.list(fit$sd_report, "Std. Error", report = TRUE)
  expect_length(rep_est$kappaS_nl, 1L)
  expect_length(rep_est$kappaT_nl, 1L)
  expect_length(rep_est$rhoT, 1L)
  expect_length(rep_est$MSDK, 1L)
  expect_length(rep_est$RMSDK, 1L)
  expect_length(rep_est$log_MSDK, 1L)
  expect_length(rep_est$log_RMSDK, 1L)
  expect_length(rep_se$kappaS_nl, 1L)
  expect_length(rep_se$kappaT_nl, 1L)
  expect_length(rep_se$rhoT, 1L)
  expect_length(rep_se$MSDK, 1L)
  expect_length(rep_se$RMSDK, 1L)
  expect_length(rep_se$log_MSDK, 1L)
  expect_length(rep_se$log_RMSDK, 1L)
  expect_gt(rep_est$MSDK, 0)
  expect_gt(rep_est$RMSDK, 0)
  expect_true(all(c(rep_se$MSDK, rep_se$RMSDK) >= 0))
  expect_null(rep_est[["MSD", exact = TRUE]])
  expect_null(rep_est[["RMSD", exact = TRUE]])
  expect_equal(
    rep_est$MSDK,
    4 / (rep_est$kappaS_nl^2 * (1 + rep_est$kappaT_nl)),
    tolerance = 1e-6
  )
  expect_equal(rep_est$RMSDK^2, rep_est$MSDK, tolerance = 1e-6)

  td_msd <- td[td$term == "MSDK[x1]", , drop = FALSE]
  td_rmsd <- td[td$term == "RMSDK[x1]", , drop = FALSE]
  expect_gt(td_msd$conf.low, 0)
  expect_gt(td_rmsd$conf.low, 0)
  expect_equal(
    td_msd$conf.low,
    as.numeric(exp(rep_est$log_MSDK - stats::qnorm(0.975) * rep_se$log_MSDK)),
    tolerance = 1e-6
  )
  expect_equal(
    td_rmsd$conf.low,
    as.numeric(exp(rep_est$log_RMSDK - stats::qnorm(0.975) * rep_se$log_RMSDK)),
    tolerance = 1e-6
  )
})

test_that("print output reports covariate diffusion structure and diagnostics", {
  skip_on_cran()
  d <- make_nl_output_data(~ diffusion(x1) + time_lag(x1))

  fit <- sdmTMB(
    y ~ x1 + x2,
    data = d$dat,
    mesh = d$mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
    nonlocal_data = d$grid,
    control = sdmTMBcontrol(newton_loops = 1)
  )

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Nonlocal formula: diffusion\\(x1\\) \\+ time_lag\\(x1\\)")
  expect_match(out, "rhoT\\[x1\\]=")
  expect_match(out, "RMSDK\\[x1\\]=")
})
