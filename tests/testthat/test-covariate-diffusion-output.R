make_nl_output_data <- function() {
  set.seed(101)
  n_t <- 5L
  n_s <- 6L
  year <- rep(seq_len(n_t), each = n_s)
  x <- rep(seq_len(n_s), times = n_t)
  y <- rep(1:2, length.out = n_t * n_s)
  x1 <- as.numeric(scale(sin(year) + x / max(x)))
  x2 <- as.numeric(scale(cos(year / 2) + y / max(y)))
  eta <- 0.2 + 0.5 * x1 - 0.1 * x2
  data.frame(
    y = eta + rnorm(length(eta), sd = 0.12),
    x1 = x1,
    x2 = x2,
    year = year,
    X = x,
    Y = y
  )
}

make_nl_output_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

make_nl_output_grid <- function(mesh, years) {
  make_nl_covariate_grid(mesh, years, c("x1", "x2"))
}

test_that("covariate diffusion fixed effects are named consistently in tidy/coef/vcov", {
  skip_on_cran()
  dat <- make_nl_output_data()
  mesh <- make_nl_output_mesh(dat)
  grid <- make_nl_output_grid(mesh, sort(unique(dat$year)))

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x2),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0)
  ))

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
  dat <- make_nl_output_data()
  mesh <- make_nl_output_mesh(dat)
  grid <- make_nl_output_grid(mesh, sort(unique(dat$year)))

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0)
  ))

  td <- tidy(fit, effects = "ran_pars", silent = TRUE)
  expected_terms <- c(
    "kappaS_nl[x1]",
    "kappaT_nl[x1]",
    "rhoT[x1]",
    "MSD[x1]",
    "RMSD[x1]"
  )
  expect_true(all(expected_terms %in% td$term), info = paste(setdiff(expected_terms, td$term), collapse = ", "))

  rep_est <- as.list(fit$sd_report, "Estimate", report = TRUE)
  rep_se <- as.list(fit$sd_report, "Std. Error", report = TRUE)
  expect_length(rep_est$kappaS_nl, 1L)
  expect_length(rep_est$kappaT_nl, 1L)
  expect_length(rep_est$rhoT, 1L)
  expect_length(rep_est$MSD, 1L)
  expect_length(rep_est$RMSD, 1L)
  expect_length(rep_se$kappaS_nl, 1L)
  expect_length(rep_se$kappaT_nl, 1L)
  expect_length(rep_se$rhoT, 1L)
  expect_length(rep_se$MSD, 1L)
  expect_length(rep_se$RMSD, 1L)
})

test_that("print output reports covariate diffusion structure and diagnostics", {
  skip_on_cran()
  dat <- make_nl_output_data()
  mesh <- make_nl_output_mesh(dat)
  grid <- make_nl_output_grid(mesh, sort(unique(dat$year)))

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
    nonlocal_data = grid,
    control = sdmTMBcontrol(newton_loops = 0)
  ))

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Nonlocal formula: diffusion\\(x1\\) \\+ time_lag\\(x1\\)")
  expect_match(out, "rhoT\\[x1\\]=")
  expect_match(out, "RMSD\\[x1\\]=")
})
