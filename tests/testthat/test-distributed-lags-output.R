make_dl_output_data <- function() {
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

make_dl_output_mesh <- function(dat) {
  make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.5)
}

test_that("distributed lag fixed effects are named consistently in tidy/coef/vcov", {
  skip_on_cran()
  dat <- make_dl_output_data()
  mesh <- make_dl_output_mesh(dat)

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1) + time(x2),
    control = sdmTMBcontrol(newton_loops = 0)
  ))

  lag_terms <- fit$distributed_lags_data$term_coef_name
  td <- tidy(fit, effects = "fixed", silent = TRUE)
  expect_true(all(lag_terms %in% td$term))
  expect_equal(length(unique(td$term)), nrow(td))

  cf <- coef(fit)
  expect_true(all(lag_terms %in% names(cf)))

  vc <- vcov(fit)
  expect_true(all(lag_terms %in% rownames(vc)))
  expect_true(all(lag_terms %in% colnames(vc)))
})

test_that("distributed lag ran_pars include lag scales and derived diagnostics", {
  skip_on_cran()
  dat <- make_dl_output_data()
  mesh <- make_dl_output_mesh(dat)

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1) + time(x1) + spacetime(x2),
    control = sdmTMBcontrol(newton_loops = 0)
  ))

  td <- tidy(fit, effects = "ran_pars", silent = TRUE)
  expected_terms <- c("kappaS_dl", "kappaT_dl", "kappaST_dl", "rhoT", "MSD", "RMSD")
  expect_true(all(expected_terms %in% td$term))
  for (nm in expected_terms) {
    expect_equal(sum(td$term == nm), 1L)
  }
})

test_that("print output reports distributed lag structure and diagnostics", {
  skip_on_cran()
  dat <- make_dl_output_data()
  mesh <- make_dl_output_mesh(dat)

  fit <- suppressWarnings(sdmTMB(
    y ~ x1 + x2,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    distributed_lags = ~ space(x1) + time(x1) + spacetime(x2),
    control = sdmTMBcontrol(newton_loops = 0)
  ))

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Distributed lags: space\\(x1\\) \\+ time\\(x1\\) \\+ spacetime\\(x2\\)")
  expect_match(out, "Distributed lag temporal persistence \\(rhoT\\):")
  expect_match(out, "Distributed lag root mean square displacement \\(RMSD\\):")
})
