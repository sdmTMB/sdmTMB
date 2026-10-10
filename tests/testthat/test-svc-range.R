test_that("range_group_labels adds SVC rows only when needed", {
  labels <- function(g, svc = c("a", "b"), spatial = "on",
                     spatiotemporal = "iid", omit = FALSE, n_m = 1L) {
    range_group_labels(n_m, rep_len(spatial, n_m), rep_len(spatiotemporal, n_m),
      rep(TRUE, n_m), g, svc = svc, omit_spatial_intercept = omit)
  }
  # Unnamed SVCs use the spatial range: no extra rows
  x <- labels(NULL)
  expect_identical(dim(x), c(2L, 1L))
  expect_identical(svc_kappa_rows(x, 2L), matrix(0L, 2L, 1L))
  x <- labels(c(spatial = "s", a = "s"))
  expect_identical(nrow(x), 2L)
  # User labels never match generated defaults
  x <- labels(c(a = ".spatial1"))
  expect_identical(nrow(x), 4L)
  x <- labels(c(a = "default:spatial1"))
  expect_identical(nrow(x), 4L)

  # An SVC with its own range
  x <- labels(c(a = "z"))
  expect_identical(rownames(x), c("spatial", "spatiotemporal", "a", "b"))
  expect_identical(get_kappa_map(x), factor(c(1, 1, 2, 1)))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(2L, 0L), 2L))

  # Two SVCs sharing a range use the first row with that label
  x <- labels(c(a = "z", b = "z"))
  expect_identical(get_kappa_map(x), factor(c(1, 1, 2, 2)))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(2L, 2L), 2L))

  # An SVC sharing the spatiotemporal range
  x <- labels(c(spatial = "s", spatiotemporal = "st", a = "st"))
  expect_identical(get_kappa_map(x), factor(c(1, 2, 2, 1)))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(1L, 0L), 2L))

  # An SVC sharing the other component's spatial range
  x <- labels(list(c(spatial = "s1"), c(spatial = "s2", a = "s1")), n_m = 2L)
  expect_identical(get_kappa_map(x),
    factor(c(1, 1, 1, 1, 2, 2, 1, 2)))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(0L, 0L, 2L, 0L), 2L))

  # spatial = "off" with SVCs: the spatial row is on only if an SVC uses it
  x <- labels(NULL, spatiotemporal = "off", omit = TRUE)
  expect_identical(attr(x, "on")[, 1], c(spatial = TRUE, spatiotemporal = FALSE))
  x <- labels(c(a = "z", b = "z"), spatiotemporal = "off", omit = TRUE)
  expect_identical(get_kappa_map(x), factor(c(NA, NA, 1, 1)))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(2L, 2L), 2L))
  x <- labels(c(a = "z", b = "z"), omit = TRUE)
  expect_identical(attr(x, "on")[1:2, 1],
    c(spatial = FALSE, spatiotemporal = TRUE))
  expect_identical(get_kappa_map(x), factor(c(1, 1, 2, 2)))

  # SVC fields enter every component, even where the spatial field is off,
  # so that component's spatial range is kept for the SVC that uses it
  x <- labels(list(c(a = "z"), c(a = "z")), spatial = c("on", "off"), n_m = 2L)
  expect_identical(x[, 2L], c(spatial = "default:spatial2",
    spatiotemporal = "default:spatial2", a = "user:z", b = "default:spatial2"))
  expect_identical(attr(x, "on")[1L, 2L], c(spatial = TRUE))
  expect_identical(svc_kappa_rows(x, 2L), matrix(c(2L, 0L, 2L, 0L), 2L))

  expect_error(labels(c(foo = "z")), "Valid names")
  expect_error(labels(c(a = "z"), svc = c("spatial", "a")), "named `spatial`")
})

svc_range_build <- function(..., data = pcod_2011, mesh = NULL) {
  if (is.null(mesh)) mesh <- make_mesh(data, c("X", "Y"), cutoff = 20)
  sdmTMB(density ~ 1, data = data, mesh = mesh, family = tweedie(),
    spatial_varying = ~ 0 + depth_scaled + depth_scaled2,
    control = sdmTMBcontrol(backend = "rtmb"), do_fit = FALSE, ...)
}

# RTMB objective of a do_fit = FALSE model at `ln_kappa`, with other
# parameters at fixed values
svc_range_nll <- function(fit, ln_kappa) {
  d <- fit$tmb_data
  d$normalize_in_r <- 0L
  p <- fit$tmb_params
  p$ln_tau_O[] <- 0.5
  p$ln_tau_Z[] <- c(1, 1.5)
  p$ln_tau_E[] <- 0.8
  p$ln_kappa[] <- ln_kappa
  p$zeta_s[] <- seq(-1, 1, length.out = length(p$zeta_s))
  p$omega_s[] <- seq(1, -1, length.out = length(p$omega_s))
  p$ln_H_input[] <- c(0.2, -0.1)
  obj <- make_sdmTMB_adfun(d, p, fit$tmb_map, random = NULL, backend = "rtmb")
  obj$fn(obj$par)
}

test_that("SVC ranges reproduce the shared range when equal", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  for (anisotropy in c(FALSE, TRUE)) {
    f0 <- svc_range_build(mesh = mesh, anisotropy = anisotropy)
    expect_identical(dim(f0$tmb_params$ln_kappa), c(2L, 1L))
    expect_identical(f0$tmb_data$svc_kappa_row, matrix(0L, 2L, 1L))
    f1 <- svc_range_build(mesh = mesh, anisotropy = anisotropy,
      range_groups = c(depth_scaled2 = "z"))
    expect_identical(dim(f1$tmb_params$ln_kappa), c(4L, 1L))
    expect_identical(f1$tmb_map$ln_kappa, factor(c(1, 1, 1, 2)))
    expect_equal(svc_range_nll(f1, c(-2, -2, -2, -2)),
      svc_range_nll(f0, c(-2, -2)))
    expect_false(isTRUE(all.equal(svc_range_nll(f1, c(-2, -2, -2, -1)),
      svc_range_nll(f0, c(-2, -2)))))
  }
})

test_that("SVC ranges can share with other fields and each other", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  build <- function(g) {
    svc_range_build(mesh = mesh, time = "year", range_groups = g)
  }
  own <- build(c(spatial = "s", spatiotemporal = "st", depth_scaled = "a",
    depth_scaled2 = "b"))
  expect_identical(own$tmb_map$ln_kappa, factor(c(1, 2, 3, 4)))
  # Tied to the spatiotemporal range
  tied <- build(c(spatial = "s", spatiotemporal = "st", depth_scaled = "st"))
  expect_identical(tied$tmb_map$ln_kappa, factor(c(1, 2, 2, 1)))
  expect_identical(tied$tmb_data$svc_kappa_row, matrix(c(1L, 0L), 2L))
  expect_equal(svc_range_nll(tied, c(-2, -1.5, -1.5, -2)),
    svc_range_nll(own, c(-2, -1.5, -1.5, -2)))
  # Two SVCs sharing one range build one precision matrix
  shared <- build(c(depth_scaled = "z", depth_scaled2 = "z"))
  expect_identical(shared$tmb_map$ln_kappa, factor(c(1, 1, 2, 2)))
  expect_identical(shared$tmb_data$svc_kappa_row, matrix(c(2L, 2L), 2L))
  expect_equal(svc_range_nll(shared, c(-2, -2, -1, -1)),
    svc_range_nll(own, c(-2, -2, -1, -1)))
})

test_that("SVC ranges can be shared across delta components", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = delta_gamma(), spatial_varying = ~ 0 + depth_scaled,
    range_groups = list(c(spatial = "s1"), c(spatial = "s2", depth_scaled = "s1")),
    control = sdmTMBcontrol(backend = "rtmb"), do_fit = FALSE)
  expect_identical(fit$tmb_map$ln_kappa, factor(c(1, 1, 1, 2, 2, 1)))
  expect_identical(fit$tmb_data$svc_kappa_row, matrix(c(0L, 2L), 1L))
})

test_that("SVC ranges require the RTMB backend", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  expect_error(sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    spatial_varying = ~ 0 + depth_scaled, family = tweedie(),
    range_groups = c(depth_scaled = "z"), do_fit = FALSE,
    control = sdmTMBcontrol(backend = "tmb")), "RTMB backend")
  f <- svc_range_build(mesh = mesh, range_groups = c(depth_scaled = "z"))
  expect_error(make_sdmTMB_adfun(f$tmb_data, f$tmb_params, f$tmb_map,
    random = f$tmb_random), "RTMB backend")
  # Default SVC fits still work on TMB
  f <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    spatial_varying = ~ 0 + depth_scaled, family = tweedie(), do_fit = FALSE)
  expect_identical(dim(f$tmb_params$ln_kappa), c(2L, 1L))
})

test_that("SVC range fits work with post-fit methods", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = tweedie(), spatial_varying = ~ 0 + depth_scaled,
    spatial = "off", range_groups = c(depth_scaled = "z"),
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_identical(fit$tmb_map$ln_kappa, factor(c(NA, NA, 1)))
  r <- fit$tmb_obj$report()
  expect_identical(dim(r$range), c(2L, 1L))
  expect_identical(dim(r$range_Z), c(1L, 1L))
  expect_equal(r$range_Z[1, 1],
    sqrt(8) / exp(fit$tmb_obj$env$parList()$ln_kappa[3, 1]))
  expect_true("log_range_Z" %in% names(fit$sd_report$value))
  expect_equal(r$sigma_Z[1, 1], exp(-fit$model$par[["ln_tau_Z"]] -
    fit$model$par[["ln_kappa"]]) / sqrt(4 * pi))

  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
    family = tweedie(), spatial_varying = ~ 0 + depth_scaled,
    spatiotemporal = "off", range_groups = c(depth_scaled = "z"),
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_length(fit$model$par[names(fit$model$par) == "ln_kappa"], 2L)
  expect_s3_class(tidy(fit, "ran_pars"), "data.frame")
  expect_output(print(fit), "depth_scaled")
  nd <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  expect_s3_class(get_index(fit, newdata = nd), "data.frame")
  s <- simulate(fit, nsim = 2)
  expect_identical(dim(s), c(nrow(pcod_2011), 2L))
  # spread_sims() uses each field's own range
  set.seed(1)
  sims <- spread_sims(fit, nsim = 300)
  expect_false("ln_kappa" %in% names(sims))
  r <- fit$tmb_obj$report()
  expect_equal(median(sims$sigma_Z), r$sigma_Z[1, 1], tolerance = 0.2)
})

test_that("tidy() and print() show SVC ranges and shared ranges", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
    family = tweedie(), spatial_varying = ~ 0 + depth_scaled + depth_scaled2,
    range_groups = c(spatial = "s", spatiotemporal = "st",
      depth_scaled = "st", depth_scaled2 = "z"),
    control = sdmTMBcontrol(backend = "rtmb"))
  expect_identical(fit$range_groups[, 1],
    c(spatial = "user:s", spatiotemporal = "user:st", depth_scaled = "user:st",
      depth_scaled2 = "user:z"))
  b <- tidy(fit, "ran_pars")
  r <- fit$tmb_obj$report()
  expect_equal(b$estimate[b$term == "range"], r$range[, 1])
  # One row per SVC, in the order of the sigma_Z rows, even when shared
  expect_equal(b$estimate[b$term == "range_Z"], r$range_Z[, 1])
  expect_equal(r$range_Z[1, 1], r$range[2, 1])
  expect_true(all(b$conf.low[b$term == "range_Z"] < r$range_Z[, 1]))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "range \\(spatiotemporal\\): [0-9.]+ \\(shared with depth_scaled\\)")
  expect_match(out, "range \\(depth_scaled\\): [0-9.]+ \\(shared with spatiotemporal\\)")
  expect_match(out, "range \\(depth_scaled2\\): [0-9.]+\n")

  # No spatial or spatiotemporal field: only the SVC range is shown
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = tweedie(), spatial_varying = ~ 0 + depth_scaled, spatial = "off",
    range_groups = c(depth_scaled = "z"),
    control = sdmTMBcontrol(backend = "rtmb"))
  b <- tidy(fit, "ran_pars")
  expect_false("range" %in% b$term)
  expect_identical(sum(b$term == "range_Z"), 1L)
  out <- capture.output(print(fit))
  expect_identical(grep("range", out, value = TRUE),
    paste0("Matérn range (depth_scaled): ",
      mround(fit$tmb_obj$report()$range_Z[1, 1], 2L)))

  # Default SVC fits show no SVC ranges
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = tweedie(), spatial_varying = ~ 0 + depth_scaled)
  expect_false("range_Z" %in% tidy(fit, "ran_pars")$term)
  expect_length(grep("range", capture.output(print(fit))), 1L)
})

test_that("print() notes ranges shared across delta components", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = delta_gamma(), range_groups = c(spatial = "a"))
  out <- capture.output(print(fit))
  expect_length(grep("range: [0-9.]+ \\(shared with model 2 spatial\\)$", out), 1L)
  expect_length(grep("range: [0-9.]+ \\(shared with model 1 spatial\\)$", out), 1L)
  expect_identical(fit$range_groups, structure(
    matrix("user:a", 2L, 2L, dimnames = list(c("spatial", "spatiotemporal"), NULL)),
    on = matrix(c(TRUE, FALSE), 2L, 2L,
      dimnames = list(c("spatial", "spatiotemporal"), NULL))))
})

test_that("matern_svc PC priors apply to SVC fields once per range group", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  pc <- pc_matern(range_gt = 10, sigma_lt = 5)
  flags <- function(priors, ...) {
    d <- svc_range_build(mesh = mesh, priors = priors, ...)$tmb_data
    list(sigma = as.vector(d$sigma_prior), range = as.vector(d$range_prior))
  }
  # Rows: spatial, spatiotemporal, depth_scaled, depth_scaled2
  expect_identical(flags(sdmTMBpriors(matern_svc = pc)),
    list(sigma = c(0L, 0L, 1L, 1L), range = c(0L, 0L, 1L, 0L)))
  expect_identical(flags(sdmTMBpriors(matern_s = pc, matern_svc = pc)),
    list(sigma = c(1L, 0L, 1L, 1L), range = c(1L, 0L, 0L, 0L)))
  expect_identical(flags(sdmTMBpriors(matern_s = pc, matern_svc = pc),
    range_groups = c(depth_scaled2 = "z")),
    list(sigma = c(1L, 0L, 1L, 1L), range = c(1L, 0L, 0L, 1L)))
  # With spatial = "off", matern_svc rather than matern_s sets the SVC range
  expect_identical(flags(sdmTMBpriors(matern_s = pc, matern_svc = pc),
    spatial = "off"), list(sigma = c(0L, 0L, 1L, 1L), range = c(0L, 0L, 1L, 0L)))
  expect_identical(flags(sdmTMBpriors(matern_s = pc), spatial = "off"),
    list(sigma = c(0L, 0L, 0L, 0L), range = c(1L, 0L, 0L, 0L)))
  # One prior for every component
  d <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    family = delta_gamma(), spatial_varying = ~ 0 + depth_scaled,
    priors = sdmTMBpriors(matern_svc = pc),
    control = sdmTMBcontrol(backend = "rtmb"), do_fit = FALSE)$tmb_data
  expect_identical(d$sigma_prior, matrix(c(0L, 0L, 1L), 3L, 2L))
  expect_identical(d$range_prior, matrix(c(0L, 0L, 1L), 3L, 2L))

  # Objective adds the expected terms, including with the Stan Jacobian
  pc_svc <- pc_matern(range_gt = 5, sigma_lt = 2)
  for (bayesian in c(FALSE, TRUE)) {
    nll <- vapply(list(sdmTMBpriors(matern_s = pc, matern_svc = pc_svc),
      sdmTMBpriors()), function(priors) {
        # The spatial and depth_scaled range priors differ (warns)
        f <- suppressWarnings(svc_range_build(mesh = mesh, priors = priors,
          bayesian = bayesian, range_groups = c(depth_scaled2 = "z")))
        svc_range_nll(f, c(-2, -2, -2, -1))
      }, numeric(1L))
    expected <- -rtmb_pc_matern(0.5, -2, pc, stan = bayesian) -
      rtmb_pc_matern(1, -2, pc_svc, include_range = FALSE, stan = bayesian) -
      rtmb_pc_matern(1.5, -1, pc_svc, stan = bayesian)
    expect_equal(nll[[1]] - nll[[2]], expected, tolerance = 1e-8)
  }
})

test_that("matern_svc requires RTMB and old fits get no SVC prior", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  pc <- pc_matern(range_gt = 10, sigma_lt = 5)
  expect_error(sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
    spatial_varying = ~ 0 + depth_scaled, family = tweedie(),
    priors = sdmTMBpriors(matern_svc = pc), do_fit = FALSE,
    control = sdmTMBcontrol(backend = "tmb")), "RTMB backend")
  # Data from before `matern_svc`: shorter priors vector, two flag rows
  f <- svc_range_build(mesh = mesh, priors = sdmTMBpriors(matern_s = pc))
  nll <- svc_range_nll(f, c(-2, -2))
  f$tmb_data$priors <- utils::head(f$tmb_data$priors, -4L)
  f$tmb_data$sigma_prior <- f$tmb_data$sigma_prior[1:2, , drop = FALSE]
  f$tmb_data$range_prior <- f$tmb_data$range_prior[1:2, , drop = FALSE]
  expect_equal(svc_range_nll(f, c(-2, -2)), nll)
})

test_that("spread_sims() returns each SVC's SD and range draws", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit_svc <- function(...) {
    sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, family = tweedie(),
      control = sdmTMBcontrol(backend = "rtmb"), ...)
  }
  # Expected transformations of the same underlying draws
  raw_draws <- function(fit, nsim) {
    set.seed(1)
    s <- rmvnorm_prec(fit$tmb_obj$env$last.par.best, fit$sd_report, nsim)
    pn <- names(c(fit$sd_report$par.fixed, fit$sd_report$par.random))
    s <- s[pn %in% names(fit$sd_report$par.fixed), , drop = FALSE]
    list(ln_kappa = s[rownames(s) == "ln_kappa", , drop = FALSE],
      ln_tau_Z = s[rownames(s) == "ln_tau_Z", , drop = FALSE])
  }
  sims <- function(fit, nsim = 5) {
    set.seed(1)
    spread_sims(fit, nsim = nsim)
  }
  sd_z <- function(ln_tau, ln_kappa) 1 / sqrt(4 * pi * exp(2 * ln_tau + 2 * ln_kappa))

  # One SVC with its own range and no other fields
  fit <- fit_svc(spatial_varying = ~ 0 + depth_scaled, spatial = "off",
    range_groups = c(depth_scaled = "z"))
  x <- sims(fit)
  r <- raw_draws(fit, 5)
  expect_false("range" %in% names(x))
  expect_false(any(grepl("ln_", names(x))))
  expect_equal(x$range_Z, sqrt(8) / exp(r$ln_kappa[1, ]), ignore_attr = TRUE)
  expect_equal(x$sigma_Z, sd_z(r$ln_tau_Z[1, ], r$ln_kappa[1, ]),
    ignore_attr = TRUE)
  g <- gather_sims(fit, nsim = 5)
  expect_setequal(unique(g$.variable), setdiff(names(x), ".iteration"))

  # Two SVCs: one shares the spatial range, one has its own
  fit <- fit_svc(spatial_varying = ~ 0 + depth_scaled + depth_scaled2,
    range_groups = c(depth_scaled2 = "z"))
  x <- sims(fit)
  r <- raw_draws(fit, 5)
  expect_false(any(grepl("ln_", names(x))))
  expect_equal(x$range, sqrt(8) / exp(r$ln_kappa[1, ]), ignore_attr = TRUE)
  expect_equal(x$range_Z_depth_scaled, x$range)
  expect_equal(x$range_Z_depth_scaled2, sqrt(8) / exp(r$ln_kappa[2, ]),
    ignore_attr = TRUE)
  expect_equal(x$sigma_Z_depth_scaled, sd_z(r$ln_tau_Z[1, ], r$ln_kappa[1, ]),
    ignore_attr = TRUE)
  expect_equal(x$sigma_Z_depth_scaled2, sd_z(r$ln_tau_Z[2, ], r$ln_kappa[2, ]),
    ignore_attr = TRUE)

  # Two SVCs sharing one range
  fit <- fit_svc(spatial_varying = ~ 0 + depth_scaled + depth_scaled2,
    spatial = "off", range_groups = c(depth_scaled = "z", depth_scaled2 = "z"))
  x <- sims(fit)
  expect_false("range" %in% names(x))
  expect_equal(x$range_Z_depth_scaled, x$range_Z_depth_scaled2)

  # Default single-SVC output is unchanged
  fit <- fit_svc(spatial_varying = ~ 0 + depth_scaled)
  x <- sims(fit)
  r <- raw_draws(fit, 5)
  expect_true(all(c("range", "sigma_O", "sigma_Z") %in% names(x)))
  expect_false(any(grepl("range_Z|ln_", names(x))))
  expect_equal(x$sigma_Z, sd_z(r$ln_tau_Z[1, ], r$ln_kappa[1, ]),
    ignore_attr = TRUE)
})
