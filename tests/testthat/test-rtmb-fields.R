rtmb_spatial_data <- function(n = 90L, seed = 72) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 5), y = runif(n, 0, 5), z = rnorm(n),
    time = rep(1:3, each = n / 3))
  d$response <- ifelse(runif(n) < 0.6, rgamma(n, 2, 1), 0)
  d$gaussian <- d$z + sin(d$x) + rnorm(n, sd = 0.3)
  d
}

# A ring of ten areas plus one isolated area, with a smooth areal effect.
rtmb_areal_data <- function() {
  g <- igraph::make_ring(10L)
  g <- igraph::add_vertices(g, 1L)
  igraph::V(g)$name <- paste0("u", seq_len(igraph::vcount(g)))
  domain <- suppressWarnings(make_areal_domain(g, space_column = "unit"))
  d <- expand.grid(unit = domain$unit_names, time = 1:3,
    stringsAsFactors = FALSE)
  set.seed(73)
  d$z <- rnorm(nrow(d))
  effect <- c(sin(2 * pi * (1:10) / 10), 0.5)
  d$y <- rpois(nrow(d),
    exp(1 + 0.2 * d$z + effect[match(d$unit, domain$unit_names)]))
  list(data = d, domain = domain)
}

test_that("RTMB isotropic SPDE matches frozen TMB values", {
  set.seed(35)
  n <- 48L
  observed <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  observed$response <- 0.4 - 0.2 * observed$z + rnorm(n, sd = 0.7)
  mesh <- make_mesh(observed, c("x", "y"), n_knots = 16L,
    type = "kmeans")
  fit <- sdmTMB(response ~ z, data = observed, mesh = mesh,
    spatial = "on", spatiotemporal = "off", family = gaussian(),
    control = sdmTMBcontrol(multiphase = FALSE), do_fit = FALSE)
  p <- fit$tmb_params
  p$b_j[] <- c(0.4, -0.2)
  p$ln_phi[] <- log(0.7)
  p$ln_tau_O[] <- 0.15
  p$ln_kappa[] <- -0.1
  p$omega_s[] <- seq(-0.1, 0.1, length.out = length(p$omega_s))
  expect_frozen_objective("fields isotropic spde", fit$tmb_data, p,
    fit$tmb_map, random = "omega_s")
})

test_that("RTMB spatial Gaussian fit predicts and simulates", {
  set.seed(36)
  n <- 80L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  d$response <- 0.3 + 0.4 * d$z + rnorm(n, sd = 0.2)
  mesh <- suppressWarnings(
    make_mesh(d, c("x", "y"), n_knots = 20L, type = "kmeans"),
    classes = "sdmTMB_mesh_outside_warning"
  )
  rt <- sdmTMB(response ~ z, data = d, mesh = mesh, spatial = "on",
    spatiotemporal = "off", family = gaussian(),
    control = sdmTMBcontrol(multiphase = FALSE))
  expect_true(rt$sd_report$pdHess)
  projection_offset <- rep(0.2, 4L)
  p <- suppressMessages(predict(rt, newdata = d[1:4, ],
    offset = projection_offset, se_fit = TRUE))
  expect_true(all(is.finite(p$est_se)))
  population <- predict(rt, newdata = d[1:4, ],
    offset = projection_offset, pop_pred = TRUE)
  expect_true(all(is.finite(population$est)))
  held <- simulate(rt, nsim = 2L, observation_error = FALSE, silent = TRUE)
  expect_equal(as.numeric(held[, 1L]),
    as.numeric(predict(rt, type = "response")$est), tolerance = 1e-6)
  draws <- simulate(rt, nsim = 2L, re_form = ~0,
    return_tmb_report = TRUE, silent = TRUE)
  expect_false(isTRUE(all.equal(draws[[1L]]$omega_s,
    draws[[2L]]$omega_s)))
  expect_equal(as.numeric(draws[[1L]]$omega_s_A),
    as.numeric(rt$tmb_data$A_st %*% draws[[1L]]$omega_s),
    tolerance = 1e-6)
})

test_that("RTMB isotropic IID, AR1, and RW fields match frozen TMB values", {
  set.seed(37)
  n <- 50L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    time = rep(1:5, each = 10L))
  d$response <- 0.3 + 0.4 * d$z + rnorm(n, sd = 0.2)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  for (field_type in c("iid", "ar1", "rw")) {
    fit <- sdmTMB(response ~ z, data = d, mesh = mesh, time = "time",
      spatial = "off", spatiotemporal = field_type, family = gaussian(),
      control = sdmTMBcontrol(multiphase = FALSE), do_fit = FALSE)
    p <- fit$tmb_params
    p$ln_tau_E[] <- 0.2
    p$ln_kappa[] <- 0.1
    p$ar1_phi[] <- 0.4
    p$epsilon_st[] <- seq(-0.12, 0.17, length.out = length(p$epsilon_st))
    expect_frozen_objective(paste("fields spatiotemporal", field_type),
      fit$tmb_data, p, fit$tmb_map, fit$tmb_random)
  }
})

test_that("RTMB temporal field simulation retains selected time steps", {
  set.seed(38)
  n <- 50L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n),
    time = rep(1:5, each = 10L))
  d$response <- rnorm(n)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 15L, type = "kmeans")
  for (field_type in c("iid", "ar1", "rw")) {
    fit <- sdmTMB(response ~ z, data = d, mesh = mesh, time = "time",
      spatial = "off", spatiotemporal = field_type, family = gaussian(),
      control = sdmTMBcontrol(multiphase = FALSE), do_fit = FALSE)
    p <- fit$tmb_params
    p$ln_tau_E[] <- 0.4
    p$ln_kappa[] <- 0.1
    p$ar1_phi[] <- 0.4
    p$epsilon_st[] <- seq(-0.1, 0.1, length.out = length(p$epsilon_st))
    data <- fit$tmb_data
    data$sim_re[[2L]] <- 1L
    data$simulate_t <- c(0L, 1L, 0L, 1L, 0L)
    obj <- make_sdmTMB_adfun(data, p, fit$tmb_map, fit$tmb_random)
    draw <- obj$simulate()
    for (t in c(1L, 3L, 5L)) {
      expect_equal(draw$epsilon_st[, t, 1L], p$epsilon_st[, t, 1L],
        info = field_type)
    }
    for (t in c(2L, 4L)) {
      expect_false(isTRUE(all.equal(draw$epsilon_st[, t, 1L],
        p$epsilon_st[, t, 1L])), info = field_type)
    }
    expect_equal(as.numeric(draw$eta_i),
      as.numeric(data$X_ij[[1L]] %*% p$b_j + data$offset_i) +
        as.numeric(draw$epsilon_st_A_vec), tolerance = 1e-7,
      info = field_type)
  }
})

test_that("RTMB anisotropic fields match frozen TMB values", {
  d <- rtmb_spatial_data()
  d$response[[5L]] <- NA
  mesh <- make_mesh(d, c("x", "y"), n_knots = 30L, type = "kmeans")
  newdata <- d[c(1, 40, 90), ]
  cases <- list(
    ar1 = list(gaussian ~ z, family = gaussian(), time = "time",
      spatiotemporal = "ar1", share_range = FALSE),
    svc = list(gaussian ~ z, family = gaussian(), spatial_varying = ~ 0 + z),
    delta = list(response ~ z, family = delta_gamma(), time = "time",
      spatiotemporal = list("iid", "off"), share_range = list(FALSE, TRUE))
  )
  for (name in names(cases)) {
    args <- c(cases[[name]], list(data = d, mesh = mesh, anisotropy = TRUE,
      do_fit = FALSE))
    fit <- do.call(sdmTMB, args)
    expect_frozen_fit_data("fields", fit, newdata, info = name)
  }

  large_mesh <- make_mesh(d, c("x", "y"), n_knots = 60L,
    type = "kmeans")
  large_fit <- sdmTMB(gaussian ~ z, data = d, mesh = large_mesh,
    anisotropy = TRUE, do_fit = FALSE)
  expect_frozen_fit_data("fields", large_fit, newdata, info = "large anisotropy")
})

test_that("RTMB areal SAR and CAR fields match frozen TMB values", {
  skip_if_not_installed("igraph")
  for (spatial_model in c("sar", "car")) {
    areal <- rtmb_areal_data()
    d <- areal$data
    d$y[[3L]] <- NA
    fit <- sdmTMB(y ~ z, data = d, mesh = areal$domain,
      spatial_model = spatial_model, family = poisson(), time = "time",
      spatiotemporal = "iid", spatial_varying = ~ 0 + z, do_fit = FALSE)
    expect_frozen_fit_data("fields", fit, d[c(1, 12, 30), ], info = spatial_model)
    fit <- sdmTMB(y ~ z, data = d, mesh = areal$domain,
      spatial_model = spatial_model,
      family = delta_gamma(type = "poisson-link"), time = "time", spatiotemporal = list("ar1", "off"), do_fit = FALSE)
    expect_frozen_fit_data("fields", fit, d[c(1, 12, 30), ],
      info = paste(spatial_model, "delta"))
  }
})

test_that("RTMB barrier fields match frozen TMB values", {
  skip_if_not_installed("sdmTMBextra")
  skip_if_not_installed("INLAspacetime")
  skip_if_not_installed("sf")
  d <- rtmb_spatial_data()
  mesh <- make_mesh(d, c("x", "y"), cutoff = 0.7)
  barrier <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_polygon(list(
    matrix(c(2.3, -0.2, 2.7, -0.2, 2.7, 5.2, 2.3, 5.2, 2.3, -0.2),
      ncol = 2, byrow = TRUE)))))
  mesh <- suppressMessages(suppressWarnings(
    sdmTMBextra::add_barrier_mesh(mesh, barrier, range_fraction = 0.2,
      plot = FALSE)))
  fit <- sdmTMB(gaussian ~ z, data = d, mesh = mesh, time = "time",
    spatiotemporal = "rw", spatial_varying = ~ 0 + z, share_range = FALSE,
    do_fit = FALSE)
  expect_frozen_fit_data("fields", fit, d[c(1, 40, 90), ], info = "barrier")

  # The comparison with INLAspacetime's barrier precision is in
  # tests/optional/, since INLAspacetime is not in Suggests.
})

test_that("RTMB restricted spatial regression matches frozen TMB values", {
  d <- rtmb_spatial_data()
  mesh <- make_mesh(d, c("x", "y"), n_knots = 30L, type = "kmeans")
  fit <- sdmTMB(gaussian ~ z, data = d, mesh = mesh, time = "time",
    spatiotemporal = "iid", spatial_varying = ~ 0 + z,
    control = sdmTMBcontrol(get_rsr = TRUE), do_fit = FALSE)
  expect_frozen_fit_data("fields", fit, d, info = "rsr", project = FALSE)
  fit <- sdmTMB(response ~ z, data = d, mesh = mesh, family = delta_gamma(),
    control = sdmTMBcontrol(get_rsr = TRUE), do_fit = FALSE)
  expect_frozen_fit_data("fields", fit, d, info = "rsr delta", project = FALSE)
})
