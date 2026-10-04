nngp_test_coords <- function(n = 30L) {
  set.seed(2)
  cbind(stats::runif(n, 0, 10), stats::runif(n, 0, 10))
}

test_that("make_nngp_data builds a reproducible predecessor graph", {
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")
  coords <- nngp_test_coords()

  set.seed(1)
  seed <- .Random.seed
  graph <- make_nngp_data(coords, 4L)
  expect_identical(.Random.seed, seed)
  expect_identical(make_nngp_data(coords, 4L), graph)

  expect_identical(sort(graph$order), seq_len(nrow(coords)))
  expect_identical(lengths(graph$neighbours), pmin(0:29, 4L))
  position <- match(seq_len(nrow(coords)), graph$order)
  for (i in seq_along(graph$order)) {
    expect_true(all(position[graph$neighbours[[i]]] < i))
  }
  expect_equal(graph$distance_to[[5]],
    sqrt(colSums((t(coords[graph$neighbours[[5]], ]) - coords[graph$order[5], ])^2)))

  expect_identical(make_nngp_data(matrix(c(1, 2), 1L), 4L)$neighbours,
    list(integer(0)))
  expect_error(make_nngp_data(coords[c(1, 1:3), ], 2L), "unique")
})

test_that("NNGP density matches the dense Gaussian it defines", {
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")
  coords <- nngp_test_coords(15L)
  n <- nrow(coords)
  kappa <- 0.6
  sigma <- 1.7
  x <- stats::rnorm(n)

  # Conditioning on all predecessors gives the (regularized) parent GP
  full <- nngp_factor(make_nngp_data(coords, n - 1L), kappa)
  R <- nngp_correlation(as.matrix(stats::dist(coords)), kappa)
  R <- (R + diag(nngp_nugget, n)) / (1 + nngp_nugget)
  expect_equal(nngp_nll(x, full, sigma),
    -mvtnorm::dmvnorm(x, sigma = sigma^2 * R, log = TRUE))

  # With few neighbours: Q = (I - B)' D^-1 (I - B) in the original order
  graph <- make_nngp_data(coords, 3L)
  f <- nngp_factor(graph, kappa)
  B <- matrix(0, n, n)
  for (i in seq_len(n)) {
    j <- graph$neighbours[[i]]
    B[graph$order[i], j] <- f$weights[i, seq_along(j)]
  }
  L <- diag(n) - B
  Q <- crossprod(L, diag(1 / (sigma * f$sd[match(seq_len(n), f$order)])^2) %*% L)
  expect_equal(nngp_nll(x, f, sigma),
    -mvtnorm::dmvnorm(x, sigma = solve(Q), log = TRUE))
  expect_equal(nngp_nll(cbind(x, 2 * x), f, sigma),
    nngp_nll(x, f, sigma) + nngp_nll(2 * x, f, sigma))

  set.seed(1)
  draws <- simulate_nngp(f, sigma, 20000L)
  expect_equal(stats::cov(t(draws)), solve(Q), tolerance = 0.05)
})

test_that("Matérn NNGP correlation matches the Matérn nu = 1 and its dense GP", {
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")
  kappa <- 0.6
  d <- c(0, 0.1, 1, 5, 20)
  expect_equal(nngp_correlation(d, kappa, "matern"),
    c(1, (kappa * d[-1]) * besselK(kappa * d[-1], 1)))
  # Correlation at the range sqrt(8) / kappa, as for the SPDE
  expect_equal(nngp_correlation(sqrt(8) / kappa, kappa, "matern"),
    sqrt(8) * besselK(sqrt(8), 1))

  coords <- nngp_test_coords(15L)
  n <- nrow(coords)
  x <- stats::rnorm(n)
  graph <- make_nngp_data(coords, n - 1L, "matern")
  expect_identical(graph$covariance, "matern")
  R <- nngp_correlation(as.matrix(stats::dist(coords)), kappa, "matern")
  R <- (R + diag(nngp_nugget, n)) / (1 + nngp_nugget)
  expect_equal(nngp_nll(x, nngp_factor(graph, kappa), 1.7),
    -mvtnorm::dmvnorm(x, sigma = 1.7^2 * R, log = TRUE))

  # The kappa gradient through besselK() matches finite differences
  f <- RTMB::MakeTape(function(k) nngp_nll(x, nngp_factor(graph, k), 1.7),
    kappa)
  h <- 1e-5
  expect_equal(as.numeric(f$jacobian(kappa)),
    (f(kappa + h) - f(kappa - h)) / (2 * h), tolerance = 1e-5)
  expect_error(sdmTMBcontrol(nngp_covariance = "gaussian"))
})

test_that("NNGP models fit, predict, and simulate with the RTMB backend", {
  skip_on_cran()
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")
  d <- subset(pcod, year >= 2013)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 25)
  fit <- sdmTMB(density ~ 1, data = d, mesh = mesh, time = "year",
    family = tweedie(), spatial_model = "nngp", spatiotemporal = "ar1",
    spatial_varying = ~ 0 + depth_scaled)
  expect_true(all(unlist(sanity(fit, silent = TRUE))))
  expect_false(is.null(fit$tmb_data$nngp_data))

  nd <- replicate_df(qcs_grid, "year", unique(d$year))
  p <- predict(fit, newdata = nd, nsim = 3L)
  expect_identical(dim(p), c(nrow(nd), 3L))
  expect_length(simulate(fit, nsim = 2L, type = "mle-mvn"), 2L * nrow(d))

  fit_matern <- sdmTMB(density ~ 1, data = d, mesh = mesh, time = "year",
    family = tweedie(), spatial_model = "nngp", spatiotemporal = "ar1",
    control = sdmTMBcontrol(nngp_covariance = "matern"))
  expect_true(all(unlist(sanity(fit_matern, silent = TRUE))))
  expect_identical(fit_matern$tmb_data$nngp_data$covariance, "matern")
  expect_length(predict(fit_matern, newdata = nd, nsim = 2L), 2L * nrow(nd))

  expect_error(
    sdmTMB(density ~ 1, data = d, mesh = mesh, spatial_model = "nngp",
      control = sdmTMBcontrol(backend = "tmb")),
    "rtmb"
  )
  expect_error(
    sdmTMB(density ~ 1, data = d, mesh = mesh, spatial_model = "nngp",
      anisotropy = TRUE),
    "not supported"
  )
})

test_that("New NNGP fields simulate with the right AR1 and SVC covariance", {
  skip_on_cran()
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")
  set.seed(3)
  n_s <- 20L
  xy <- cbind(stats::runif(n_s, 0, 10), stats::runif(n_s, 0, 10))
  d <- data.frame(x = xy[, 1], y = xy[, 2], year = rep(1:3, each = n_s),
    covariate = stats::rnorm(3L * n_s), z = 0)
  mesh <- make_mesh(d, c("x", "y"), mesh = fmesher::fm_delaunay_2d(xy))
  fit <- sdmTMB(z ~ 1, data = d, mesh = mesh, time = "year",
    spatial_model = "nngp", spatiotemporal = "ar1",
    spatial_varying = ~ 0 + covariate, do_fit = FALSE)

  kappa <- 0.5
  sigma_O <- 0.5
  sigma_E <- 0.8
  sigma_Z <- 0.4
  rho <- 0.6
  ln_tau <- function(sigma) -log(sigma) - log(kappa) - log(4 * pi) / 2
  par <- fit$tmb_obj$env$last.par
  par[names(par) == "b_j"] <- 0
  par[names(par) == "ln_kappa"] <- log(kappa)
  par[names(par) == "ln_tau_O"] <- ln_tau(sigma_O)
  par[names(par) == "ln_tau_E"] <- ln_tau(sigma_E)
  par[names(par) == "ln_tau_Z"] <- ln_tau(sigma_Z)
  par[names(par) == "ar1_phi"] <- stats::qlogis((rho + 1) / 2)
  nsim <- 2000L
  sims <- simulate(fit, nsim = nsim, seed = 1, re_form = ~0,
    observation_error = FALSE, silent = TRUE,
    mcmc_samples = matrix(par, length(par), nsim))

  # NNGP covariance at the observations from Q = (I - B)' D^-1 (I - B)
  graph <- fit$tmb_data$nngp_data
  f <- nngp_factor(graph, kappa)
  B <- matrix(0, n_s, n_s)
  for (i in seq_len(n_s)) {
    j <- graph$neighbours[[i]]
    B[graph$order[i], j] <- f$weights[i, seq_along(j)]
  }
  L <- diag(n_s) - B
  D <- diag(1 / f$sd[match(seq_len(n_s), f$order)]^2)
  A <- as.matrix(fit$tmb_data$A_st)
  C <- A %*% solve(crossprod(L, D %*% L)) %*% t(A)

  # Recover the variance multiplying C at each time lag and for the SVC
  lag <- abs(outer(d$year, d$year, "-"))
  X <- cbind(lag0 = c((lag == 0) * C), lag1 = c((lag == 1) * C),
    lag2 = c((lag == 2) * C), svc = c(outer(d$covariate, d$covariate) * C))
  estimate <- stats::coef(stats::lm(c(stats::cov(t(sims))) ~ 0 + X))
  expected <- c(sigma_O^2 + sigma_E^2 * rho^(0:2), sigma_Z^2)
  expect_equal(unname(estimate), expected, tolerance = 0.05)
})

test_that("NNGP fits match tinyVAST's exponential NNGP", {
  skip_on_cran()
  skip_if_not_installed("tinyVAST")
  skip_if_not_installed("sf")
  skip_if_not_installed("GpGp")
  skip_if_not_installed("GPvecchia")

  # tinyVAST assigns observations to polygons and uses their centroids as
  # NNGP vertices. With a mesh whose vertices are those centroids and
  # observations at the vertices, the A matrix is one-hot and the two models
  # are the same.
  set.seed(1)
  points <- sf::st_as_sf(data.frame(x = stats::runif(80, 0, 10),
    y = stats::runif(80, 0, 10)), coords = c("x", "y"))
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10)))
  cells <- sf::st_intersection(sf::st_collection_extract(
    sf::st_voronoi(sf::st_union(points), box)), box)
  xy <- sf::st_coordinates(sf::st_centroid(cells))
  d <- data.frame(x = xy[, 1], y = xy[, 2])
  field <- simulate_nngp(nngp_factor(make_nngp_data(xy, 20L), 0.7), 0.8)
  d$z <- 1 + as.vector(field) + stats::rnorm(nrow(d), 0, 0.3)

  tv <- tinyVAST::tinyVAST(z ~ 1, data = d, space_term = "",
    spatial_domain = tinyVAST::make_nngp_domain(cells, nn = 6),
    space_columns = c("x", "y"),
    control = tinyVAST::tinyVASTcontrol(gmrf_parameterization = "projection"))
  mesh <- make_mesh(d, c("x", "y"), mesh = fmesher::fm_delaunay_2d(xy))
  fit <- sdmTMB(z ~ 1, data = d, mesh = mesh, spatial_model = "nngp",
    control = sdmTMBcontrol(nngp_neighbours = 6L))

  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(tv)),
    tolerance = 1e-7)
  # tinyVAST's correlation is exp(-d / e), with e-folding distance e
  # = exp(log_kappa); sdmTMB's range is 2 e
  pars <- tidy(fit, "ran_pars")
  expect_equal(pars$estimate[pars$term == "range"] / 2,
    exp(tv$opt$par[["log_kappa"]]), tolerance = 1e-4)
  expect_equal(pars$estimate[pars$term == "sigma_O"],
    abs(tv$opt$par[["theta_z"]]), tolerance = 1e-4)
})
