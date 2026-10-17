test_that("pc_anisotropy() parses its arguments", {
  expect_equal(pc_anisotropy(3), c(3, 0.05), ignore_attr = TRUE)
  expect_error(pc_anisotropy(0.5))
  expect_error(pc_anisotropy(3, 1.5))
  expect_error(pc_anisotropy(Inf))
  expect_error(pc_anisotropy(3, NA))
})

test_that("sdmTMBpriors() keeps its positional arguments", {
  p <- sdmTMBpriors(pc_matern(10, 5), pc_matern(20, 5), halfnormal(0, 1))
  expect_equal(p$phi, halfnormal(0, 1))
  expect_true(is.na(p$aniso[[1L]]))
})

test_that("Anisotropy PC prior is a normalized density of ln_H_input", {
  aniso <- rtmb_aniso_prior_inputs(list(ratio_lt = 3, ratio_prob = 0.05),
    list(matern_s = pc_matern(2, 1)))
  g <- seq(-8, 8, length.out = 401)
  h <- diff(g)[[1L]]
  grid <- expand.grid(x1 = g, x2 = g)
  sinh2 <- rtmb_aniso_sinh2(list(grid$x1, grid$x2))
  dens <- exp(rtmb_pc_aniso(sinh2, aniso$lambda_v))
  expect_equal(sum(dens) * h^2, 1, tolerance = 1e-3)
  # ratio of major to minor axis ranges is exp(|v|)
  expect_equal(sum(dens[exp(asinh(sqrt(sinh2))) > 3]) * h^2, 0.05,
    tolerance = 0.01)

  # marginal P(range < range_gt) = range_prob
  set.seed(1)
  f_v <- pc_aniso_f(0) - log(stats::runif(1e5)) / aniso$lambda_v
  range <- aniso$range_s * f_v / stats::rexp(1e5)
  expect_equal(mean(range < 2), 0.05, tolerance = 0.05)
})

test_that("Anisotropy PC prior calibrates both tail probabilities", {
  f0 <- sqrt(4 * pi / 3)
  # f(r) - f(0) from the paper's f(r), or its leading term pi * r^2 / f(0)
  # where the direct difference cancels
  df <- function(r) {
    if (r < 1e-4) return(pi * r^2 / f0)
    sqrt(pi / 3 * (3 * cosh(2 * r) + 1)) - f0
  }
  grid <- expand.grid(ratio_lt = c(10, 3, 1.5, 1 + 1e-5, 1 + 1e-7, 1 + 1e-8),
    ratio_prob = c(1e-8, 0.05, 0.999), range_prob = c(1e-8, 0.05, 0.999))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    aniso <- rtmb_aniso_prior_inputs(list(ratio_lt = g$ratio_lt,
      ratio_prob = g$ratio_prob), list(matern_s = pc_matern(2, 1,
      range_prob = g$range_prob), matern_st = pc_matern(5, 1,
      range_prob = g$range_prob)))
    expect_true(all(is.finite(unlist(aniso[c("lambda_v", "range_s",
      "range_st")]))))
    # P(ratio > ratio_lt) = exp(-lambda_v * (f(log(ratio_lt)) - f(0)))
    expect_equal(-aniso$lambda_v * df(log(g$ratio_lt)), log(g$ratio_prob),
      tolerance = 1e-8)
    # marginal P(range < range_gt) = exp(-q * f0) / (1 + q / lambda_v)
    for (x in list(c(aniso$range_s, 2), c(aniso$range_st, 5))) {
      q <- x[[1L]] / x[[2L]]
      expect_gt(q, 0)
      expect_equal(-q * f0 - log1p(q / aniso$lambda_v), log(g$range_prob),
        tolerance = 1e-10)
    }
  }
})

test_that("Anisotropy PC prior reduces to the isotropic range prior", {
  rate <- function(ratio_lt) {
    rtmb_aniso_prior_inputs(list(ratio_lt = ratio_lt, ratio_prob = 0.05),
      list(matern_s = pc_matern(2, 1)))$range_s * pc_aniso_f(0)
  }
  expect_equal(rate(1.01), -log(0.05) * 2, tolerance = 1e-3)
  expect_equal(rate(1 + 1e-8), -log(0.05) * 2, tolerance = 1e-6)
})

test_that("Anisotropy PC prior matches the paper's density of v", {
  lambda_v <- 1.6
  f <- function(r) sqrt(pi / 3 * (3 * cosh(2 * r) + 1))
  # Eq. 16 of Llamazares-Elias et al. (2024), marginal in v
  dens_v <- function(v) {
    r <- sqrt(sum(v^2))
    df <- sqrt(pi) * sinh(2 * r) / sqrt(cosh(2 * r) + 1 / 3)
    lambda_v * df / (2 * pi * r) * exp(-lambda_v * (f(r) - f(0)))
  }
  # v from H = MakeH(x), via H = cosh(r) I + sinh(r) / r * [v1 v2; v2 -v1]
  x_to_v <- function(x) {
    w <- c((exp(x[[1L]]) - (1 + x[[2L]]^2) * exp(-x[[1L]])) / 2, x[[2L]])
    asinh(sqrt(sum(w^2))) * w / sqrt(sum(w^2))
  }
  H_v <- function(v) {
    r <- sqrt(sum(v^2))
    cosh(r) * diag(2) + sinh(r) / r * matrix(c(v[1], v[2], v[2], -v[1]), 2)
  }
  for (x in list(c(0.3, -0.2), c(-1.5, 0.8), c(2, 3), c(0.01, 0))) {
    v <- x_to_v(x)
    H <- matrix(c(exp(x[1]), x[2], x[2], (1 + x[2]^2) * exp(-x[1])), 2)
    expect_equal(H_v(v), H)
    # Jacobian of x -> v by central differences
    h <- 1e-6
    J <- vapply(1:2, function(j) {
      e <- replace(numeric(2), j, h)
      (x_to_v(x + e) - x_to_v(x - e)) / (2 * h)
    }, numeric(2))
    expect_equal(rtmb_pc_aniso(rtmb_aniso_sinh2(x), lambda_v),
      log(dens_v(v)) + log(abs(det(J))), tolerance = 1e-7)
  }
})

test_that("Anisotropy PC prior has finite AD derivatives near isotropy", {
  f0 <- pc_aniso_f(0)
  for (ratio_lt in c(3, 1 + 1e-7)) {
    lambda_v <- rtmb_aniso_prior_inputs(list(ratio_lt = ratio_lt,
      ratio_prob = 0.05), list())$lambda_v
    obj <- RTMB::MakeADFun(function(p) {
      -rtmb_pc_aniso(rtmb_aniso_sinh2(p$x), lambda_v)
    }, list(x = c(0, 0)), silent = TRUE)
    # at isotropy, sinh2 ~ |x|^2, so the Hessian is
    # (2 * pi * lambda_v / f0 + 1 / 2) I
    expect_equal(obj$gr(c(0, 0)), matrix(0, 1, 2))
    expect_equal(obj$he(c(0, 0)), (2 * pi * lambda_v / f0 + 0.5) * diag(2))
    scale <- 1 / sqrt(lambda_v)
    for (x in list(c(1, -2) * scale, c(-0.5, 0.1) * scale, c(0.3, 0.2))) {
      expect_true(is.finite(obj$fn(x)))
      expect_true(all(is.finite(obj$gr(x))))
      expect_true(all(is.finite(obj$he(x))))
      # the gradient matches central differences of the objective
      h <- 1e-4 * scale
      num <- vapply(1:2, function(j) {
        e <- replace(numeric(2), j, h)
        (obj$fn(x + e) - obj$fn(x - e)) / (2 * h)
      }, numeric(1))
      expect_equal(as.vector(obj$gr(x)), num, tolerance = 1e-5)
    }
  }
})

test_that("Anisotropy PC prior fits and regularizes the anisotropy", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = tweedie(),
      anisotropy = TRUE, ...)
  }
  ratio <- function(m) {
    e <- eigen(m$tmb_obj$report()$H)$values
    sqrt(e[[1L]] / e[[2L]])
  }
  m0 <- fit()
  m1 <- fit(priors = sdmTMBpriors(matern_s = pc_matern(10, 5),
    aniso = pc_anisotropy(3)))
  expect_equal(m1$model$convergence, 0L)
  expect_lt(max(abs(m1$gradients)), 1e-4)
  expect_lt(ratio(m1), ratio(m0))
  expect_lt(ratio(m1), 3)
  expect_true(sanity(m1, silent = TRUE)$all_ok)

  mb <- fit(priors = sdmTMBpriors(matern_s = pc_matern(10, 5),
    aniso = pc_anisotropy(3)), bayesian = TRUE)
  expect_equal(mb$model$convergence, 0L)
})

test_that("Anisotropy PC prior counts a shared H once", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = delta_gamma(),
      anisotropy = TRUE, priors = sdmTMBpriors(aniso = pc_anisotropy(3)),
      do_fit = FALSE, ...)
  }
  expect_equal(fit()$tmb_data$priors_aniso$components, 1L)
  separate <- fit(control = sdmTMBcontrol(map = list(ln_H_input = factor(1:4))))
  expect_equal(separate$tmb_data$priors_aniso$components, 1:2)
})

test_that("Anisotropy PC prior checks its requirements", {
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 30)
  fit <- function(...) {
    sdmTMB(density ~ 1, data = pcod, mesh = mesh, family = tweedie(),
      priors = sdmTMBpriors(aniso = pc_anisotropy(3)), do_fit = FALSE, ...)
  }
  expect_error(fit(), "anisotropy = TRUE")
  expect_error(fit(anisotropy = TRUE,
    control = sdmTMBcontrol(backend = "tmb")), "rtmb")
})
