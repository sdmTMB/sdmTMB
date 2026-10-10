rtmb_ctl <- function(...) sdmTMBcontrol(backend = "rtmb", ...)

custom_prior_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), f = factor(rep(c("a", "b", "c"), n / 3)),
    g = factor(rep(seq_len(6), each = n / 6)))
  u <- stats::rnorm(6, 0, 0.5)
  d$y <- 1 + 0.5 * d$x + c(0, 0.3, -0.2)[d$f] + u[d$g] + stats::rnorm(n, 0, 0.7)
  d
}

test_that("custom prior arguments are validated", {
  expect_error(sdmTMBpriors(custom = 1), "function")
  expect_error(sdmTMBpriors(custom_log_jacobian = function(par, theta) 0),
    "requires `custom`")
  expect_null(sdmTMBpriors()$custom)
  d <- custom_prior_data()
  expect_error(
    sdmTMB(y ~ x, data = d, spatial = "off", do_fit = FALSE,
      control = sdmTMBcontrol(backend = "tmb"),
      priors = sdmTMBpriors(custom = function(par, theta) 0)),
    "rtmb"
  )
  fit <- function(custom) {
    sdmTMB(y ~ x, data = d, spatial = "off", do_fit = FALSE,
      control = rtmb_ctl(), priors = sdmTMBpriors(custom = custom))
  }
  expect_error(fit(function(par, theta) stop("oops")), "`custom` prior function failed")
  expect_error(fit(function(par, theta) "a"), "numeric log densities")
  expect_error(fit(function(par, theta) numeric(0)), "numeric log densities")
  expect_error(fit(function(par, theta) log(0)), "non-finite")
  # arrays are summed like vectors
  expect_s3_class(fit(function(par, theta) matrix(0, 2, 2)), "sdmTMB")
  expect_error(
    sdmTMB(y ~ x, data = d, spatial = "off", do_fit = FALSE, control = rtmb_ctl(),
      bayesian = TRUE, priors = sdmTMBpriors(custom = function(par, theta) 0,
        custom_log_jacobian = function(par, theta) stop("oops"))),
    "`custom_log_jacobian` prior function failed"
  )
  # an unused Jacobian isn't evaluated
  expect_s3_class(
    sdmTMB(y ~ x, data = d, spatial = "off", do_fit = FALSE, control = rtmb_ctl(),
      priors = sdmTMBpriors(custom = function(par, theta) 0,
        custom_log_jacobian = function(par, theta) stop("oops"))),
    "sdmTMB"
  )
  # stats::dnorm() can't take AD values
  expect_error(
    fit(function(par, theta) stats::dnorm(par$b_j[1], log = TRUE)),
    "`custom` prior function failed"
  )
})

test_that("NULL custom priors leave the objective unchanged", {
  d <- custom_prior_data()
  m <- sdmTMB(y ~ x, data = d, spatial = "off", control = rtmb_ctl(),
    do_fit = FALSE)
  expect_null(m$tmb_data$priors_custom)
  expect_false("custom" %in% names(rtmb_prepare(m$tmb_data)$priors))
})

test_that("a custom normal prior on b matches the built-in prior", {
  d <- custom_prior_data()
  fit <- function(priors, do_fit = TRUE) {
    sdmTMB(y ~ x + f, data = d, spatial = "off", control = rtmb_ctl(),
      priors = priors, do_fit = do_fit)
  }
  builtin <- sdmTMBpriors(b = normal(c(0, 0, 0, 0), c(1, 2, 0.5, 0.5)))
  custom <- sdmTMBpriors(custom = function(par, theta) {
    RTMB::dnorm(par$b_j, 0, c(1, 2, 0.5, 0.5), log = TRUE)
  })
  b1 <- fit(builtin, FALSE)$tmb_obj
  c1 <- fit(custom, FALSE)$tmb_obj
  p <- rtmb_perturb(b1$par)
  expect_equal(c1$fn(p), b1$fn(p))
  expect_equal(c1$gr(p), b1$gr(p))
  expect_equal(c1$he(p), b1$he(p))

  b2 <- fit(builtin)
  c2 <- fit(custom)
  expect_equal(c2$model$par, b2$model$par, tolerance = 1e-6)
  expect_equal(c2$model$objective, b2$model$objective, tolerance = 1e-8)
  expect_sdreport_matches(c2, b2)

  terms <- get_prior_densities(c2)
  expect_equal(nrow(terms), 4L)
  expect_equal(terms$type, rep("density", 4))
  expect_equal(sum(terms$log_density),
    sum(stats::dnorm(c2$model$par[1:4], 0, c(1, 2, 0.5, 0.5), log = TRUE)))
})

test_that("custom priors select delta model components", {
  d <- custom_prior_data()
  d$y <- ifelse(d$y > 1, exp(d$y), 0)
  fit <- function(priors) {
    sdmTMB(y ~ x, data = d, spatial = "off", control = rtmb_ctl(),
      family = delta_gamma(), priors = priors, do_fit = FALSE)$tmb_obj
  }
  b <- fit(sdmTMBpriors(b = normal(c(0, 0), c(1, 2))))
  c <- fit(sdmTMBpriors(custom = function(par, theta) {
    c(RTMB::dnorm(par$b_j, 0, c(1, 2), log = TRUE),
      RTMB::dnorm(par$b_j2, 0, c(1, 2), log = TRUE))
  }))
  p <- rtmb_perturb(b$par)
  expect_equal(c$fn(p), b$fn(p))
  expect_equal(c$gr(p), b$gr(p))
})

test_that("transformed custom priors add a Jacobian only if bayesian", {
  d <- custom_prior_data()
  fit <- function(priors = sdmTMBpriors(), bayesian = FALSE) {
    sdmTMB(y ~ x, data = d, spatial = "off", control = rtmb_ctl(),
      priors = priors, bayesian = bayesian, do_fit = FALSE)$tmb_obj
  }
  priors <- sdmTMBpriors(
    custom = function(par, theta) RTMB::dnorm(theta$phi[1], 0, 5, log = TRUE),
    custom_log_jacobian = function(par, theta) par$ln_phi[1]
  )
  base <- fit()
  p <- rtmb_perturb(base$par)
  ln_phi <- p[["ln_phi"]]
  prior <- stats::dnorm(exp(ln_phi), 0, 5, log = TRUE)
  expect_equal(fit(priors)$fn(p), base$fn(p) - prior)
  expect_equal(fit(priors, bayesian = TRUE)$fn(p),
    fit(bayesian = TRUE)$fn(p) - prior - ln_phi)
  # derivative of -(log density + Jacobian) with respect to ln_phi
  phi <- exp(ln_phi)
  expect_equal(
    fit(priors, bayesian = TRUE)$gr(p)[, names(p) == "ln_phi"] -
      fit(bayesian = TRUE)$gr(p)[, names(p) == "ln_phi"],
    phi^2 / 25 - 1
  )
})

test_that("custom priors can penalize contrasts and use names", {
  d <- custom_prior_data()
  fit <- function(priors = sdmTMBpriors()) {
    sdmTMB(y ~ x + f, data = d, spatial = "off", control = rtmb_ctl(),
      priors = priors, do_fit = FALSE)$tmb_obj
  }
  m <- fit(sdmTMBpriors(custom = function(par, theta) {
    c(contrast = RTMB::dnorm(par$b_j[3] - par$b_j[4], 0, 0.1, log = TRUE),
      slope = RTMB::dnorm(par$b_j[2], 0.5, 1, log = TRUE))
  }))
  base <- fit()
  p <- rtmb_perturb(base$par)
  b <- unname(p[names(p) == "b_j"])
  expect_equal(m$fn(p), base$fn(p) -
    stats::dnorm(b[3] - b[4], 0, 0.1, log = TRUE) -
    stats::dnorm(b[2], 0.5, 1, log = TRUE))
})

test_that("custom priors on random effects enter the Laplace approximation", {
  d <- custom_prior_data()
  # An extra N(0, 1) density on each random intercept u ~ N(0, sigma^2) gives
  # u ~ N(0, s^2), s^2 = sigma^2 / (sigma^2 + 1), times N(0; 0, sigma^2 + 1).
  # The Gaussian model's Laplace approximation is exact.
  m <- sdmTMB(y ~ x + (1 | g), data = d, spatial = "off", control = rtmb_ctl(),
    do_fit = FALSE, priors = sdmTMBpriors(custom = function(par, theta) {
      RTMB::dnorm(par$re_b_pars[seq_len(6), 1], 0, 1, log = TRUE)
    }))
  obj <- m$tmb_obj
  p <- obj$par
  p[] <- c(0.8, 0.4, log(0.6), log(0.5))[seq_along(p)]
  expect_equal(names(p), c("b_j", "b_j", "ln_phi", "re_cov_pars"))
  b <- p[1:2]
  phi <- exp(p[[3]])
  sigma <- exp(p[[4]])
  s2 <- sigma^2 / (sigma^2 + 1)
  X <- cbind(1, d$x)
  Z <- stats::model.matrix(~ 0 + g, d)
  V <- s2 * tcrossprod(Z) + diag(phi^2, nrow(d))
  r <- d$y - X %*% b
  R <- chol(V)
  loglik <- -0.5 * sum(backsolve(R, r, transpose = TRUE)^2) -
    sum(log(diag(R))) - 0.5 * nrow(d) * log(2 * pi)
  expected <- -(loglik + 6 * stats::dnorm(0, 0, sqrt(sigma^2 + 1), log = TRUE))
  expect_equal(obj$fn(p), expected, tolerance = 1e-8, ignore_attr = TRUE)

  expect_false(any(get_prior_parameters(m)$status == "random"))
  re <- get_prior_parameters(m, random = TRUE)
  expect_equal(re$expression[re$status == "random"],
    paste0("par$re_b_pars[", 1:6, ", 1]"))
})

test_that("custom prior parameters are listed with labels and map status", {
  d <- custom_prior_data()
  m <- sdmTMB(y ~ x + f, data = d, spatial = "off", control = rtmb_ctl(
    map = list(b_j = factor(c(1, NA, 2, 2))),
    start = list(b_j = c(0, 0.5, 0, 0))), do_fit = FALSE)
  p <- get_prior_parameters(m)
  b <- p[p$name == "b_j", ]
  # fixed elements are dropped
  expect_equal(b$expression, paste0("par$b_j[", c(1, 3, 4), "]"))
  expect_equal(b$label, c("(Intercept)", "fb", "fc"))
  expect_equal(b$status, c("estimated", "shared", "shared"))
  phi <- p[p$expression == "theta$phi[1]", ]
  expect_equal(phi$status, "derived")
  expect_equal(phi$value, exp(p$value[p$expression == "par$ln_phi[1]"]))
  # no spatial field, so sigma_O is zero
  expect_equal(p$value[p$expression == "theta$sigma_O[1, 1]"], 0)
  expect_false(any(p$status == "fixed"))
  expect_error(get_prior_densities(m), "no custom priors")
  expect_error(get_prior_parameters(
    sdmTMB(y ~ x, data = d, spatial = "off", do_fit = FALSE,
      control = sdmTMBcontrol(backend = "tmb"))), "rtmb")
})

test_that("custom priors persist through post-fit methods", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  fit <- function(priors = sdmTMBpriors(), ...) {
    sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
      family = tweedie(), control = rtmb_ctl(...), priors = priors)
  }
  # a strong prior on the spatial SD changes the fit
  priors <- sdmTMBpriors(custom = function(par, theta) {
    RTMB::dnorm(theta$log_sigma_O[1, 1], log(0.5), 0.05, log = TRUE)
  })
  m0 <- fit()
  m <- fit(priors)
  expect_equal(tidy(m, "ran_pars")$estimate[tidy(m, "ran_pars")$term == "sigma_O"],
    0.5, tolerance = 0.1)
  expect_false(isTRUE(all.equal(m$model$objective, m0$model$objective)))

  nd <- replicate_df(qcs_grid, "year", 2011)
  nd <- nd[seq(1, nrow(nd), by = 20), ]
  ind <- get_index(m, newdata = nd, bias_correct = FALSE)
  ind0 <- get_index(m0, newdata = nd, bias_correct = FALSE)
  expect_false(isTRUE(all.equal(ind$se, ind0$se)))

  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(m, f)
  m2 <- readRDS(f)
  expect_equal(get_index(m2, newdata = nd, bias_correct = FALSE), ind)
  expect_equal(predict(m2, newdata = nd, se_fit = TRUE),
    predict(m, newdata = nd, se_fit = TRUE))

  # the first multiphase phase skips custom priors, which may refer to fields
  # that are off in that phase; the final fit includes them
  m3 <- fit(priors, multiphase = TRUE)
  expect_equal(m3$model$objective, m$model$objective, tolerance = 1e-6)

  # extra optimization rebuilds the objective with the prior
  m4 <- run_extra_optimization(m, nlminb_loops = 1, newton_loops = 1)
  expect_equal(m4$model$objective, m$model$objective, tolerance = 1e-6)
})

test_that("cross-validation refits keep custom priors but scores exclude them", {
  skip_on_cran()
  d <- custom_prior_data()
  d$fold <- rep(1:2, length.out = nrow(d))
  d$X <- stats::runif(nrow(d))
  d$Y <- stats::runif(nrow(d))
  mesh <- make_mesh(d, c("X", "Y"), n_knots = 10, type = "kmeans")
  cv <- function(priors = sdmTMBpriors()) {
    sdmTMB_cv(y ~ x, data = d, mesh = mesh, spatial = "off", fold_ids = d$fold,
      control = rtmb_ctl(), priors = priors, predictive = "mle-eb",
      save_models = TRUE)
  }
  # a constant density changes the objective but not the estimates
  constant <- cv(sdmTMBpriors(custom = function(par, theta) -100))
  base <- cv()
  expect_equal(constant$models[[1]]$model$objective,
    base$models[[1]]$model$objective + 100)
  expect_equal(constant$data$cv_loglik, base$data$cv_loglik, tolerance = 1e-6)
})

test_that("theta$sigma_E matches the reported spatiotemporal SD", {
  m <- sdmTMB(density ~ 1, data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie(), time = "year", spatiotemporal = "iid",
    control = rtmb_ctl(), do_fit = FALSE)
  obj <- m$tmb_obj
  p <- rtmb_perturb(obj$env$last.par)
  theta <- rtmb_transform(obj$env$parList(par = p), rtmb_prepare(m$tmb_data))
  sigma_E <- obj$report(p)$sigma_E
  expect_equal(c(sigma_E), rep(theta$sigma_E[1, 1], nrow(sigma_E)))
})

test_that("theta$sigma_E is zero without a spatiotemporal field", {
  fit <- function(priors = sdmTMBpriors()) {
    sdmTMB(density ~ 1, data = pcod_2011, mesh = pcod_mesh_2011,
      family = tweedie(), control = rtmb_ctl(), priors = priors,
      do_fit = FALSE)$tmb_obj
  }
  base <- fit()
  m <- fit(sdmTMBpriors(custom = function(par, theta) {
    RTMB::dnorm(theta$sigma_E[1, 1], 0, 0.1, log = TRUE)
  }))
  p <- rtmb_perturb(base$par)
  # a constant density, so the spatial range is unaffected
  expect_equal(m$fn(p) - base$fn(p), -stats::dnorm(0, 0, 0.1, log = TRUE),
    ignore_attr = TRUE)
  expect_equal(m$gr(p), base$gr(p))
  # sigma_E is still reported, as with the TMB backend
  expect_true(all(base$report(base$env$last.par)$sigma_E > 0))
})

test_that("get_prior_densities() evaluates the Jacobian only if bayesian", {
  d <- custom_prior_data()
  fit <- function(bayesian, jacobian) {
    sdmTMB(y ~ x, data = d, spatial = "off", control = rtmb_ctl(),
      bayesian = bayesian, do_fit = FALSE, priors = sdmTMBpriors(
        custom = function(par, theta) RTMB::dnorm(theta$phi[1], 0, 5, log = TRUE),
        custom_log_jacobian = jacobian))
  }
  ok <- function(par, theta) par$ln_phi[1]
  # an unused Jacobian isn't evaluated, so a broken one doesn't matter
  terms <- get_prior_densities(fit(FALSE, function(par, theta) stop("oops")))
  expect_equal(terms$type, "density")
  terms <- get_prior_densities(fit(TRUE, ok))
  expect_equal(terms$type, c("density", "log_jacobian"))
  expect_equal(terms$log_density[2], 0)
})
