# Fixes to model code in both backends, checked against direct calculations.

# Reports at exactly `parameters`, with no inner optimization.
report_both_backends <- function(fit, parameters) {
  lapply(c(tmb = "tmb", rtmb = "rtmb"), function(backend) {
    obj <- make_sdmTMB_adfun(fit$tmb_data, parameters, fit$tmb_map,
      random = NULL, backend = backend)
    obj$report()
  })
}

test_that("each time-varying term enters the linear predictor once", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(x = rnorm(60), year = rep(1:3, each = 20))
  d$y <- rnorm(60, 1 + 0.5 * d$x)
  fit <- sdmTMB(y ~ 0, data = d, time = "year", spatial = "off",
    spatiotemporal = "off", time_varying = ~ 1 + x, do_fit = FALSE)
  p <- fit$tmb_params
  p$b_rw_t[] <- seq(-1, 1, length.out = length(p$b_rw_t))
  expected <- p$b_rw_t[d$year, 1, 1] + d$x * p$b_rw_t[d$year, 2, 1]
  for (r in report_both_backends(fit, p)) {
    expect_equal(r$eta_rw_i[, 1], expected)
    expect_equal(r$eta_i[, 1], expected)
  }

  nd <- d[c(1, 25, 60), ]
  projection <- predict(fit, newdata = nd, return_tmb_data = TRUE)
  fit$tmb_data <- projection
  for (r in report_both_backends(fit, p)) {
    expect_equal(r$proj_eta[, 1], expected[c(1, 25, 60)])
  }
})

test_that("sigma_V priors apply once per component in delta models", {
  skip_on_cran()
  set.seed(2)
  d <- data.frame(x = rnorm(80), year = rep(1:4, each = 20))
  d$y <- ifelse(runif(80) < 0.6, rgamma(80, 2, 2 / exp(1 + 0.3 * d$x)), 0)
  fit_args <- list(y ~ 1, data = d, time = "year", spatial = "off",
    spatiotemporal = "off", time_varying = ~ 0 + x, family = delta_gamma(),
    do_fit = FALSE)
  prior <- gamma_cv(0.2, 0.5)
  fit0 <- do.call(sdmTMB, fit_args)
  fit1 <- do.call(sdmTMB, c(fit_args,
    list(priors = sdmTMBpriors(sigma_V = prior))))
  p <- fit0$tmb_params
  p$ln_tau_V[] <- c(-1, -0.5)
  sigma_V <- exp(p$ln_tau_V[1, ])
  expected <- -sum(stats::dgamma(sigma_V, shape = prior[[1]],
    scale = prior[[2]], log = TRUE))
  for (backend in c("tmb", "rtmb")) {
    nll <- vapply(list(fit0, fit1), function(fit) {
      obj <- make_sdmTMB_adfun(fit$tmb_data, p, fit$tmb_map, NULL,
        backend = backend)
      obj$fn(obj$par)
    }, numeric(1))
    expect_equal(nll[[2]] - nll[[1]], expected, info = backend)
  }
})

test_that("binomial deviance residuals allow more than one trial", {
  skip_on_cran()
  set.seed(3)
  n <- sample(1:8, 40, replace = TRUE)
  d <- data.frame(x = rnorm(40), n = n)
  d$successes <- rbinom(40, n, stats::plogis(0.3 + d$x))
  d$y <- d$successes / d$n
  fit <- sdmTMB(y ~ x, data = d, weights = d$n, spatial = "off",
    family = binomial(), do_fit = FALSE)
  p <- fit$tmb_params
  p$b_j[] <- c(0.2, 0.8)
  prob <- stats::plogis(p$b_j[[1]] + p$b_j[[2]] * d$x)
  expected <- sign(d$y - prob) *
    sqrt(stats::binomial()$dev.resids(d$y, prob, d$n))
  for (r in report_both_backends(fit, p)) {
    expect_equal(r$devresid[, 1], expected)
  }
})
