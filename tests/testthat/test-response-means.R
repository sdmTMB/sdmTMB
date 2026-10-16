test_that("Means without observation error are response expectations", {
  cases <- list(
    list(family = betabinomial(), formula = cbind(y, 8 - y) ~ 1, eta = 0,
      expected = 4),
    list(family = truncated_nbinom2(), eta = log(0.1), expected = 1.1),
    list(family = truncated_nbinom1(), eta = log(0.5),
      expected = 0.5 / (1 - 2^-0.5)),
    list(family = gamma_mix(), eta = log(2), expected = 2.4),
    list(family = ordbeta(), eta = 0, y = c(0, 0.3, 0.5, 1),
      expected = plogis(-2) + (1 - plogis(-1) - plogis(-2)) * 0.5)
  )
  for (case in cases) for (backend in c("tmb", "rtmb")) {
    d <- data.frame(y = if (is.null(case$y)) 1:4 else case$y)
    fit <- sdmTMB(if (is.null(case$formula)) y ~ 1 else case$formula,
      data = d, family = case$family, spatial = "off", do_fit = FALSE,
      control = sdmTMBcontrol(backend = backend, multiphase = FALSE))
    p <- fit$tmb_params
    p$b_j[] <- case$eta
    p$ln_phi[] <- 0
    p$logit_p_extreme[] <- qlogis(0.1)
    p$log_ratio_mix[] <- log(2) # ratio 3
    if (length(p$psi)) p$psi[] <- c(-1, log(3)) # cutpoints -1 and 2
    td <- fit$tmb_data
    td$sim_obs <- 0L
    obj <- make_sdmTMB_adfun(td, p, fit$tmb_map, fit$tmb_random,
      backend = backend)
    expect_equal(obj$simulate()$y_i[1], case$expected,
      label = paste(case$family$family, backend))
  }
})

test_that("Mixture predictions are means for any link, including population", {
  d <- data.frame(y = 1:4)
  for (link in c("log", "identity")) for (backend in c("tmb", "rtmb")) {
    fit <- sdmTMB(y ~ 1, data = d, family = gamma_mix(link = link),
      spatial = "off", do_fit = FALSE,
      control = sdmTMBcontrol(backend = backend, multiphase = FALSE))
    p <- fit$tmb_params
    p$b_j[] <- stats::make.link(link)$linkfun(2)
    p$logit_p_extreme[] <- qlogis(0.1)
    p$log_ratio_mix[] <- log(2)
    td <- predict(fit, newdata = d, return_tmb_data = TRUE)
    r <- make_sdmTMB_adfun(td, p, fit$tmb_map, fit$tmb_random,
      backend = backend)$report()
    inv <- stats::make.link(link)$linkinv
    expect_equal(inv(r$proj_eta[1, 1]), 2.4, label = paste(link, backend))
    expect_equal(inv(r$proj_fe[1, 1]), 2.4, label = paste(link, backend))
  }
})

test_that("gengamma() rejects parameters without a finite mean", {
  for (backend in c("tmb", "rtmb")) {
    obj <- sdmTMB(y ~ 1, data = data.frame(y = c(0.3, 1, 3)),
      family = gengamma(), spatial = "off", do_fit = FALSE,
      control = sdmTMBcontrol(backend = backend, multiphase = FALSE))$tmb_obj
    p <- obj$par
    p[names(p) == "b_j"] <- 0
    p[names(p) == "ln_phi"] <- 0
    p[names(p) == "gengamma_Q"] <- -2 # 1 + sigma * Q < 0
    expect_true(is.nan(obj$fn(p)), label = backend)
    p[names(p) == "gengamma_Q"] <- -0.5
    expect_true(is.finite(obj$fn(p)), label = backend)
  }
})

test_that("Response-scale draws use each draw's ordbeta cutpoints and NB dispersion", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(y = c(0, 0.2, 0.5, 0.7, 1, 0.4, 0.6, 0.3))
  fit_ob <- sdmTMB(y ~ 1, data = d, family = ordbeta(), spatial = "off")
  d2 <- data.frame(y = c(1, 2, 4, 7, 3, 2, 5, 9))
  fit_nb <- sdmTMB(y ~ 1, data = d2, family = truncated_nbinom2(),
    spatial = "off")
  # ordered beta: vary only the cutpoints
  vals <- list(c(-3, log(2)), c(0, log(2)))
  best <- fit_ob$tmb_obj$env$last.par.best
  m <- vapply(vals, function(v) {
    x <- best
    x[names(x) == "psi"] <- v
    x
  }, numeric(length(best)))
  p <- predict(fit_ob, newdata = d[1, , drop = FALSE], type = "response",
    mcmc_samples = m)
  eta <- fit_ob$parlist$b_j
  expected <- vapply(vals, function(v) {
    cp <- c(v[1], v[1] + exp(v[2]))
    p0 <- plogis(cp[1] - eta)
    p1 <- plogis(eta - cp[2])
    p1 + (1 - p0 - p1) * plogis(eta)
  }, numeric(1))
  expect_equal(as.numeric(p[1, ]), expected)
  expect_gt(abs(diff(expected)), 0.01)

  # truncated NB2: vary only phi
  best <- fit_nb$tmb_obj$env$last.par.best
  lphi <- log(c(0.2, 20))
  m <- vapply(lphi, function(v) {
    x <- best
    x[names(x) == "ln_phi"] <- v
    x
  }, numeric(length(best)))
  p <- predict(fit_nb, newdata = d2[1, , drop = FALSE], type = "response",
    mcmc_samples = m)
  mu <- exp(fit_nb$parlist$b_j)
  expected <- vapply(exp(lphi), function(phi) {
    mu / (1 - dnbinom(0, size = phi, mu = mu))
  }, numeric(1))
  expect_equal(as.numeric(p[1, ]), expected)
  expect_gt(abs(diff(expected)), 0.1)
})
