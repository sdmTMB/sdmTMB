# Simulated length-at-age with ages read for up to `n_per_bin` fish per unit
# and 2 cm length bin. Lengths are lognormal, gengamma with shape `Q`, or Gamma.
lsas_sim <- function(n_per_bin, seed = 1, Q = NULL, gamma = FALSE) {
  set.seed(seed)
  ages <- 1:12
  mu_a <- 60 * (1 - exp(-0.25 * (ages + 0.5)))
  pop <- do.call(rbind, lapply(1:3, function(u) {
    age <- sample(ages, 3000, replace = TRUE, prob = exp(-0.35 * ages))
    data.frame(unit = u, age = age,
      length = exp(rnorm(length(age), log(mu_a[age]) - 0.1^2 / 2, 0.1)))
  }))
  if (!is.null(Q)) {
    pop$length <- rtmb_obs_families$gengamma$simulate(mu_a[pop$age],
      list(phi = 0.1, Q = Q))
  }
  if (gamma) pop$length <- rgamma(nrow(pop), shape = 100, scale = mu_a[pop$age] / 100)
  breaks <- seq(0, 80, by = 2)
  pop$bin <- .lsas_bin(pop$length, breaks)
  pop$aged <- FALSE
  for (g in split(seq_len(nrow(pop)), list(pop$unit, pop$bin), drop = TRUE)) {
    pop$aged[g[sample.int(length(g), min(n_per_bin, length(g)))]] <- TRUE
  }
  counts <- stats::aggregate(cbind(n_measured = 1, n_aged = aged) ~ unit + bin,
    data = pop, FUN = sum)
  counts$length <- breaks[counts$bin]
  list(pop = pop, aged = pop[pop$aged, ], mu_a = mu_a,
    design = lsas("unit", breaks, counts))
}

lsas_fit <- function(data, family = lognormal(),
  control = sdmTMBcontrol(backend = "rtmb"), ...) {
  sdmTMB(length ~ 0 + factor(age), data = data, family = family,
    spatial = "off", control = control, ...)
}

test_that("length-stratified age sampling corrects length-at-age", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = 5)
  naive <- exp(tidy(lsas_fit(sim$aged))$estimate)
  fit <- lsas_fit(sim$aged, length_stratified = sim$design)
  corrected <- exp(tidy(fit)$estimate)
  old <- 6:12
  expect_gt(mean(naive[old] - sim$mu_a[old]), 3)
  expect_lt(abs(mean(corrected[old] - sim$mu_a[old])), 1)
})

test_that("ageing every fish leaves the likelihood unchanged when all bins are occupied", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = Inf)
  # Use occupied bins: Candy excludes empty bins even in a census.
  sim$design <- lsas("unit", c(0, 40, 80), sim$design$counts)
  fit <- lsas_fit(sim$aged, length_stratified = sim$design)
  expect_equal(logLik(fit), logLik(lsas_fit(sim$aged)))
})

test_that("the selection probability sums bin probabilities times fractions", {
  pi <- c(1, 0.5, 0.2)
  design <- list(pi = matrix(pi, 1), cuts = c(10, 20), unit_i = c(1L, 1L))
  s <- list(mu = c(12, 18), phi = 4)
  expected <- vapply(s$mu, function(mu) {
    sum(pi * diff(c(0, stats::pnorm(design$cuts, mu, s$phi), 1)))
  }, numeric(1))
  expect_equal(rtmb_lsas_log_selection(s, "gaussian", design, 1:2), log(expected))
})

test_that("the length CDFs match their densities", {
  s <- list(mu = 50, phi = 0.15)
  densities <- list(
    gaussian = function(x) dnorm(x, 50, 8),
    lognormal = function(x) dlnorm(x, log(50) - 0.15^2 / 2, 0.15)
  )
  expect_equal(rtmb_lsas_cdf(45, list(mu = 50, phi = 8), "gaussian"),
    integrate(densities$gaussian, -Inf, 45)$value, tolerance = 1e-6)
  expect_equal(rtmb_lsas_cdf(45, s, "lognormal"),
    integrate(densities$lognormal, 0, 45)$value, tolerance = 1e-6)
  expect_equal(rtmb_lsas_cdf(45, list(mu = 50, phi = 100), "Gamma"),
    integrate(function(x) dgamma(x, 100, scale = 0.5), 0, 45)$value,
    tolerance = 1e-6)
  for (Q in c(-0.8, 0.3, 1.2)) {
    s$Q <- Q
    expect_equal(rtmb_lsas_cdf(45, s, "gengamma"),
      integrate(function(x) exp(rtmb_dgengamma(x, 50, 0.15, Q)), 0, 45)$value,
      tolerance = 1e-6)
  }
})

test_that("gengamma() corrects length-at-age with shape away from zero", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = 5, Q = -0.6)
  expect_message(
    # The provisional pgamma path is unstable near Q = 0. Fix the shape to
    # test the sampling correction separately from estimating that parameter.
    fit <- lsas_fit(sim$aged, family = gengamma(), length_stratified = sim$design,
      control = sdmTMBcontrol(backend = "rtmb", start = list(gengamma_Q = -0.6),
        map = list(gengamma_Q = factor(NA)))),
    "provisional")
  corrected <- exp(tidy(fit)$estimate)
  old <- 6:12
  expect_lt(abs(mean(corrected[old] - sim$mu_a[old])), 1)
})

test_that("Gamma() also corrects length-at-age", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = 5, gamma = TRUE)
  # From the default start, the provisional CDF-difference path can step to
  # where every occupied bin's probability underflows (NaN gradient). Start
  # from the naive fit, a reasonable start in practice.
  naive <- lsas_fit(sim$aged, family = Gamma(link = "log"))$model$par
  start <- list(b_j = naive[names(naive) == "b_j"], ln_phi = naive[["ln_phi"]])
  expect_message(
    fit <- lsas_fit(sim$aged, family = Gamma(link = "log"),
      length_stratified = sim$design,
      control = sdmTMBcontrol(backend = "rtmb", start = start)), "provisional")
  corrected <- exp(tidy(fit)$estimate)
  old <- 6:12
  expect_lt(abs(mean(corrected[old] - sim$mu_a[old])), 1)
})

test_that("length_stratified checks its inputs", {
  sim <- lsas_sim(n_per_bin = 5)
  aged <- sim$aged[1:200, ]
  expect_error(sdmTMB(length ~ 1, data = aged, family = lognormal(),
    spatial = "off", control = sdmTMBcontrol(backend = "tmb"),
    length_stratified = sim$design), "rtmb")
  expect_error(sdmTMB(length ~ 1, data = aged, family = student(df = 5),
    spatial = "off", control = sdmTMBcontrol(backend = "rtmb"),
    length_stratified = sim$design), "gaussian")
  counts <- sim$design$counts
  counts$n_aged <- 0
  expect_error(lsas_fit(aged,
    length_stratified = lsas("unit", sim$design$breaks, counts)), "no aged fish")
  expect_error(lsas("unit", c(2, 1), counts), "increasing")
})

test_that("missing and invalid sampling counts cannot bypass correction", {
  counts <- data.frame(unit = "a", length = 12, n_measured = 100, n_aged = 5)
  data <- data.frame(unit = "a", length = 12)
  prepare <- function(counts, data, breaks = c(0, 20, 40)) {
    .lsas_tmb_data(lsas("unit", breaks, counts), data, data$length,
      list(family = list(family = "gaussian")), "rtmb")
  }
  expect_error(prepare(counts, transform(data, unit = "b")), "must occur")
  expect_error(prepare(counts, transform(data, unit = NA_character_)), "missing")
  expect_error(prepare(transform(counts, unit = NA_character_), data), "missing")
  expect_error(prepare(counts, transform(data, length = 22)), "no aged fish")
  expect_error(prepare(transform(counts, n_measured = 0, n_aged = 0), data),
    "no aged fish")
  for (bad in c(-1, NA, Inf, 1.5)) {
    expect_error(prepare(transform(counts, n_aged = bad), data), "nonnegative integers")
    expect_error(prepare(transform(counts, n_measured = bad), data), "nonnegative integers")
  }
  expect_error(prepare(transform(counts, length = NA_real_), data), "finite numbers")
  expect_error(prepare(counts, data, c(0, NA, 40)), "increasing")
})

test_that("Candy excludes empty strata and aggregates counts before division", {
  counts <- data.frame(unit = c("a", "a", "b"), length = c(12, 14, 25),
    n_measured = c(10, 90, 20), n_aged = c(2, 3, 4))
  data <- data.frame(unit = c("b", "a"), length = c(25, 12))
  design <- .lsas_tmb_data(lsas("unit", c(0, 20, 40), counts), data, data$length,
    list(family = list(family = "gaussian")), "rtmb")
  expect_equal(design$pi, matrix(c(0, 0.05, 0.2, 0), 2))
  s <- list(mu = c(25, 12), phi = 4)
  expected <- c(log(0.2) + pnorm(20, 25, 4, lower.tail = FALSE, log.p = TRUE),
    log(0.05) + pnorm(20, 12, 4, log.p = TRUE))
  expect_equal(rtmb_lsas_log_selection(s, "gaussian", design, 1:2), expected)
})

test_that("tail bin probabilities and their derivatives remain finite", {
  for (family in c("gaussian", "lognormal")) {
    for (side in c(-1, 1)) {
      # A finite interval deep in either tail, plus open-ended tail coverage.
      bounds <- sort(side * c(10, 11))
      q <- if (family == "lognormal") exp(bounds) else bounds
      state <- function(mu) list(mu = if (family == "lognormal") exp(mu + 0.5) else mu,
        phi = 1)
      objective <- function(p) rtmb_lsas_log_interval(q[1], q[2], state(p$mu), family)
      obj <- RTMB::MakeADFun(objective, list(mu = 0), silent = TRUE)
      expected <- if (side > 0) {
        pnorm(10, lower.tail = FALSE, log.p = TRUE) +
          log1p(-exp(pnorm(11, lower.tail = FALSE, log.p = TRUE) -
            pnorm(10, lower.tail = FALSE, log.p = TRUE)))
      } else {
        pnorm(-10, log.p = TRUE) + log1p(-exp(pnorm(-11, log.p = TRUE) -
          pnorm(-10, log.p = TRUE)))
      }
      expect_equal(obj$fn(0), expected, tolerance = 1e-10)
      # Evaluate the same tape across both tail branches.
      for (mu in c(-20, 0, 20)) {
        gradient <- (obj$fn(mu + 1e-5) - obj$fn(mu - 1e-5)) / 2e-5
        expect_equal(as.numeric(obj$gr(mu)), gradient, tolerance = 1e-6)
        expect_true(all(is.finite(obj$he(mu))))
      }
    }
  }
  design <- list(pi = matrix(c(0, 1), 1), cuts = 10, unit_i = 1L)
  obj <- RTMB::MakeADFun(function(p) rtmb_lsas_log_selection(
    list(mu = p$mu, phi = 1), "gaussian", design, 1L), list(mu = 0), silent = TRUE)
  expect_equal(obj$fn(0), pnorm(10, lower.tail = FALSE, log.p = TRUE))
  expect_true(all(is.finite(obj$gr(0))))
})

test_that("unsupported LSAS diagnostics fail explicitly", {
  sim <- lsas_sim(n_per_bin = 5)
  fit <- lsas_fit(sim$aged, length_stratified = sim$design)
  for (type in c("mle-mvn", "mle-eb", "mle-mcmc", "response", "pearson", "deviance")) {
    expect_error(residuals(fit, type = type), "not yet supported.*length_stratified")
  }
  expect_error(simulate(fit), "not yet supported.*length_stratified")
  expect_error(dharma_residuals(NULL, fit), "not yet supported.*length_stratified")
  expect_s3_class(predict(fit, newdata = data.frame(age = 1:12)), "data.frame")
})
