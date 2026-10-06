# Simulated lognormal (or gaussian) length-at-age with ages read for up to
# `n_per_bin` fish per unit and 2 cm length bin.
lsas_sim <- function(n_per_bin, seed = 1, gaussian = FALSE) {
  set.seed(seed)
  ages <- 1:12
  mu_a <- 60 * (1 - exp(-0.25 * (ages + 0.5)))
  pop <- do.call(rbind, lapply(1:3, function(u) {
    age <- sample(ages, 3000, replace = TRUE, prob = exp(-0.35 * ages))
    data.frame(unit = u, age = age,
      length = exp(rnorm(length(age), log(mu_a[age]) - 0.1^2 / 2, 0.1)))
  }))
  if (gaussian) pop$length <- rnorm(nrow(pop), mu_a[pop$age], 4)
  breaks <- seq(0, 80, by = 2)
  bin <- .lsas_bin(pop$length, breaks)
  pop$aged <- FALSE
  for (g in split(seq_len(nrow(pop)), list(pop$unit, bin), drop = TRUE)) {
    pop$aged[g[sample.int(length(g), min(n_per_bin, length(g)))]] <- TRUE
  }
  list(pop = pop, aged = pop[pop$aged, ], mu_a = mu_a,
    design = length_strata(pop, unit = "unit", breaks = breaks))
}

lsas_fit <- function(data, family = lognormal(), ...) {
  sdmTMB(length ~ 0 + factor(age), data = data, family = family,
    spatial = "off", control = sdmTMBcontrol(backend = "rtmb"), ...)
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

test_that("gaussian() also corrects length-at-age", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = 5, gaussian = TRUE)
  naive <- tidy(lsas_fit(sim$aged, family = gaussian()))$estimate
  fit <- lsas_fit(sim$aged, family = gaussian(), length_stratified = sim$design)
  old <- 6:12
  expect_gt(mean(naive[old] - sim$mu_a[old]), 2)
  expect_lt(abs(mean(tidy(fit)$estimate[old] - sim$mu_a[old])), 1)
})

test_that("ageing every fish leaves the likelihood unchanged when all bins are occupied", {
  skip_on_cran()
  sim <- lsas_sim(n_per_bin = Inf)
  # Use occupied bins: Candy excludes empty bins even in a census.
  design <- length_strata(sim$pop, "unit", c(0, 40, 80))
  fit <- lsas_fit(sim$aged, length_stratified = design)
  expect_equal(logLik(fit), logLik(lsas_fit(sim$aged)))
})

test_that("the selection probability sums bin probabilities times fractions", {
  frac <- c(1, 0.5, 0.2)
  design <- list(frac = matrix(frac, 1), breaks = c(0, 10, 20, 30), unit_i = c(1L, 1L))
  s <- list(mu = c(12, 18), phi = 4)
  expected <- vapply(s$mu, function(mu) {
    sum(frac * diff(c(0, stats::pnorm(c(10, 20), mu, s$phi), 1)))
  }, numeric(1))
  expect_equal(rtmb_lsas_log_selection(s, "gaussian", design, 1:2), log(expected))
  s <- list(mu = c(12, 18), phi = 0.3)
  meanlog <- log(s$mu) - s$phi^2 / 2
  expected <- vapply(meanlog, function(m) {
    sum(frac * diff(c(0, stats::plnorm(c(10, 20), m, s$phi), 1)))
  }, numeric(1))
  expect_equal(rtmb_lsas_log_selection(s, "lognormal", design, 1:2), log(expected))
})

test_that("length_strata() tabulates sampling fractions, excluding empty strata", {
  measured <- data.frame(
    unit = c(rep("a", 100), rep("b", 20)),
    length = c(rep(12, 10), rep(14, 90), rep(25, 20)),
    aged = c(rep(TRUE, 2), rep(FALSE, 8), rep(TRUE, 3), rep(FALSE, 87),
      rep(TRUE, 4), rep(FALSE, 16))
  )
  design <- length_strata(measured, "unit", c(0, 20, 40))
  expect_equal(design$frac, matrix(c(0.05, 0, 0, 0.2), 2))
  data <- data.frame(unit = c("b", "a"), length = c(25, 12))
  prepared <- .lsas_tmb_data(design, data, data$length,
    list(family = list(family = "gaussian")), "rtmb")
  s <- list(mu = c(25, 12), phi = 4)
  expected <- c(log(0.2) + pnorm(20, 25, 4, lower.tail = FALSE, log.p = TRUE),
    log(0.05) + pnorm(20, 12, 4, log.p = TRUE))
  expect_equal(rtmb_lsas_log_selection(s, "gaussian", prepared, 1:2), expected)
})

test_that("length_stratified checks its inputs", {
  sim <- lsas_sim(n_per_bin = 5)
  aged <- sim$aged[1:200, ]
  expect_error(sdmTMB(length ~ 1, data = aged, family = lognormal(),
    spatial = "off", control = sdmTMBcontrol(backend = "tmb"),
    length_stratified = sim$design), "rtmb")
  expect_error(lsas_fit(aged, family = Gamma(link = "log"),
    length_stratified = sim$design), "gaussian")
  expect_error(lsas_fit(aged, length_stratified = list()), "length_strata")

  measured <- data.frame(unit = "a", length = c(12, 13), aged = c(TRUE, FALSE))
  data <- data.frame(unit = "a", length = 12)
  prepare <- function(measured, data, breaks = c(0, 20, 40)) {
    .lsas_tmb_data(length_strata(measured, "unit", breaks), data, data$length,
      list(family = list(family = "gaussian")), "rtmb")
  }
  expect_error(prepare(measured, transform(data, unit = "b")), "must occur")
  expect_error(prepare(measured, transform(data, unit = NA_character_)), "missing")
  expect_error(prepare(transform(measured, unit = NA_character_), data), "missing")
  expect_error(prepare(measured, transform(data, length = 22)), "no aged fish")
  expect_error(prepare(transform(measured, aged = FALSE), data), "no aged fish")
  expect_error(prepare(transform(measured, aged = c(NA, TRUE)), data), "TRUE")
  expect_error(prepare(transform(measured, aged = 1), data), "TRUE")
  expect_error(prepare(transform(measured, length = c(NA, 1)), data), "finite")
  expect_error(prepare(measured, data, c(0, NA, 40)), "increasing")
  expect_error(prepare(measured, data, c(20, 0)), "increasing")
  expect_error(length_strata(measured, "year", c(0, 20)), "year")
  expect_error(length_strata(as.list(measured), "unit", c(0, 20)), "data frame")
})

test_that("tail bin probabilities and their derivatives remain finite", {
  # A finite interval deep in either tail of a standard normal at mu = 0.
  for (bounds in list(c(10, 11), c(-11, -10))) {
    objective <- function(p) rtmb_lsas_log_interval(bounds[1], bounds[2], p$mu, 1)
    obj <- RTMB::MakeADFun(objective, list(mu = 0), silent = TRUE)
    expected <- if (bounds[1] > 0) {
      pnorm(10, lower.tail = FALSE, log.p = TRUE) +
        log1p(-exp(pnorm(11, lower.tail = FALSE, log.p = TRUE) -
          pnorm(10, lower.tail = FALSE, log.p = TRUE)))
    } else {
      pnorm(-10, log.p = TRUE) + log1p(-exp(pnorm(-11, log.p = TRUE) -
        pnorm(-10, log.p = TRUE)))
    }
    expect_equal(obj$fn(0), expected, tolerance = 1e-10)
    # The same tape far into both tails and at the interval's centre.
    for (mu in c(-40, -20, 0, sum(bounds) / 2, 20, 40)) {
      gradient <- (obj$fn(mu + 1e-5) - obj$fn(mu - 1e-5)) / 2e-5
      expect_equal(as.numeric(obj$gr(mu)), gradient, tolerance = 1e-6)
      expect_true(all(is.finite(obj$he(mu))))
    }
  }
  # Open-ended bins, e.g., a gaussian fish far from its default start.
  design <- list(frac = matrix(c(0, 1), 1), breaks = c(0, 10, 20), unit_i = 1L)
  obj <- RTMB::MakeADFun(function(p) rtmb_lsas_log_selection(
    list(mu = p$mu, phi = 1), "gaussian", design, 1L), list(mu = 0), silent = TRUE)
  expect_equal(obj$fn(0), pnorm(10, lower.tail = FALSE, log.p = TRUE))
  expect_equal(obj$fn(-40), pnorm(50, lower.tail = FALSE, log.p = TRUE))
  expect_true(all(is.finite(obj$gr(-40))))
})

test_that("unsupported length_stratified diagnostics fail explicitly", {
  sim <- lsas_sim(n_per_bin = 5)
  fit <- lsas_fit(sim$aged, length_stratified = sim$design)
  for (type in c("mle-mvn", "mle-eb", "mle-mcmc", "response", "pearson", "deviance")) {
    expect_error(residuals(fit, type = type), "not yet supported.*length_stratified")
  }
  expect_error(simulate(fit), "not yet supported.*length_stratified")
  expect_error(dharma_residuals(NULL, fit), "not yet supported.*length_stratified")
  expect_s3_class(predict(fit, newdata = data.frame(age = 1:12)), "data.frame")
})
