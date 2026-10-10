# Extracted from test-rtmb-families.R:199

# prequel ----------------------------------------------------------------------
rtmb_family_data <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  mu <- exp(0.5 + 0.4 * d$z)
  d$tw <- ifelse(runif(n) < 0.3, 0, rgamma(n, 2, 2 / mu))
  d$cont <- 1 + d$z + 0.5 * rt(n, 4)
  d$prop <- rbeta(n, 2, 3)
  d$ordprop <- d$prop
  d$ordprop[seq_len(n / 10)] <- 0
  d$ordprop[n / 10 + seq_len(n / 10)] <- 1
  d$count <- pmax(rnbinom(n, size = 2, mu = mu), 1)
  d$count0 <- rnbinom(n, size = 2, mu = mu)
  d$pos <- rgamma(n, 2, 2 / mu)
  d$mixpos <- ifelse(runif(n) < 0.1, 5, 1) * d$pos
  d$trials <- rep(c(5, 10, 3), length.out = n)
  d$succ <- rbinom(n, d$trials, rbeta(n, 2, 3))
  d$hurdle_count <- ifelse(runif(n) < 0.4, 0, d$count)
  d$hurdle_prop <- ifelse(runif(n) < 0.4, 0, d$prop)
  d$hurdle_pos <- ifelse(runif(n) < 0.4, 0, d$pos)
  d
}
censpois_oracle <- function(lambda, L, U) {
  lp <- if (is.infinite(U)) {
    if (L == 0) 0 else ppois(L - 1, lambda, lower.tail = FALSE, log.p = TRUE)
  } else {
    ld <- dpois(L:U, lambda, log = TRUE)
    max(ld) + log(sum(exp(ld - max(ld))))
  }
  r <- function(k) {
    if (k < 0 || is.infinite(k)) 0 else exp(dpois(k, lambda, log = TRUE) - lp)
  }
  g <- r(L - 1) - r(U)
  dg <- (r(L - 2) - r(L - 1)) - (r(U - 1) - r(U)) - g^2
  c(value = lp, gradient = lambda * g, hessian = lambda * g + lambda^2 * dg)
}

# test -------------------------------------------------------------------------
cases <- rbind(
    c(1, 100, Inf), c(1, 100, 102), c(1000, 0, 2), c(3, 0, Inf), c(3, 0, 5),
    c(3, 2, 6), c(1e-3, 50, Inf), c(50, 3, 3), c(5, 10, 500),
    c(1000, 10, 200), c(300, 200, 400), c(1e-8, 2, Inf), c(400, 2, Inf))
expect_equal(rtmb_dcenspois(cases[, 2], cases[, 1], cases[, 3]),
    apply(cases, 1, function(x) censpois_oracle(x[1], x[2], x[3])[[1]]),
    tolerance = 1e-12)
for (k in seq_len(nrow(cases))) {
    lambda <- cases[k, 1]
    L <- cases[k, 2]
    d <- data.frame(y = L, o = log(lambda))
    fit <- sdmTMB(y ~ 1, offset = "o", data = d, spatial = "off",
      do_fit = FALSE, family = censored_poisson(),
      censored_upper = upr[k])
    for (backend in c("tmb", "rtmb")) {
      obj <- make_sdmTMB_adfun(fit$tmb_data, fit$tmb_params, fit$tmb_map,
        backend = backend)
      for (b in c(0, 0.3, -0.4)) {
        expected <- censpois_oracle(lambda * exp(b), L, cases[k, 3])
        got <- -c(obj$fn(b), obj$gr(b), obj$he(b))
        expect_lt(max(abs(got - expected) / pmax(1, abs(expected))), 1e-8,
          label = paste(backend, k, b))
      }
    }
  }
