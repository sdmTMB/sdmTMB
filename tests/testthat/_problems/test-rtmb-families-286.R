# Extracted from test-rtmb-families.R:286

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
censbetabinom_oracle <- function(L, U, n, a, b) {
  k <- L:U
  l <- lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
  max(l) + log(sum(exp(l - max(l))))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("numDeriv")
set.seed(81)
n <- 30L
d <- data.frame(z = rnorm(n), hooks = sample(c(20, 50), n, replace = TRUE))
p <- 1 - exp(-exp(-1.5 + 0.4 * d$z))
d$y <- rbinom(n, d$hooks, rbeta(n, p * 10, (1 - p) * 10))
type <- seq_len(n) %% 3
upr <- ifelse(type == 0, Inf, ifelse(type == 1, d$y, pmin(d$y + 4, d$hooks)))
fit <- sdmTMB(y ~ z, data = d, weights = d$hooks, spatial = "off",
    do_fit = FALSE, family = censored_betabinomial(link = "cloglog"),
    censored_upper = upr)
expect_equal(fit$tmb_data$upr, ifelse(is.na(upr), d$hooks, upr))
