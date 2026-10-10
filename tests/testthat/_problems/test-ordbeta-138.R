# Extracted from test-ordbeta.R:138

# test -------------------------------------------------------------------------
skip_on_cran()
set.seed(3)
n <- 300
x <- rnorm(n)
eta <- 0.2 + 0.7 * x
psi <- c(-1.0, 1.2)
p0 <- plogis(psi[1] - eta)
p1 <- plogis(eta - psi[2])
u <- runif(n)
y <- rbeta(n, plogis(eta) * 5, (1 - plogis(eta)) * 5)
y[u < p0] <- 0
y[u > 1 - p1] <- 1
d <- data.frame(y = y, x = x, year = 1L)
for (backend in c("tmb", "rtmb")) {
    fit <- sdmTMB(y ~ x, data = d, family = ordbeta(), spatial = "off",
      time = "year", control = sdmTMBcontrol(backend = backend))
    cuts <- ordbeta_cutpoints(get_pars(fit)$psi)
    eta_hat <- predict(fit)$est
    mean_y <- plogis(eta_hat - cuts[2]) +
      (1 - plogis(cuts[1] - eta_hat) - plogis(eta_hat - cuts[2])) *
      plogis(eta_hat)
    p <- predict(fit, type = "response")
    expect_equal(p$est, mean_y, tolerance = 1e-6)
    expect_equal(residuals(fit, type = "response"), y - mean_y,
      tolerance = 1e-6)
    p <- predict(fit, newdata = d, return_tmb_object = TRUE)
    ind <- get_index(p, area = 1 / n, bias_correct = FALSE)
    expect_equal(ind$est, mean(mean_y), tolerance = 1e-6)
  }
