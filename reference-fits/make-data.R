#!/usr/bin/env Rscript

# One-time generator for the frozen reference-fit fixtures. Run from the
# package root:
#
#   Rscript reference-fits/make-data.R
#
# Responses are simulated in base R from the model each case fits, with
# clearly identifiable parameters. The CSV files, not this script, are the
# fixtures: rerunning this script (e.g., under a different R version) may not
# reproduce them, and that is fine. Only rerun it deliberately, then re-record
# the reference.

out_dir <- file.path("reference-fits", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_fixture <- function(x, file, digits = 5L) {
  num <- vapply(x, is.double, logical(1))
  x[num] <- lapply(x[num], signif, digits = digits)
  utils::write.csv(x, file.path(out_dir, file), row.names = FALSE)
}

inv_link <- function(eta, link) {
  switch(link,
    identity = eta,
    log = exp(eta),
    logit = stats::plogis(eta),
    inverse = 1 / eta,
    cloglog = 1 - exp(-exp(eta))
  )
}

# Draws matching sdmTMB's parameterizations (see the `simulate` entries in
# R/rtmb-obs-families.R).
r_family <- function(family, mu, p = list()) {
  n <- length(mu)
  lognormal <- function(m, s) exp(stats::rnorm(n, log(m) - s^2 / 2, s))
  truncated_nb <- function(size, m) {
    log_nz <- stats::pnbinom(0, size = size, mu = m, lower.tail = FALSE, log.p = TRUE)
    stats::qnbinom(log(stats::runif(n)) + log_nz, size = size, mu = m,
      lower.tail = FALSE, log.p = TRUE)
  }
  mixture <- function(draw) {
    ifelse(stats::rbinom(n, 1, p$p_extreme) == 1, draw(mu * p$ratio), draw(mu))
  }
  gengamma <- function(m) {
    k <- p$Q^-2
    beta <- p$Q / p$sigma
    log_theta <- log(m) - lgamma((k * beta + 1) / beta) + lgamma(k)
    w <- log(stats::rgamma(n, k, 1))
    exp(w / beta + log_theta)
  }
  switch(family,
    gaussian = stats::rnorm(n, mu, p$sd),
    student = mu + p$sigma * stats::rt(n, p$df),
    Gamma = stats::rgamma(n, shape = p$shape, scale = mu / p$shape),
    lognormal = lognormal(mu, p$sdlog),
    gengamma = gengamma(mu),
    poisson = stats::rpois(n, mu),
    nbinom2 = stats::rnbinom(n, size = p$size, mu = mu),
    nbinom1 = stats::rnbinom(n, size = mu / p$phi, mu = mu),
    truncated_nbinom2 = truncated_nb(p$size, mu),
    truncated_nbinom1 = truncated_nb(mu / p$phi, mu),
    tweedie = {
      count <- stats::rpois(n, mu^(2 - p$power) / (p$phi * (2 - p$power)))
      stats::rgamma(n, shape = count * (2 - p$power) / (p$power - 1),
        scale = p$phi * (p$power - 1) * mu^(p$power - 1))
    },
    Beta = pmin(pmax(stats::rbeta(n, mu * p$phi, (1 - mu) * p$phi), 1e-6), 1 - 1e-6),
    gamma_mix = mixture(function(m) stats::rgamma(n, shape = p$shape, scale = m / p$shape)),
    lognormal_mix = mixture(function(m) lognormal(m, p$sdlog)),
    nbinom2_mix = mixture(function(m) stats::rnbinom(n, size = p$size, mu = m)),
    stop("No generator for ", family)
  )
}

# Positive-part parameters shared by non-spatial and delta generators.
pars <- list(
  gaussian = list(sd = 0.3),
  student = list(sigma = 0.3, df = 3),
  Gamma = list(shape = 4),
  lognormal = list(sdlog = 0.4),
  gengamma = list(sigma = 0.4, Q = -0.6),
  nbinom2 = list(size = 2),
  nbinom1 = list(phi = 1.5),
  truncated_nbinom2 = list(size = 2),
  truncated_nbinom1 = list(phi = 1.5),
  tweedie = list(power = 1.5, phi = 2),
  Beta = list(phi = 10),
  gamma_mix = list(shape = 6, p_extreme = 0.1, ratio = 5),
  lognormal_mix = list(sdlog = 0.3, p_extreme = 0.1, ratio = 5),
  nbinom2_mix = list(size = 5, p_extreme = 0.1, ratio = 5)
)

# Linear predictors by link that keep means valid over x in [-1, 1].
eta_for <- function(link, x) {
  switch(link,
    identity = 2 + 0.7 * x,
    log = 0.6 + 0.5 * x,
    inverse = 0.5 + 0.15 * x,
    logit = -0.2 + 1.2 * x,
    cloglog = -0.6 + 1.0 * x
  )
}

# Non-spatial family fixture ---------------------------------------------------

set.seed(20260930)
n <- 300L
fam <- data.frame(
  x = stats::runif(n, -1, 1),
  off = log(stats::runif(n, 0.6, 1.6))
)
col <- function(family, link, offset = FALSE) {
  paste0("y_", family, "_", link, if (offset) "_off")
}

# Families and links with a standard generator.
standard <- list(
  gaussian = c("identity", "log", "inverse"),
  student = c("identity", "log", "inverse"),
  Gamma = c("inverse", "identity", "log"),
  lognormal = c("identity", "log", "inverse"),
  gengamma = c("identity", "log", "inverse"),
  poisson = c("log", "identity"),
  nbinom2 = "log",
  nbinom1 = "log",
  truncated_nbinom2 = "log",
  truncated_nbinom1 = "log",
  tweedie = c("log", "identity"),
  gamma_mix = c("log", "identity", "inverse"),
  lognormal_mix = c("log", "identity", "inverse"),
  nbinom2_mix = "log"
)
for (f in names(standard)) {
  for (l in standard[[f]]) {
    mu <- inv_link(eta_for(l, fam$x), l)
    y <- r_family(f, mu, pars[[f]])
    # sdmTMB rejects negative responses with non-identity links; redraw them.
    while (f == "student" && l != "identity" && any(y <= 0)) {
      i <- y <= 0
      y[i] <- r_family(f, mu[i], pars[[f]])
    }
    fam[[col(f, l)]] <- y
  }
}
offset_cases <- list(gaussian = "identity", poisson = "log", Gamma = "inverse")
for (f in names(offset_cases)) {
  l <- offset_cases[[f]]
  # The inverse link needs a larger intercept to keep eta positive.
  eta <- eta_for(l, fam$x) + fam$off + if (l == "inverse") 0.7 else 0
  fam[[col(f, l, TRUE)]] <- r_family(f, inv_link(eta, l), pars[[f]])
}

# Binomial-type families: 0/1 responses, and successes out of `size` trials.
fam$size <- 5L
for (l in c("logit", "cloglog")) {
  p <- inv_link(eta_for(l, fam$x), l)
  fam[[col("binomial", l)]] <- stats::rbinom(n, 1, p)
  p_off <- inv_link(eta_for(l, fam$x) + fam$off, l)
  fam[[col("binomial", l, TRUE)]] <- stats::rbinom(n, 1, p_off)
  shape_phi <- 5
  pb <- stats::rbeta(n, p * shape_phi, (1 - p) * shape_phi)
  fam[[col("betabinomial", l)]] <- stats::rbinom(n, fam$size, pb)
}
fam$y_binomial_trials <- stats::rbinom(n, fam$size, inv_link(eta_for("logit", fam$x), "logit"))

# Beta and ordered beta.
mu <- inv_link(eta_for("logit", fam$x), "logit")
fam$y_Beta_logit <- r_family("Beta", mu, pars$Beta)
psi <- c(-1.8, 1.4)
u <- stats::runif(n)
eta <- eta_for("logit", fam$x)
p0 <- stats::plogis(psi[1] - eta)
p1 <- stats::plogis(eta - psi[2])
fam$y_ordbeta_logit <- ifelse(u < p0, 0, ifelse(u > 1 - p1, 1,
  r_family("Beta", mu, list(phi = 10))))

# Censored Poisson: right-censored at 8 (upper bound NA), otherwise exact.
cp <- stats::rpois(n, exp(1.5 + 0.5 * fam$x))
fam$y_censored_poisson_log <- pmin(cp, 8L)
fam$upr_censored_poisson <- ifelse(cp >= 8L, NA, cp)

# Delta models: encounter x positive. Standard delta offsets enter only the
# positive component; Poisson-link offsets enter both means.
delta_positive <- c("Gamma", "lognormal", "gengamma", "truncated_nbinom2",
  "truncated_nbinom1", "gamma_mix", "lognormal_mix", "Beta")
delta_std <- function(family, link1 = "logit", link2 = "log", offset = 0) {
  eta1 <- if (link1 == "logit") 0.2 + 0.8 * fam$x else -0.4 + 0.7 * fam$x
  eta2 <- switch(link2,
    log = 0.8 + 0.4 * fam$x,
    inverse = 0.4 + 0.1 * fam$x,
    identity = 2.5 + 0.8 * fam$x,
    logit = -0.3 + 0.8 * fam$x
  ) + offset
  present <- stats::rbinom(n, 1, inv_link(eta1, link1))
  present * r_family(family, inv_link(eta2, link2), pars[[family]])
}
for (f in delta_positive) {
  link2 <- if (f == "Beta") "logit" else "log"
  fam[[paste0("y_delta_", f)]] <- delta_std(f, link2 = link2)
}
fam$y_delta_Gamma_cloglog <- delta_std("Gamma", link1 = "cloglog")
fam$y_delta_Gamma_inverse <- delta_std("Gamma", link2 = "inverse")
fam$y_delta_Gamma_identity <- delta_std("Gamma", link2 = "identity")
fam$y_delta_lognormal_cloglog <- delta_std("lognormal", link1 = "cloglog")
for (f in c("Gamma", "lognormal", "truncated_nbinom2")) {
  fam[[paste0("y_delta_", f, "_off")]] <- delta_std(f, offset = fam$off)
}

delta_pl <- function(family, offset = 0) {
  n_dens <- exp(0.3 + 0.8 * fam$x + offset)
  w <- exp(0.2 + 0.3 * fam$x)
  p <- 1 - exp(-n_dens)
  stats::rbinom(n, 1, p) * r_family(family, n_dens * w / p, pars[[family]])
}
for (f in c("Gamma", "lognormal", "gengamma", "lognormal_mix")) {
  fam[[paste0("y_deltapl_", f)]] <- delta_pl(f)
}
for (f in c("Gamma", "lognormal")) {
  fam[[paste0("y_deltapl_", f, "_off")]] <- delta_pl(f, offset = fam$off)
}

# Censored NB: right-censored at 4 (upper bound Inf), otherwise exact. Drawn
# after the fixtures above so their draws are unchanged.
mu <- inv_link(eta_for("log", fam$x), "log")
for (f in c("nbinom2", "nbinom1")) {
  cnb <- r_family(f, mu, pars[[f]])
  fam[[col(paste0("censored_", f), "log")]] <- pmin(cnb, 4L)
  fam[[paste0("upr_censored_", f)]] <- ifelse(cnb >= 4L, Inf, cnb)
}

# Censored (beta-)binomial: catch out of `hooks`, censored when fewer than 3
# hooks are left (upper bound Inf, capped at `hooks`), otherwise exact.
fam$hooks <- 20L
censor_hooks <- function(y) ifelse(y >= fam$hooks - 2L, Inf, y)
for (l in c("logit", "cloglog")) {
  p <- inv_link(eta_for(l, fam$x) + 0.5, l)
  y <- stats::rbinom(n, fam$hooks, p)
  fam[[col("censored_binomial", l)]] <- y
  fam[[paste0("upr_censored_binomial_", l)]] <- censor_hooks(y)
  pb <- stats::rbeta(n, p * 10, (1 - p) * 10)
  y <- stats::rbinom(n, fam$hooks, pb)
  fam[[col("censored_betabinomial", l)]] <- y
  fam[[paste0("upr_censored_betabinomial_", l)]] <- censor_hooks(y)
}
# Per-hook Poisson-lognormal catch rate with at most one fish per hook.
eta <- eta_for("cloglog", fam$x) + 0.5 + stats::rnorm(n, 0, 0.5)
y <- stats::rbinom(n, fam$hooks, inv_link(eta, "cloglog"))
fam$y_censored_binomial_cloglog_re <- y
fam$upr_censored_binomial_cloglog_re <- censor_hooks(y)

write_fixture(fam, "family.csv")
write_fixture(data.frame(
  x = seq(-1, 1, length.out = 41),
  off = seq(-0.4, 0.4, length.out = 41),
  size = 5L
), "family-newdata.csv")

# Spatial fixture -----------------------------------------------------------------

set.seed(20261001)
years <- 2011:2015
n_per_year <- 100L
sp <- data.frame(
  X = stats::runif(n_per_year * length(years), 0, 10),
  Y = stats::runif(n_per_year * length(years), 0, 10),
  year = rep(years, each = n_per_year)
)
n_sp <- nrow(sp)

# A deterministic covariate surface, so nonlocal (diffusion) cases can
# evaluate it anywhere.
x_surface <- function(X, Y, year) {
  sin(X / 2.5) + cos(Y / 3) * 0.8 + 0.25 * (year - 2013) * cos(X / 4)
}
sp$x <- x_surface(sp$X, sp$Y, sp$year)
sp$z <- stats::runif(n_sp, -1, 1)
sp$off <- log(stats::runif(n_sp, 0.6, 1.6))
n_g <- 20L
sp$g <- factor(sprintf("g%02d", sample(n_g, n_sp, replace = TRUE)))

# Matern (nu = 1) field at the observation locations, optionally anisotropic
# (rotate by `angle`, then stretch one axis by `ratio`).
rfield <- function(range, sd, angle = 0, ratio = 1) {
  rot <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
  xy <- as.matrix(sp[, c("X", "Y")]) %*% rot
  xy[, 2] <- xy[, 2] * ratio
  d <- as.matrix(stats::dist(xy))
  kappa <- sqrt(8) / range
  kd <- kappa * d
  S <- sd^2 * ifelse(d == 0, 1, kd * besselK(kd, 1))
  as.vector(t(chol(S + diag(1e-8, n_sp))) %*% stats::rnorm(n_sp))
}
# Spatiotemporal field with AR1 (rho) correlation across years; each row
# takes its own year's value.
st_field <- function(range, sd, rho = 0) {
  eps <- rfield(range, sd)
  out <- numeric(n_sp)
  for (t in seq_along(years)) {
    if (t > 1) eps <- rho * eps + sqrt(1 - rho^2) * rfield(range, sd)
    i <- sp$year == years[t]
    out[i] <- eps[i]
  }
  out
}
yr_effect <- c(0, 0.3, -0.2, 0.4, 0.1)[match(sp$year, years)]

omega <- rfield(range = 4, sd = 0.6)
eps_ar1 <- st_field(range = 3, sd = 0.4, rho = 0.6)
sp$y_gauss <- 0.5 + yr_effect + 0.4 * sp$x + omega + eps_ar1 +
  stats::rnorm(n_sp, 0, 0.3)

omega_aniso <- rfield(range = 6, sd = 0.7, angle = pi / 6, ratio = 2)
sp$y_aniso <- 0.5 + 0.4 * sp$x + omega_aniso + stats::rnorm(n_sp, 0, 0.3)

eps_short <- st_field(range = 1.5, sd = 0.5)
sp$y_range <- 0.5 + yr_effect + 0.4 * sp$x + rfield(range = 7, sd = 0.6) +
  eps_short + stats::rnorm(n_sp, 0, 0.2)

zeta <- rfield(range = 5, sd = 0.4)
sp$y_svc <- 0.5 + (0.4 + zeta) * sp$x + omega + stats::rnorm(n_sp, 0, 0.3)

b0_t <- c(0.2, 0.6, 0.4, 0.9, 0.7)[match(sp$year, years)]
b1_t <- c(0.8, 0.5, 0.3, 0.6, 0.9)[match(sp$year, years)]
sp$y_tv <- b0_t + b1_t * sp$z + omega + stats::rnorm(n_sp, 0, 0.3)

re <- MASS::mvrnorm(n_g, c(0, 0), matrix(c(0.5^2, 0.5 * 0.5 * 0.3, 0.5 * 0.5 * 0.3, 0.3^2), 2))
gi <- as.integer(sp$g)
sp$y_re <- 0.5 + re[gi, 1] + (0.4 + re[gi, 2]) * sp$x + omega + stats::rnorm(n_sp, 0, 0.3)

sp$y_smooth <- 0.5 + sin(2.5 * sp$z) + omega + stats::rnorm(n_sp, 0, 0.3)
sp$y_breakpt <- 0.5 + 1.2 * pmin(sp$z, 0.2) + omega + stats::rnorm(n_sp, 0, 0.25)
sp$y_logistic <- 0.2 + 1.5 * stats::plogis((sp$z - 0) / 0.15) + omega +
  stats::rnorm(n_sp, 0, 0.25)

eps_count <- st_field(range = 3, sd = 0.4, rho = 0.6)
sp$y_nb2 <- r_family("nbinom2", exp(0.5 + 0.4 * sp$x + omega + eps_count + sp$off),
  list(size = 2))
sp$y_gamma <- r_family("Gamma", exp(0.5 + 0.4 * sp$x + omega), list(shape = 4))

# Delta: separate fields per component with different ranges (and different
# spatial and spatiotemporal ranges in component 1), and an offset in the
# positive component. Variants add the one structure a case estimates.
omega1 <- rfield(range = 5, sd = 0.7)
eps1 <- st_field(range = 3, sd = 0.4, rho = 0.5)
omega2 <- rfield(range = 2.5, sd = 0.5)
eps2 <- st_field(range = 2.5, sd = 0.3, rho = 0.6)
r_delta <- function(eta1, eta2) {
  stats::rbinom(n_sp, 1, stats::plogis(eta1)) *
    r_family("Gamma", exp(eta2 + sp$off), list(shape = 4))
}
sp$y_delta <- r_delta(0.2 + 0.5 * sp$x + omega1 + eps1, 1 + 0.3 * sp$z + omega2 + eps2)
sp$y_delta_aniso <- r_delta(
  0.4 + 0.5 * sp$x + rfield(range = 4, sd = 1.2, angle = pi / 6, ratio = 2),
  1 + 0.3 * sp$z + rfield(range = 4, sd = 0.8, angle = pi / 6, ratio = 2))
sp$y_delta_svc <- r_delta(
  0.2 + (0.5 + rfield(range = 5, sd = 1.0)) * sp$x + omega1,
  1 + 0.3 * sp$z + rfield(range = 5, sd = 0.3) * sp$x + omega2)
b0_t1 <- c(0.3, 0.9, 0.1, 0.6, -0.2)[match(sp$year, years)]
b0_t2 <- c(1.0, 1.4, 0.8, 1.2, 1.5)[match(sp$year, years)]
sp$y_delta_tv <- r_delta(b0_t1 + 0.5 * sp$x + omega1, b0_t2 + 0.3 * sp$z + omega2)
re1 <- stats::rnorm(n_g, 0, 0.6)
re2 <- stats::rnorm(n_g, 0, 0.4)
sp$y_delta_re <- r_delta(0.2 + 0.5 * sp$x + re1[gi] + omega1,
  1 + 0.3 * sp$z + re2[gi] + omega2)
n_dens <- exp(0.3 + 0.5 * sp$x + omega1 + eps1)
w <- exp(0.2 + 0.3 * sp$z + omega2)
p <- 1 - exp(-n_dens)
sp$y_deltapl <- stats::rbinom(n_sp, 1, p) *
  r_family("lognormal", n_dens * w / p, list(sdlog = 0.4))

# Multi-family: Gaussian and Poisson rows sharing one spatial field.
sp$dist <- rep(c("gauss", "pois"), length.out = n_sp)
sp$y_mf <- ifelse(sp$dist == "gauss",
  stats::rnorm(n_sp, 1 + 0.4 * sp$x + omega, 0.3),
  stats::rpois(n_sp, exp(0.5 + 0.4 * sp$x + omega)))

# Nonlocal: responds to a spatially smoothed (Gaussian kernel) covariate and
# to the previous year's covariate.
kernel_x <- vapply(seq_len(n_sp), function(i) {
  g <- expand.grid(X = seq(0, 10, 0.5), Y = seq(0, 10, 0.5))
  wt <- exp(-((g$X - sp$X[i])^2 + (g$Y - sp$Y[i])^2) / (2 * 1.5^2))
  sum(wt * x_surface(g$X, g$Y, sp$year[i])) / sum(wt)
}, numeric(1))
x_lag <- x_surface(sp$X, sp$Y, pmax(sp$year - 1, 2011))
sp$y_nl <- 0.5 + 0.2 * sp$x + 0.6 * kernel_x + 0.4 * x_lag + stats::rnorm(n_sp, 0, 0.2)

write_fixture(sp, "spatial.csv")

# The nonlocal covariate grid is built in cases.R from x_surface() (keep the
# two definitions in sync).

nd <- expand.grid(X = seq(0.5, 9.5, length.out = 8), Y = seq(0.5, 9.5, length.out = 8),
  year = years)
nd$x <- x_surface(nd$X, nd$Y, nd$year)
nd$z <- rep(seq(-0.9, 0.9, length.out = 16), length.out = nrow(nd))
nd$off <- rep(log(c(0.8, 1, 1.25)), length.out = nrow(nd))
nd$g <- sprintf("g%02d", rep(seq_len(n_g), length.out = nrow(nd)))
nd$dist <- rep(c("gauss", "pois"), length.out = nrow(nd))
nd$area <- 1
write_fixture(nd, "spatial-newdata.csv")

# Frozen mesh: vertices and triangles, rebuilt with fmesher in the cases.
devtools::load_all(quiet = TRUE)
mesh <- make_mesh(sp, c("X", "Y"), cutoff = 0.7)
loc <- mesh$mesh$loc
utils::write.csv(data.frame(X = sprintf("%.17g", loc[, 1]), Y = sprintf("%.17g", loc[, 2])),
  file.path(out_dir, "mesh-loc.csv"), row.names = FALSE, quote = FALSE)
tv <- mesh$mesh$graph$tv
utils::write.csv(data.frame(v1 = tv[, 1], v2 = tv[, 2], v3 = tv[, 3]),
  file.path(out_dir, "mesh-tv.csv"), row.names = FALSE)
cat("Mesh vertices:", mesh$mesh$n, "\n")

# Areal adjacency: frozen Ohio county edge list (used with ohio_df).
domain <- make_areal_domain(ohio_sf, id_column = "county")
W <- as.matrix(domain$W_raw)
idx <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
utils::write.csv(data.frame(
  from = domain$unit_names[idx[, 1]],
  to = domain$unit_names[idx[, 2]]
), file.path(out_dir, "ohio-edges.csv"), row.names = FALSE)
