#############################################################
# Simulation / comparison tests of RTMB vs sdmTMB on dev-preferential
#
# A square grid with delta-lognormal catches and Conn-style areal sampling,
# simulated *within* the joint model class that sdmTMB fits:
#   catch fields      omega_m, epsilon_mt (SPDE GMRFs on the estimation mesh)
#   encounter         logit p_gt   = beta1_t + omega1_g + eps1_gt
#   positive mean     log mu_gt    = beta2_t + omega2_g + eps2_gt
#   expected catch    E[y_gt]      = p_gt * mu_gt   (lognormal bias-corrected)
#   sampling          logit P(R_gt = 1) = gamma_t + b_t * h_gt + xi_g
#   target            h_gt = log p_gt - log plogis(beta1_t) + omega2_g + eps2_gt
#                     (= sdmTMB target = "fields" for a delta model)
#   catches           n_hauls lognormal-delta draws in each sampled cell-year
# so the true index and the true b_t are known exactly, and b_t is on the
# same scale as sdmTMB's estimate.
#
# Tests:
#   T1: Fits the delta-lognormal model with no sampling component in both sdmTMB
#   and the hand-coded RTMB model, and checks that they give the same likelihood,
#   parameters and index
#
#   T2: Adds the sampling model with one preference coefficient b for all years,
#   and checks that both implementations give the same joint likelihood, b,
#   sampling intercepts, sampling-field parameters and index.
#
#   T3:  Lets the preference coefficient follow a random walk over years, and
#   checks that both implementations give the same b_t path, sigma_b and index
#
#   T4: Fixes b at 0, so the sampling locations carry no information about catch,
#   and checks that both implementations reproduce the catch-only (T1) estimates exactly.
#
#   T5: Starts from your original RTMB model and changes one assumption at a time until
#   it matches sdmTMB, all on the same data set, to show which choices change b and the
#   index. In order, the changes are:
#   - annual sampling intercepts
#   - sdmTMB’s target, which adds the encounter fields
#   - a separate kappa for each field
#   - the bias-corrected lognormal
#   - first-year b as a fixed parameter, with your priors removed
#
#   T6: Simulates replicate data sets under three scenarios: no preference,
#   constant preference, and declining preference. It fits each one with sdmTMB naive,
#   sdmTMB constant b, sdmTMB random walk, and your original RTMB model. It measures
#   how well each recovers the true b_t, the absolute index and the index trend,
#   plus index coverage for the sdmTMB fits.
#############################################################

library(RTMB)
library(sdmTMB)
library(fmesher)
library(Matrix)
library(dplyr)
library(tidyr)

#############################################################
# Settings
#############################################################
out_dir <- "simple_sim"
dir.create(out_dir, showWarnings = FALSE)

run_tests <- c("T1", "T2", "T3", "T4", "T5", "T6")

# fast = TRUE: 20 x 20 grid, 6 years, ~160 mesh vertices (2,400 frame rows).
# fast = FALSE: 30 x 30 grid, 10 years, ~330 vertices (9,000 frame rows).
fast <- TRUE # The `fast` version still takes an hour or so
n_side <- if (fast) 20 else 30 # grid is n_side x n_side cells (area 1)
n_t <- if (fast) 6 else 10
verbose_fits <- FALSE # print sdmTMB optimizer progress
years <- 2011:(2011 + n_t - 1)

# catch model truth
true <- list(
  beta1 = qlogis(0.6) + seq(0, -0.4, length.out = n_t), # encounter, declining
  beta2 = log(5) + seq(0, -0.5, length.out = n_t), # positive mean, declining
  range = c(10, 10), # cells (comp 1, comp 2)
  sigma_O = c(0.8, 0.6),
  sigma_E = c(0.3, 0.3),
  phi = 0.8 # SD of log(positive catch)
)

# sampling truth
samp <- list(
  p_sampled = if (fast) 0.25 else 0.15, # baseline proportion of cells sampled per year
  effort = "constant", # "constant" or "rising" (0.08 -> 0.25)
  range_xi = 8,
  sigma_xi = 0.5,
  n_hauls = if (fast) 2 else 1 # catch observations per sampled cell-year
)

# b_t scenarios: one function of t per scenario
b_scenarios <- list(
  b0       = function(t) rep(0, length(t)),
  const1   = function(t) rep(1, length(t)),
  decline  = function(t) 1.5 * (n_t - t) / (n_t - 1) # 1.5 -> 0, as in 00_
)

scenario_for_equiv <- "decline"
seed_for_equiv <- 1
t6_scenarios <- names(b_scenarios)
t6_reps <- 1:5

tol_nll <- 1e-3
tol_par <- 2e-3
tol_lin <- 5e-3
tol_idx <- 2e-3

#############################################################
# Grid and mesh
#############################################################
cells <- expand.grid(X = seq(0.5, n_side - 0.5, 1), Y = seq(0.5, n_side - 0.5, 1))
cells$cell_id <- seq_len(nrow(cells))
cells$area <- 1
cells <- cells[, c("cell_id", "X", "Y", "area")]
year_lu <- data.frame(year = years, t = seq_len(n_t))

fm <- fmesher::fm_mesh_2d_inla(
  loc = as.matrix(cells[, c("X", "Y")]),
  max.edge = c(2, 6), offset = c(1, 6), cutoff = 1.6
)
mesh_cells <- sdmTMB::make_mesh(cells, c("X", "Y"), mesh = fm)
spde <- fmesher::fm_fem(fm)
A_grid <- fmesher::fm_basis(fm, loc = as.matrix(cells[, c("X", "Y")]))
cat("Cells:", nrow(cells), " mesh vertices:", fm$n, "\n")

#############################################################
# Simulation
#############################################################
# Draw from an SPDE GMRF with Matern range and marginal SD
rgmrf <- function(n, range, sigma) {
  kappa <- sqrt(8) / range
  tau2 <- 1 / (4 * pi * kappa^2 * sigma^2)
  Q <- tau2 * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  L <- Matrix::Cholesky(Matrix::forceSymmetric(Q), LDL = FALSE, perm = TRUE)
  z <- matrix(rnorm(fm$n * n), fm$n, n)
  x <- Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt")
  as.matrix(x)
}

simulate_data <- function(b_fun, seed) {
  set.seed(seed)
  b_t <- b_fun(seq_len(n_t))

  # catch fields at mesh vertices -> cells
  om1 <- as.vector(A_grid %*% rgmrf(1, true$range[1], true$sigma_O[1]))
  om2 <- as.vector(A_grid %*% rgmrf(1, true$range[2], true$sigma_O[2]))
  ep1 <- as.matrix(A_grid %*% rgmrf(n_t, true$range[1], true$sigma_E[1]))
  ep2 <- as.matrix(A_grid %*% rgmrf(n_t, true$range[2], true$sigma_E[2]))
  xi <- as.vector(A_grid %*% rgmrf(1, samp$range_xi, samp$sigma_xi))

  eta1 <- sweep(om1 + ep1, 2, true$beta1, "+") # n_cells x n_t
  eta2 <- sweep(om2 + ep2, 2, true$beta2, "+")
  dens <- plogis(eta1) * exp(eta2) # expected catch per unit area
  h <- log(plogis(eta1)) - matrix(log(plogis(true$beta1)), nrow(eta1), n_t, byrow = TRUE) +
    om2 + ep2

  p_t <- if (samp$effort == "rising") seq(0.08, 0.25, length.out = n_t) else
    rep(samp$p_sampled, n_t)
  gamma_t <- qlogis(p_t)
  lp <- sweep(sweep(h, 2, b_t, "*") + xi, 2, gamma_t, "+")
  R <- matrix(rbinom(length(lp), 1, plogis(lp)), nrow(lp), n_t)

  idx <- which(R == 1, arr.ind = TRUE)
  idx <- idx[rep(seq_len(nrow(idx)), each = samp$n_hauls), , drop = FALSE]
  e1 <- eta1[idx]
  e2 <- eta2[idx]
  z <- rbinom(nrow(idx), 1, plogis(e1))
  pos <- rlnorm(nrow(idx), e2 - true$phi^2 / 2, true$phi)
  dat <- data.frame(
    cell_id = cells$cell_id[idx[, 1]],
    X = cells$X[idx[, 1]], Y = cells$Y[idx[, 1]],
    t = idx[, 2], year = years[idx[, 2]],
    total_catch_wt_kg = z * pos,
    log_effort = 0
  )
  list(
    dat = dat, R = R, b_true = b_t, gamma_true = gamma_t,
    true_index = colSums(dens * cells$area),
    h_sd = apply(h, 2, sd)
  )
}

prep <- function(sim) {
  frame <- tidyr::crossing(cells, year_lu) |> as.data.frame()
  frame$sampled <- sim$R[cbind(match(frame$cell_id, cells$cell_id), frame$t)]
  mesh <- sdmTMB::make_mesh(sim$dat, c("X", "Y"), mesh = fm)
  c(sim, list(frame = frame, mesh = mesh))
}

#############################################################
# RTMB model
#############################################################
spec_original <- list(kappa = "shared", lnorm_bc = FALSE, pref = TRUE,
                      target = "pos", intercept = "pooled", b = "rw_random", priors = TRUE)
spec_sdmtmb <- list(kappa = "separate", lnorm_bc = TRUE, pref = TRUE,
                    target = "combined", intercept = "annual", b = "rw_fixed", priors = FALSE)

rtmb_data <- function(p) {
  y <- p$dat$total_catch_wt_kg
  pos <- which(y > 0)
  list(
    enc = as.numeric(y > 0), logy_pos = log(y[pos]), pos = pos,
    t_obs = p$dat$t, n_obs = nrow(p$dat), n_t = n_t,
    A_obs = fmesher::fm_basis(fm, loc = as.matrix(p$dat[, c("X", "Y")])),
    A_grid = A_grid,
    c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
    R_mat = p$R,
    obs_rows = lapply(seq_len(n_t), function(t) which(!is.na(p$R[, t]))),
    area = cells$area
  )
}

make_nll <- function(d, o) {
  function(p) {
    getAll(d, p)
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    Qk <- function(k) k^4 * c0 + 2 * k^2 * g1 + g2
    tau2 <- function(k, s) 1 / (4 * pi * k^2 * s^2)
    kap <- exp(log_kappa)
    k1 <- kap[1]
    k2 <- kap[length(kap)]
    sO <- exp(log_sigma_O)
    sE <- exp(log_sigma_E)
    s_obs <- exp(log_sigma_obs)
    Q1 <- Qk(k1)
    Q2 <- Qk(k2)
    nll <- 0
    nll <- nll - dgmrf(omega1, 0, Q1 * tau2(k1, sO[1]), log = TRUE)
    nll <- nll - dgmrf(omega2, 0, Q2 * tau2(k2, sO[2]), log = TRUE)
    for (t in 1:n_t) {
      nll <- nll - dgmrf(eps1[, t], 0, Q1 * tau2(k1, sE[1]), log = TRUE)
      nll <- nll - dgmrf(eps2[, t], 0, Q2 * tau2(k2, sE[2]), log = TRUE)
    }
    ii <- cbind(seq_len(n_obs), t_obs)
    eta1 <- beta1[t_obs] + as.vector(A_obs %*% omega1) + (A_obs %*% eps1)[ii]
    eta2 <- beta2[t_obs] + as.vector(A_obs %*% omega2) + (A_obs %*% eps2)[ii]
    nll <- nll - sum(dbinom_robust(enc, 1, eta1, log = TRUE))
    mlog <- eta2[pos]
    if (o$lnorm_bc) mlog <- mlog - s_obs^2 / 2
    nll <- nll - sum(dnorm(logy_pos, mlog, s_obs, log = TRUE)) + sum(logy_pos)

    G1 <- as.vector(A_grid %*% omega1)
    G2 <- as.vector(A_grid %*% omega2)
    GE1 <- A_grid %*% eps1
    GE2 <- A_grid %*% eps2
    log_index <- numeric(n_t)
    for (t in 1:n_t) {
      e1 <- beta1[t] + G1 + GE1[, t]
      e2 <- beta2[t] + G2 + GE2[, t]
      log_index[t] <- log(sum(plogis(e1) * exp(e2) * area))
    }
    REPORT(log_index)

    if (o$pref) {
      kx <- if (o$kappa == "separate") exp(log_kappa_xi) else k1
      nll <- nll - dgmrf(xi, 0, Qk(kx) * tau2(kx, exp(log_sigma_xi)), log = TRUE)
      xig <- as.vector(A_grid %*% xi)
      bt <- numeric(n_t)
      if (o$b == "zero") {
        bt <- 0 * bt
      } else if (o$b == "constant") {
        for (t in 1:n_t) bt[t] <- b
      } else if (o$b == "rw_fixed") {
        bt[1] <- b
        for (t in 2:n_t) bt[t] <- bt[t - 1] + b_dev[t - 1]
        nll <- nll - sum(dnorm(b_dev, 0, exp(log_sigma_b), log = TRUE))
      } else if (o$b == "rw_random") {
        bt <- b_path
        nll <- nll - dnorm(b_path[1], 0, 2, log = TRUE)
        for (t in 2:n_t) {
          nll <- nll - dnorm(b_path[t], b_path[t - 1], exp(log_sigma_b), log = TRUE)
        }
      }
      for (t in 1:n_t) {
        h <- G2 + GE2[, t]
        if (o$target == "combined") {
          e1 <- beta1[t] + G1 + GE1[, t]
          h <- h + log(plogis(e1)) - log(plogis(beta1[t]))
        }
        g <- if (o$intercept == "annual") gamma[t] else gamma[1]
        lp <- g + bt[t] * h + xig
        r <- obs_rows[[t]]
        nll <- nll - sum(dbinom_robust(R_mat[r, t], 1, lp[r], log = TRUE))
      }
      REPORT(bt)
      if (o$priors) {
        if (o$intercept == "pooled") nll <- nll - dnorm(gamma[1], 0, 1, log = TRUE)
        nll <- nll - dnorm(log_sigma_xi, log(0.5), 0.5, log = TRUE)
        if (o$b %in% c("rw_fixed", "rw_random")) {
          nll <- nll - dnorm(log_sigma_b, log(0.3), 0.5, log = TRUE)
        }
        if (o$b == "constant") nll <- nll - dnorm(b, 0, 2, log = TRUE)
      }
    }
    if (o$priors) {
      nll <- nll - sum(dnorm(log_sigma_O, log(0.5), 0.5, log = TRUE))
      nll <- nll - sum(dnorm(log_sigma_E, log(0.5), 0.5, log = TRUE))
      nll <- nll - dnorm(log_sigma_obs, log(0.5), 0.5, log = TRUE)
    }
    nll
  }
}

make_rtmb_pars <- function(d, o) {
  nv <- ncol(d$A_obs)
  p <- list(
    beta1 = rep(qlogis(mean(d$enc)), n_t),
    beta2 = rep(mean(d$logy_pos), n_t),
    log_sigma_O = rep(log(0.3), 2),
    log_sigma_E = rep(log(0.3), 2),
    log_kappa = rep(log(sqrt(8) / (n_side / 3)), if (o$kappa == "separate") 2 else 1),
    log_sigma_obs = log(0.5),
    omega1 = rep(0, nv), omega2 = rep(0, nv),
    eps1 = matrix(0, nv, n_t), eps2 = matrix(0, nv, n_t)
  )
  if (o$pref) {
    p$gamma <- rep(qlogis(mean(d$R_mat, na.rm = TRUE)),
                   if (o$intercept == "annual") n_t else 1)
    p$log_sigma_xi <- log(0.5)
    if (o$kappa == "separate") p$log_kappa_xi <- log(sqrt(8) / (n_side / 3))
    p$xi <- rep(0, nv)
    if (o$b %in% c("constant", "rw_fixed")) p$b <- 0
    if (o$b %in% c("rw_fixed", "rw_random")) p$log_sigma_b <- log(0.3)
    if (o$b == "rw_fixed") p$b_dev <- rep(0, n_t - 1)
    if (o$b == "rw_random") p$b_path <- rep(0, n_t)
  }
  p
}

rtmb_random <- function(o) {
  c("omega1", "omega2", "eps1", "eps2",
    if (o$pref) "xi",
    if (isTRUE(o$pref) && o$b == "rw_fixed") "b_dev",
    if (isTRUE(o$pref) && o$b == "rw_random") "b_path")
}

fit_rtmb <- function(d, o) {
  pars <- make_rtmb_pars(d, o)
  obj <- MakeADFun(make_nll(d, o), pars, random = rtmb_random(o), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 2000, eval.max = 4000))
  for (i in 1:2) {
    g <- max(abs(obj$gr(opt$par)))
    if (opt$convergence == 0 && g < 1e-3) break
    opt <- nlminb(opt$par, obj$fn, obj$gr,
                  control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-12))
  }
  obj$fn(opt$par)
  rp <- obj$report()
  list(obj = obj, opt = opt, rep = rp, pl = obj$env$parList(opt$par), o = o,
       max_grad = max(abs(obj$gr(opt$par))))
}

rtmb_summary <- function(h) {
  pl <- h$pl
  kap <- exp(pl$log_kappa)
  out <- list(
    nll = h$opt$objective,
    beta1 = pl$beta1, beta2 = pl$beta2,
    range1 = sqrt(8) / kap[1], range2 = sqrt(8) / kap[length(kap)],
    sigma_O1 = exp(pl$log_sigma_O[1]), sigma_O2 = exp(pl$log_sigma_O[2]),
    sigma_E1 = exp(pl$log_sigma_E[1]), sigma_E2 = exp(pl$log_sigma_E[2]),
    phi = exp(pl$log_sigma_obs),
    log_index = h$rep$log_index
  )
  if (h$o$pref) {
    out$gamma <- pl$gamma
    out$sigma_xi <- exp(pl$log_sigma_xi)
    out$range_xi <- if (!is.null(pl$log_kappa_xi)) sqrt(8) / exp(pl$log_kappa_xi) else out$range1
    out$b_t <- h$rep$bt
    if (!is.null(pl$b)) out$b_pref <- pl$b
    if (!is.null(pl$log_sigma_b)) out$sigma_b <- exp(pl$log_sigma_b)
  }
  out
}

#############################################################
# sdmTMB
#############################################################
fit_sdm <- function(p, coefficient = NULL, sampling_formula = sampled ~ 0 + factor(t)) {
  pref <- if (is.null(coefficient)) NULL else
    preferential_sampling(sampling_formula, data = p$frame, target = "fields",
                          coefficient = coefficient, spatial = "on")
  sdmTMB(
    total_catch_wt_kg ~ 0 + factor(t),
    data = p$dat, mesh = p$mesh, time = "t",
    family = delta_lognormal(),
    spatial = "on", spatiotemporal = list("iid", "iid"),
    preferential = pref,
    control = sdmTMBcontrol(backend = "rtmb"),
    silent = !verbose_fits
  )
}

sdm_summary <- function(fit, se = FALSE) {
  pick <- function(df, term) df$estimate[df$term == term][1]
  f1 <- tidy(fit, "fixed", model = 1)
  f2 <- tidy(fit, "fixed", model = 2)
  r1 <- tidy(fit, "ran_pars", model = 1)
  r2 <- tidy(fit, "ran_pars", model = 2)
  yt <- paste0("factor(t)", seq_len(n_t))
  nd <- tidyr::crossing(cells, year_lu) |> as.data.frame()
  ind <- get_index(fit, newdata = nd, area = nd$area, bias_correct = FALSE)
  ind <- ind[order(ind$t), ]
  out <- list(
    nll = fit$model$objective,
    beta1 = f1$estimate[match(yt, f1$term)],
    beta2 = f2$estimate[match(yt, f2$term)],
    range1 = pick(r1, "range"), range2 = pick(r2, "range"),
    sigma_O1 = pick(r1, "sigma_O"), sigma_O2 = pick(r2, "sigma_O"),
    sigma_E1 = pick(r1, "sigma_E"), sigma_E2 = pick(r2, "sigma_E"),
    phi = pick(r2, "phi"),
    log_index = ind$log_est
  )
  if (se) out$log_index_se <- ind$se
  if (!is.null(fit$preferential)) {
    fs <- tidy(fit, model = "sampling")
    rs <- tidy(fit, "ran_pars", model = "sampling")
    g_terms <- if (any(fs$term %in% yt)) yt else "(Intercept)"
    out$gamma <- fs$estimate[match(g_terms, fs$term)]
    out$b_pref <- pick(fs, "b_pref")
    out$sigma_xi <- pick(rs, "sigma_xi")
    out$range_xi <- pick(rs, "range_xi")
    if ("sigma_b_pref" %in% rs$term) {
      out$sigma_b <- pick(rs, "sigma_b_pref")
      rv <- tidy(fit, "ran_vals", model = "sampling")
      rv <- rv[rv$term == "b_pref_t", ]
      out$b_t <- rv$estimate[order(rv$time)]
      if (se) out$b_t_se <- rv$std.error[order(rv$time)]
    } else {
      out$b_t <- rep(out$b_pref, n_t)
      if (se) out$b_t_se <- rep(fs$std.error[fs$term == "b_pref"][1], n_t)
    }
  }
  out
}

cross_eval <- function(h, s) {
  pl <- h$pl
  pl$beta1 <- s$beta1
  pl$beta2 <- s$beta2
  pl$log_kappa <- log(sqrt(8) / c(s$range1, s$range2))[seq_along(pl$log_kappa)]
  pl$log_sigma_O <- log(c(s$sigma_O1, s$sigma_O2))
  pl$log_sigma_E <- log(c(s$sigma_E1, s$sigma_E2))
  pl$log_sigma_obs <- log(s$phi)
  if (h$o$pref) {
    pl$gamma <- s$gamma
    pl$log_sigma_xi <- log(s$sigma_xi)
    if (!is.null(pl$log_kappa_xi)) pl$log_kappa_xi <- log(sqrt(8) / s$range_xi)
    if (!is.null(pl$b)) pl$b <- s$b_pref
    if (!is.null(pl$log_sigma_b)) pl$log_sigma_b <- log(s$sigma_b)
  }
  fixed <- setdiff(names(pl), rtmb_random(h$o))
  h$obj$fn(unlist(pl[fixed]))
}

compare <- function(test, h, s, nll_at_sdm = NA_real_) {
  lin <- c("beta1", "beta2", "gamma")
  qs <- intersect(names(s), names(h))
  rows <- lapply(qs, function(q) {
    a <- as.numeric(h[[q]])
    b <- as.numeric(s[[q]])
    if (length(a) != length(b)) return(NULL)
    abs_d <- max(abs(a - b))
    rel_d <- max(abs(a - b) / pmax(abs(b), 0.1))
    tol <- if (q == "nll") tol_nll else if (q == "log_index") tol_idx else
      if (q %in% lin) tol_lin else tol_par
    stat <- if (q %in% c("nll", "log_index", lin)) abs_d else rel_d
    data.frame(test = test, quantity = q, n = length(a),
               max_abs_diff = abs_d, max_rel_diff = rel_d, pass = stat < tol)
  })
  out <- dplyr::bind_rows(rows)
  if (!is.na(nll_at_sdm)) {
    out <- rbind(out, data.frame(test = test, quantity = "rtmb_nll_at_sdmTMB_est",
                                 n = 1, max_abs_diff = abs(nll_at_sdm - s$nll), max_rel_diff = NA,
                                 pass = abs(nll_at_sdm - s$nll) < tol_nll))
  }
  out
}

run_equiv <- function(test, p, d, o, coefficient) {
  cat("\n==", test, "==\n")
  t0 <- Sys.time()
  fit <- fit_sdm(p, coefficient)
  t1 <- Sys.time()
  h <- fit_rtmb(d, o)
  t2 <- Sys.time()
  cat(sprintf("sdmTMB %.1f s, rtmb RTMB %.1f s; sanity all_ok = %s; rtmb max|grad| = %.2g\n",
              as.numeric(t1 - t0, units = "secs"), as.numeric(t2 - t1, units = "secs"),
              isTRUE(sanity(fit, silent = TRUE)$all_ok), h$max_grad))
  s <- sdm_summary(fit)
  t3 <- Sys.time()
  cat(sprintf("sdmTMB summary + get_index %.1f s\n", as.numeric(t3 - t2, units = "secs")))
  hs <- rtmb_summary(h)
  ce <- tryCatch(cross_eval(h, s), error = function(e) {
    message("cross_eval failed: ", conditionMessage(e)); NA_real_
  })
  res <- compare(test, hs, s, ce)
  print(res[, c("quantity", "n", "max_abs_diff", "max_rel_diff", "pass")], row.names = FALSE)
  list(res = res, fit = fit, rtmb = h, s = s, hs = hs)
}

# Recovery metrics against the known truth
summarise_fit <- function(label, s, sim) {
  li <- s$log_index
  lt <- log(sim$true_index)
  scaled_err <- (li - li[1]) - (lt - lt[1])
  has_b <- !is.null(s$b_t)
  data.frame(
    model = label,
    nll = s$nll,
    b_mean = if (has_b) mean(s$b_t) else NA,
    b_bias = if (has_b) mean(s$b_t - sim$b_true) else NA,
    b_rmse = if (has_b) sqrt(mean((s$b_t - sim$b_true)^2)) else NA,
    b_first = if (has_b) s$b_t[1] else NA,
    b_last = if (has_b) s$b_t[n_t] else NA,
    sigma_b = if (!is.null(s$sigma_b)) s$sigma_b else NA,
    log_index_bias = mean(li - lt), # absolute level
    trend_err_last = scaled_err[n_t], # log(est_T/est_1) - log(true_T/true_1)
    trend_rmse = sqrt(mean(scaled_err^2)),
    index_cover = if (!is.null(s$log_index_se))
      mean(abs(li - lt) < 1.96 * s$log_index_se) else NA
  )
}

#############################################################
# Quick look at one simulated data set
#############################################################
sim1 <- simulate_data(b_scenarios[[scenario_for_equiv]], seed_for_equiv)
p1 <- prep(sim1)
d1 <- rtmb_data(p1)
cat("Scenario", scenario_for_equiv, ": sampled cells per year:",
    paste(colSums(sim1$R), collapse = " "), "\n")
cat("Catch obs:", nrow(sim1$dat), " proportion positive:",
    round(mean(sim1$dat$total_catch_wt_kg > 0), 2), "\n")
cat("SD of target h by year:", paste(round(sim1$h_sd, 2), collapse = " "), "\n")

results <- list()

#############################################################
# T1-T4: equivalence
#############################################################
catch_only <- modifyList(spec_sdmtmb, list(pref = FALSE))
if ("T1" %in% run_tests) results$T1 <- run_equiv("T1 catch-only", p1, d1, catch_only, NULL)
if ("T2" %in% run_tests) results$T2 <- run_equiv("T2 joint, constant b", p1, d1,
                                                 modifyList(spec_sdmtmb, list(b = "constant")), "constant")
if ("T3" %in% run_tests) results$T3 <- run_equiv("T3 joint, RW b_t", p1, d1, spec_sdmtmb, "rw")
if ("T4" %in% run_tests) {
  results$T4 <- run_equiv("T4 joint, b = 0", p1, d1,
                          modifyList(spec_sdmtmb, list(b = "zero")), "zero")
  if (!is.null(results$T1)) {
    catch_q <- c("beta1", "beta2", "range1", "range2", "sigma_O1", "sigma_O2",
                 "sigma_E1", "sigma_E2", "phi", "log_index")
    z <- rbind(
      compare("T4 zero vs catch-only (sdmTMB)", results$T4$s[catch_q], results$T1$s[catch_q]),
      compare("T4 zero vs catch-only (rtmb)", results$T4$hs[catch_q], results$T1$hs[catch_q]))
    print(z[, c("test", "quantity", "max_abs_diff", "pass")], row.names = FALSE)
    results$T4$res <- rbind(results$T4$res, z)
  }
}

#############################################################
# T5: ladder from the original model to the sdmTMB model
#############################################################
if ("T5" %in% run_tests) {
  cat("\n== T5 ladder (scenario", scenario_for_equiv, ") ==\n")
  base <- spec_original
  ladder <- list(
    "0 original"                 = base,
    "1 + annual intercepts"      = modifyList(base, list(intercept = "annual")),
    "2 + combined-fields target" = modifyList(base, list(intercept = "annual", target = "combined")),
    "3 + separate kappas"        = modifyList(base, list(intercept = "annual", target = "combined",
                                                         kappa = "separate")),
    "4 + bias-corrected lnorm"   = modifyList(base, list(intercept = "annual", target = "combined",
                                                         kappa = "separate", lnorm_bc = TRUE)),
    "5 + RW b1 fixed, no priors" = spec_sdmtmb
  )
  t5 <- dplyr::bind_rows(lapply(names(ladder), function(lab) {
    h <- tryCatch(fit_rtmb(d1, ladder[[lab]]), error = function(e) NULL)
    if (is.null(h)) return(data.frame(model = lab))
    summarise_fit(lab, rtmb_summary(h), sim1)
  }))
  s_rw <- if (!is.null(results$T3)) results$T3$s else sdm_summary(fit_sdm(p1, "rw"))
  s_naive <- if (!is.null(results$T1)) results$T1$s else sdm_summary(fit_sdm(p1))
  t5 <- dplyr::bind_rows(t5,
                         summarise_fit("sdmTMB rw", s_rw, sim1),
                         summarise_fit("sdmTMB naive", s_naive, sim1))
  cat("True b_t:", paste(round(sim1$b_true, 2), collapse = " "), "\n")
  print(t5, row.names = FALSE, digits = 3)
  results$T5 <- t5
}

#############################################################
# T6: recovery over scenarios x replicates
#############################################################
if ("T6" %in% run_tests) {
  cat("\n== T6 recovery ==\n")
  t6 <- list()
  for (sc in t6_scenarios) {
    for (r in t6_reps) {
      sim <- simulate_data(b_scenarios[[sc]], seed = 1000 * match(sc, names(b_scenarios)) + r)
      p <- prep(sim)
      d <- rtmb_data(p)
      one <- function(label, expr) {
        s <- tryCatch(expr, error = function(e) {
          message(sc, " rep ", r, " ", label, ": ", conditionMessage(e)); NULL
        })
        if (is.null(s)) return(NULL)
        cbind(scenario = sc, rep = r, summarise_fit(label, s, sim))
      }
      t6[[length(t6) + 1]] <- one("sdmTMB naive", sdm_summary(fit_sdm(p), se = TRUE))
      t6[[length(t6) + 1]] <- one("sdmTMB constant", sdm_summary(fit_sdm(p, "constant"), se = TRUE))
      t6[[length(t6) + 1]] <- one("sdmTMB rw", sdm_summary(fit_sdm(p, "rw"), se = TRUE))
      t6[[length(t6) + 1]] <- one("RTMB original", rtmb_summary(fit_rtmb(d, spec_original)))
      cat(sc, "rep", r, "done\n")
    }
  }
  t6 <- dplyr::bind_rows(t6)
  saveRDS(t6, file.path(out_dir, "t6_recovery.rds"))
  cat("\nMeans by scenario x model:\n")
  print(t6 |> dplyr::group_by(scenario, model) |>
          dplyr::summarise(dplyr::across(c(b_mean, b_bias, b_rmse, sigma_b,
                                           log_index_bias, trend_err_last, trend_rmse, index_cover),
                                         ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop") |>
          as.data.frame(), digits = 3, row.names = FALSE)
  results$T6 <- t6
}

#############################################################
# Summary
#############################################################
equiv <- dplyr::bind_rows(lapply(results[intersect(names(results), c("T1", "T2", "T3", "T4"))],
                                 `[[`, "res"))
if (nrow(equiv)) {
  cat("\nEquivalence tests:", sum(equiv$pass), "of", nrow(equiv), "checks passed\n")
  if (any(!equiv$pass)) print(equiv[!equiv$pass, ], row.names = FALSE)
}
saveRDS(list(equiv = equiv, ladder = results$T5, recovery = results$T6,
             true = true, samp = samp),
        file.path(out_dir, "rtmb_vs_sdmtmb_simple.rds"))
