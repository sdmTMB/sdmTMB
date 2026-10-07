# Linear predictor pieces for fitted or projected rows (`rows` comes from
# `rtmb_row_inputs()`). Fitting and projection share these helpers so each
# term is translated once. Each piece is an n-by-component matrix, except
# `zeta` (n by SVC by component).
rtmb_linear_predictors <- function(par, theta, effects, prepared, rows) {
  n <- nrow(rows$X[[1L]])
  diffusion <- if (prepared$diffusion$n_terms > 0L) {
    rtmb_diffusion_terms(par, prepared$diffusion, rows)
  } else {
    matrix(0, n, 0L)
  }
  parts <- lapply(seq_len(prepared$n_m), function(m) {
    rtmb_linear_predictor(par, theta, effects, prepared, rows, m, diffusion)
  })
  out <- lapply(stats::setNames(nm = setdiff(names(parts[[1L]]), "zeta")),
    function(name) unname(do.call(cbind, lapply(parts, `[[`, name))))
  n_z <- if (prepared$svc) ncol(rows$z) else 0L
  out$zeta <- array(do.call(cbind, lapply(parts, `[[`, "zeta")),
    dim = c(n, n_z, prepared$n_m))
  out$diffusion <- diffusion
  out
}

rtmb_linear_predictor <- function(par, theta, effects, prepared, rows, m,
                                  diffusion) {
  n <- nrow(rows$X[[m]])
  zero <- rep(0, n)
  b <- if (m == 1L) par$b_j else par$b_j2

  fixed <- rtmb_product(rows$X[[m]], b)
  # Diffusion coefficients are the last `b` entries; their `X` columns are 0.
  n_terms <- ncol(diffusion)
  for (term in seq_len(n_terms)) {
    fixed <- fixed + b[[length(b) - n_terms + term]] * diffusion[, term]
  }
  if (prepared$threshold != "none") {
    fixed <- fixed + rtmb_threshold(rows$X_threshold, theta$threshold,
      prepared$threshold, m)
  }
  smooth <- zero
  if (prepared$smooths) {
    for (s in seq_along(rows$Zs)) {
      smooth <- smooth + rtmb_product(rows$Zs[[s]],
        effects$b_smooth[prepared$smooth_index[[s]], m])
    }
    smooth <- smooth + rtmb_product(rows$Xs, par$bs[, m])
  }
  rw <- zero
  if (prepared$time_varying) {
    for (k in seq_len(ncol(rows$X_rw))) {
      rw <- rw + rows$X_rw[, k] * effects$b_rw_t[rows$time, k, m]
    }
  }
  iid <- zero
  if (prepared$iid_re[[m]] && rows$include_iid) {
    iid <- rtmb_product(Matrix::t(rows$Zt[[m]]), effects$re_b_pars[, m])
  }
  omega <- zero
  if (prepared$spatial[[m]]) {
    omega <- rtmb_product(rows$A_station, effects$omega_s[, m])[
      rows$station_index]
  }
  epsilon <- zero
  if (prepared$temporal[[m]]) {
    epsilon <- rtmb_station_values(rows,
      lapply(seq_len(prepared$n_t), function(t) effects$epsilon_st[, t, m]))
  }
  zeta <- matrix(0, n, 0L)
  svc <- zero
  if (prepared$svc) {
    zeta <- do.call(cbind, lapply(seq_len(ncol(rows$z)), function(z) {
      rtmb_product(rows$A_station, effects$zeta_s[, z, m])[
        rows$station_index]
    }))
    for (z in seq_len(ncol(rows$z))) svc <- svc + zeta[, z] * rows$z[, z]
  }
  fe <- fixed + rows$offset * rows$offset_applies[, m] + smooth + rw + iid
  rf <- omega + epsilon + svc
  list(fixed = fixed, smooth = smooth, rw = rw, iid = iid,
    omega = omega, epsilon = epsilon, zeta = zeta, svc = svc,
    fe = fe, rf = rf, eta = fe + rf)
}

# Plain (AD or numeric) vector from a matrix result. Sparse products of plain
# numbers return Matrix objects, for example when reporting.
rtmb_as_vector <- function(out) {
  if (methods::is(out, "Matrix")) out <- as.matrix(out)
  unname(drop(out))
}

rtmb_product <- function(A, x) rtmb_as_vector(A %*% x)

# Project per-time vertex values to stations, then index by station and time.
# `values` is a list of one vector per time step, projected in one product.
rtmb_station_values <- function(rows, values) {
  at_station <- rows$A_station %*% do.call(cbind, values)
  if (methods::is(at_station, "Matrix")) at_station <- as.matrix(at_station)
  at_station[rows$station_index + nrow(at_station) * (rows$time - 1L)]
}

# Threshold covariate effect for component `m`, as in the C++
# `linear_threshold()` and `logistic_threshold()`.
rtmb_threshold <- function(x, threshold, type, m) {
  s <- lapply(threshold, `[[`, m)
  switch(type,
    # slope * min(x, cut)
    breakpt = s$s_slope * (x + s$s_cut - abs(x - s$s_cut)) / 2,
    logistic = s$s_max * RTMB::plogis(log(19) * (x - s$s50) / (s$s95 - s$s50)))
}

# Combined link- and response-scale projections for two-component models,
# following the C++ `combined_link_value()` and `combined_response_value()`.
rtmb_combined_projection <- function(projected, theta, prepared) {
  "[<-" <- RTMB::ADoverload("[<-")
  rows <- prepared$proj
  n <- nrow(projected$eta)
  fe <- eta <- response <- rep(0, n)
  response_eta <- if (prepared$pop_pred) projected$fe else projected$eta
  for (f in unique(rows$family_id)) {
    family <- prepared$families[[f]]
    i <- which(rows$family_id == f)
    link <- function(x1, x2) {
      switch(family$combine,
        single = x1,
        delta = rtmb_link(
          rtmb_inverse_link(x1, family$link[[1L]]) *
            rtmb_inverse_link(x2, family$link[[2L]]),
          family$link[[2L]]),
        poisson_link = x1 + x2)
    }
    fe[i] <- link(projected$fe[i, 1L], projected$fe[i, 2L])
    eta[i] <- link(projected$eta[i, 1L], projected$eta[i, 2L])
    response[i] <- switch(family$combine,
      single = rtmb_component_mean(response_eta[i, 1L], family, 1L, theta),
      delta = rtmb_component_mean(response_eta[i, 1L], family, 1L, theta) *
        rtmb_component_mean(response_eta[i, 2L], family, 2L, theta),
      poisson_link = exp(response_eta[i, 1L] + response_eta[i, 2L]))
  }
  list(fe = fe, eta = eta, response = response)
}

# Mixture families project the mean of both components: the positive
# component's mean mu becomes (1 - p) mu + p mu ratio, in both the full
# (`eta`) and population-level (`fe`) predictions.
rtmb_mixture_projection <- function(projected, theta, prepared) {
  "[<-" <- RTMB::ADoverload("[<-")
  m <- prepared$n_m
  rows <- prepared$proj
  p <- theta$p_extreme
  scale <- 1 - p + p * theta$mix_ratio
  for (f in unique(rows$family_id)) {
    i <- which(rows$active[, m] & rows$family_id == f)
    link <- prepared$families[[f]]$link[[m]]
    for (x in c("fe", "eta")) {
      projected[[x]][i, m] <- rtmb_link(
        rtmb_inverse_link(projected[[x]][i, m], link) * scale, link)
    }
  }
  projected
}

# Area-weighted totals, weighted averages, and effective area occupied (EAO)
# by projected time step, following the C++ derived-quantity block. The
# `eps_index` term makes the gradient of the marginal objective with respect
# to `eps_index` the bias-corrected output (Thorson and Kristensen 2016).
rtmb_derived_indices <- function(par, theta, prepared, projected) {
  "[<-" <- RTMB::ADoverload("[<-")
  rows <- prepared$proj
  requested <- prepared$derived
  if (prepared$n_m > 1L) {
    mu <- projected$combined$response
  } else {
    mu <- rep(0, nrow(projected$eta))
    for (f in unique(rows$family_id)) {
      family <- prepared$families[[f]]
      family$link[[1L]] <- prepared$index$link
      i <- which(rows$family_id == f)
      mu[i] <- rtmb_component_mean(projected$eta[i, 1L], family, 1L, theta)
      # Prototype: for cloglog (censored) betabinomial, the mean per-hook
      # catch rate E(-log(1 - p)) = digamma(phi) - digamma(phi * (1 - pbar))
      # instead of exp(eta). It matched an R-side check and works with bias
      # correction, but the bias-corrected index was about 4x slower on the
      # hook-competition article's grid, so it isn't exposed.
      # See the article for an R-side calculation from predict().
      # phi <- theta$phi[[family$phi]]
      # b <- phi * exp(-exp(projected$eta[i, 1L]))
      # mu[i] <- rtmb_digamma(phi) - rtmb_digamma(b)
    }
  }
  n_t <- prepared$n_t
  by_time <- lapply(seq_len(n_t), function(t) which(rows$time == t))
  time_sum <- function(x) {
    out <- rep(0, n_t)
    for (t in seq_len(n_t)) {
      if (length(by_time[[t]])) out[t] <- sum(x[by_time[[t]]])
    }
    out
  }
  # Time steps that are excluded or have no weighted rows report 0, as in C++.
  include <- prepared$index$time_include
  has_area <- vapply(by_time, function(i) any(rows$area[i] != 0), logical(1))
  out <- list(nll = 0)
  out$total <- time_sum(mu * rows$area)
  out$link_total <- log(out$total)
  # Bias-correction terms tag the per-time sums, whose Hessians are sparse,
  # and form any ratio after the tags so the sparse-plus-low-rank inner
  # Hessian stays sparse (see `Tag()` in TMB's newton.hpp).
  tag <- utils::getFromNamespace("Tag", "RTMB")
  eps_t <- if (length(par$eps_index)) which(include) else integer(0)
  if (requested[["total"]]) {
    for (t in eps_t) out$nll <- out$nll + par$eps_index[t] * tag(out$total[t])
  }
  if (requested[["weighted_avg"]]) {
    # One column per weighted vector (x and y for COG), stacked by time within
    # column; `eps_index` follows the same order.
    weight <- as.matrix(rows$weight)
    keep <- which(include & has_area)
    weighted_avg <- rep(0, n_t * ncol(weight))
    for (j in seq_len(ncol(weight))) {
      weighted_sum <- time_sum(weight[, j] * mu * rows$area)
      k <- (j - 1L) * n_t
      weighted_avg[k + keep] <- weighted_sum[keep] / out$total[keep]
      for (t in intersect(eps_t, keep)) {
        out$nll <- out$nll +
          par$eps_index[k + t] * tag(weighted_sum[t]) / tag(out$total[t])
      }
    }
    out$weighted_avg <- weighted_avg
  }
  if (requested[["eao"]]) {
    # eao = total / mean_dens = sum(area * mu)^2 / sum(area * mu^2)
    sum_dens2 <- time_sum(rows$area * mu * mu)
    keep <- include & has_area
    mean_dens <- eao <- log_eao <- rep(0, n_t)
    mean_dens[has_area] <- sum_dens2[has_area] / out$total[has_area]
    eao[keep] <- out$total[keep] / mean_dens[keep]
    log_eao[keep] <- log(eao[keep])
    out$mean_dens <- mean_dens
    out$eao <- eao
    out$log_eao <- log_eao
    for (t in intersect(eps_t, which(keep))) {
      out$nll <- out$nll + par$eps_index[t] *
        tag(out$total[t]) * tag(out$total[t]) / tag(sum_dens2[t])
    }
  }
  out
}

# AD-safe digamma for x > 0 (for the commented-out mean rate above) from
# psi(x) = psi(x + 6) - sum 1 / (x + j) and the asymptotic series at x + 6
# (relative error about 1e-12).
# rtmb_digamma <- function(x) {
#   shift <- 0
#   for (j in 0:5) shift <- shift + 1 / (x + j)
#   z <- x + 6
#   w <- 1 / (z * z)
#   log(z) - 0.5 / z -
#     w * (1 / 12 - w * (1 / 120 - w * (1 / 252 - w * (1 / 240 - w / 132)))) -
#     shift
# }
