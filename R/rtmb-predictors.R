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
    omega <- rtmb_product(rows$A_rows, effects$omega_s[, m])
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
      rtmb_product(rows$A_rows, effects$zeta_s[, z, m])
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
# `values` is a list of one vector per time step.
rtmb_station_values <- function(rows, values) {
  at_station <- do.call(cbind, lapply(values, function(v) {
    rtmb_product(rows$A_station, v)
  }))
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
    fe[i] <- rtmb_combined_link(projected$fe[i, 1L], projected$fe[i, 2L],
      family)
    eta[i] <- rtmb_combined_link(projected$eta[i, 1L], projected$eta[i, 2L],
      family)
    response[i] <- switch(family$combine,
      single = rtmb_component_mean(response_eta[i, 1L], family, 1L, theta),
      delta = rtmb_component_mean(response_eta[i, 1L], family, 1L, theta) *
        rtmb_component_mean(response_eta[i, 2L], family, 2L, theta),
      poisson_link = exp(response_eta[i, 1L] + response_eta[i, 2L]))
  }
  list(fe = fe, eta = eta, response = response)
}

# Link-scale value of both components' combined mean, following the C++
# `combined_link_value()`. With a log positive link, log(p mu) is computed as
# log(p) + eta2, which stays finite when p underflows. `x2` is unused for a
# single-component family.
rtmb_combined_link <- function(x1, x2, family) {
  switch(family$combine,
    single = x1,
    delta = if (identical(family$link[[2L]], "log")) {
      rtmb_log_inverse_link(x1, family$link[[1L]]) + x2
    } else {
      rtmb_link(
        rtmb_inverse_link(x1, family$link[[1L]]) *
          rtmb_inverse_link(x2, family$link[[2L]]),
        family$link[[2L]])
    },
    poisson_link = x1 + x2)
}

# Mixture families project the mean of both components: the positive
# component's link-scale prediction becomes log((1 - p) mu + p mu ratio).
rtmb_mixture_eta <- function(eta, theta, prepared) {
  "[<-" <- RTMB::ADoverload("[<-")
  m <- prepared$n_m
  i <- which(prepared$proj$active[, m])
  p <- theta$p_extreme
  eta[i, m] <- log((1 - p) * exp(eta[i, m]) +
    p * exp(eta[i, m]) * theta$mix_ratio)
  eta
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
  eps_values <- list()
  if (requested[["total"]]) eps_values$total <- out$total
  if (requested[["weighted_avg"]]) {
    weighted_avg <- time_sum(rows$weight * mu * rows$area)
    keep <- include & has_area
    weighted_avg[keep] <- weighted_avg[keep] / out$total[keep]
    weighted_avg[!keep] <- 0
    out$weighted_avg <- eps_values$weighted_avg <- weighted_avg
  }
  if (requested[["eao"]]) {
    sum_dens <- time_sum(mu)
    mean_dens <- rep(0, n_t)
    for (t in seq_len(n_t)) {
      i <- by_time[[t]]
      if (length(i)) mean_dens[t] <- sum(mu[i] * mu[i] / sum_dens[t])
    }
    eao <- log_eao <- rep(0, n_t)
    keep <- include & lengths(by_time) > 0L
    eao[keep] <- out$total[keep] / mean_dens[keep]
    log_eao[keep] <- log(eao[keep])
    out$mean_dens <- mean_dens
    out$eao <- eps_values$eao <- eao
    out$log_eao <- log_eao
  }
  if (length(par$eps_index)) {
    for (values in eps_values) {
      out$nll <- out$nll + sum(par$eps_index[include] * values[include])
    }
  }
  out
}
