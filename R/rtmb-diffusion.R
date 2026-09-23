# Covariate diffusion for the RTMB objective, following
# `src/covariate-diffusion.h`. Each term transforms one covariate on mesh
# vertices by time:
#
# - spatial: (M0 + kappaS^-2 M1) z_t = M0 x_t
# - time: z_t = (x_t + kappaT z_{t-1}) / (1 + kappaT), with z_{-1} = 0
# - joint: ((1 + kappaT) M0 + kappaS^-2 M1) z_t = M0 x_t + kappaT M0 z_{t-1}
#
# The formula parser combines `diffusion(x)` and `time_lag(x)` into one joint
# term, so each covariate has at most one term.

# Diffusion terms by name ("space", "time", or "joint") and one-based
# covariate, with the SPDE matrices they solve against.
rtmb_diffusion_inputs <- function(data) {
  cd <- data$covariate_diffusion
  list(n_terms = cd$n_terms, n_covariates = cd$n_covariates,
    type = c("space", "time", "joint")[cd$term_component + 1L],
    covariate = cd$term_covariate + 1L,
    M0 = data$spde$M0, M1 = data$spde$M1)
}

# Diffused term values at fitted or projected rows (n by term).
rtmb_diffusion_terms <- function(par, diffusion, rows) {
  x <- rows$diffusion_x
  M0 <- diffusion$M0
  M1 <- diffusion$M1
  terms <- lapply(seq_len(diffusion$n_terms), function(term) {
    type <- diffusion$type[[term]]
    k <- diffusion$covariate[[term]]
    X <- matrix(x[, , k], dim(x)[1L], dim(x)[2L])
    kappaT <- par$kappaT_nl_raw[[k]]
    if (type == "time") {
      z <- vector("list", ncol(X))
      z[[1L]] <- X[, 1L] / (1 + kappaT)
      for (t in seq_len(ncol(X))[-1L]) {
        z[[t]] <- (X[, t] + kappaT * z[[t - 1L]]) / (1 + kappaT)
      }
    } else {
      system <- (if (type == "joint") 1 + kappaT else 1) * M0 +
        exp(-2 * par$log_kappaS_nl[[k]]) * M1
      rhs <- as.matrix(M0 %*% X)
      if (type == "space") {
        solved <- as.matrix(Matrix::solve(system, rhs))
        z <- lapply(seq_len(ncol(X)), function(t) solved[, t])
      } else {
        z <- vector("list", ncol(X))
        for (t in seq_len(ncol(X))) {
          b <- rhs[, t]
          if (t > 1L) b <- b + kappaT * rtmb_product(M0, z[[t - 1L]])
          z[[t]] <- rtmb_rhs_solve(system, b)
        }
      }
    }
    rtmb_station_values(rows, z)
  })
  do.call(cbind, terms)
}

rtmb_rhs_solve <- function(system, b) rtmb_as_vector(Matrix::solve(system, b))

# Derived diffusion scales, as reported by C++, for covariates with spatial
# or temporal terms: kappaS_nl, the mean squared displacement MSDK =
# 4 / kappaS^2 (times 1 - rhoT for joint terms) and its root and logs,
# kappaT_nl, and rhoT = kappaT / (1 + kappaT).
rtmb_diffusion_scales <- function(par, diffusion) {
  has <- function(types) {
    vapply(seq_len(diffusion$n_covariates), function(k) {
      any(diffusion$type[diffusion$covariate == k] %in% types)
    }, logical(1L))
  }
  space <- which(has(c("space", "joint")))
  time <- which(has(c("time", "joint")))
  out <- list()
  if (length(space)) {
    kappaS_nl <- exp(par$log_kappaS_nl[space])
    kappaT <- par$kappaT_nl_raw[space]
    MSDK <- 4 / kappaS_nl^2 *
      (1 - as.numeric(space %in% time) * kappaT / (1 + kappaT))
    out <- c(out, list(kappaS_nl = kappaS_nl, MSDK = MSDK,
      RMSDK = sqrt(MSDK), log_MSDK = log(MSDK), log_RMSDK = log(sqrt(MSDK))))
  }
  if (length(time)) {
    kappaT_nl <- par$kappaT_nl_raw[time]
    out <- c(out,
      list(kappaT_nl = kappaT_nl, rhoT = kappaT_nl / (1 + kappaT_nl)))
  }
  out
}
