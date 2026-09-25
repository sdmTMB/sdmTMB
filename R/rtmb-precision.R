# Field precision matrices for the RTMB objective. `rtmb_precision_inputs()`
# builds constant sparse templates once (type "none" without SPDE fields);
# `rtmb_precision()` combines them with AD coefficients. RTMB differentiates
# linear combinations of constant sparse matrices, but not products of AD
# sparse matrices, so every product is expanded into constant terms here.

rtmb_precision_inputs <- function(data) {
  general <- function(x) methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix")
  if (data$spatial_model == 1L) {
    W <- general(data$W_ss)
    # (I - rho W)'(I - rho W) = I - rho (W + W') + rho^2 W'W
    return(list(type = "sar", I = general(Matrix::Diagonal(nrow(W))),
      W_sum = general(W + Matrix::t(W)), WtW = general(Matrix::crossprod(W))))
  }
  if (data$spatial_model == 2L) {
    W <- general(data$W_ss)
    # Isolated areas get a unit diagonal.
    D <- pmax(Matrix::rowSums(W), 1)
    return(list(type = "car", D = general(Matrix::Diagonal(x = D)), W = W))
  }
  if (data$no_spatial == 1L) return(list(type = "none"))
  if (data$barrier == 1L) {
    b <- data$spde_barrier
    fraction <- if (length(data$barrier_scaling) > 1L)
      data$barrier_scaling[[2L]] else 0.1
    C <- b$C0 + b$C1 * fraction^2
    D <- b$D0 + b$D1 * fraction^2
    iC <- Matrix::Diagonal(x = 1 / C)
    ICD <- Matrix::t(b$I) %*% iC %*% D
    return(list(type = "barrier",
      ICI = general(Matrix::t(b$I) %*% iC %*% b$I),
      ICD_DCI = general(ICD + Matrix::t(ICD)),
      DCD = general(Matrix::t(D) %*% iC %*% D)))
  }
  if (data$anisotropy == 1L) {
    a <- data$spde_aniso
    G1 <- rtmb_aniso_bases(a)
    G2 <- lapply(G1, function(left) lapply(G1, function(right) {
      general(left %*% a$G0_inv %*% right)
    }))
    return(list(type = "anisotropic", G0 = general(a$G0), G1 = G1, G2 = G2))
  }
  list(type = "spde", M0 = general(data$spde$M0),
    M1 = general(data$spde$M1), M2 = general(data$spde$M2))
}

# The anisotropic G1 matrix is linear in the entries of adj(H):
# G1 = H[2, 2] Gxx - H[1, 2] Gxy + H[1, 1] Gyy. Each basis sums
# E_a[, u] E_b[, v] / (4 area) over triangle edge pairs, as in TMB's
# `R_inla::Q_spde()` for anisotropic meshes.
rtmb_aniso_bases <- function(a) {
  edges <- list(a$E0, a$E1, a$E2)
  vertices <- a$TV + 1L
  pairs <- expand.grid(a = 1:3, b = 1:3)
  basis <- function(term) {
    x <- unlist(lapply(seq_len(nrow(pairs)), function(p) {
      term(edges[[pairs$a[[p]]]], edges[[pairs$b[[p]]]]) / (4 * a$Tri_Area)
    }))
    Matrix::sparseMatrix(i = as.vector(vertices[, pairs$a]),
      j = as.vector(vertices[, pairs$b]), x = x, dims = c(a$n_s, a$n_s))
  }
  list(
    xx = basis(function(e1, e2) e1[, 1L] * e2[, 1L]),
    xy = basis(function(e1, e2) e1[, 1L] * e2[, 2L] + e1[, 2L] * e2[, 1L]),
    yy = basis(function(e1, e2) e1[, 2L] * e2[, 2L])
  )
}

# Precision for component `m`; `r = 1` for spatial and SVC fields and `r = 2`
# for spatiotemporal fields. Areal precisions have no range parameter.
rtmb_precision <- function(inputs, theta, r, m) {
  switch(inputs$type,
    spde = {
      kappa2 <- theta$kappa[r, m]^2
      kappa2^2 * inputs$M0 + 2 * kappa2 * inputs$M1 + inputs$M2
    },
    anisotropic = {
      H <- theta$H[[m]]
      kappa2 <- theta$kappa[r, m]^2
      # adj(H) entries, in the order of the `rtmb_aniso_bases()` bases.
      coefficients <- list(H[2L, 2L], -H[1L, 2L], H[1L, 1L])
      Q <- kappa2^2 * inputs$G0
      for (a in 1:3) {
        Q <- Q + 2 * kappa2 * coefficients[[a]] * inputs$G1[[a]]
        for (b in 1:3) {
          Q <- Q + coefficients[[a]] * coefficients[[b]] * inputs$G2[[a]][[b]]
        }
      }
      Q
    },
    barrier = {
      # INLAspacetime barrier precision with unit marginal SD.
      range2 <- theta$range[r, m]^2
      (2 / (pi * range2)) * inputs$ICI + (2 / (8 * pi)) * inputs$ICD_DCI +
        (range2 * 2 / (64 * pi)) * inputs$DCD
    },
    sar = {
      rho <- theta$rho_sar[[m]]
      inputs$I - rho * inputs$W_sum + rho^2 * inputs$WtW
    },
    car = inputs$D - theta$alpha_car[[m]] * inputs$W
  )
}

# Log SD parameter of a field. For SPDE-based precisions (SPDE, anisotropic,
# barrier) this is the Matérn marginal SD, 1 / (tau kappa sqrt(4 pi)). For
# areal precisions it is the scale 1 / tau, which need not equal each area's
# marginal SD.
rtmb_log_field_sd <- function(ln_tau, ln_kappa, inputs) {
  if (rtmb_areal(inputs)) -ln_tau else -ln_tau - ln_kappa - log(4 * pi) / 2
}

# `dgmrf()` scale for a field with log SD parameter `log_sd`. The barrier
# precision already has unit marginal SD and areal fields are scaled by
# 1 / tau directly, so both use the SD itself. SPDE and anisotropic
# precisions are unscaled by tau, so their scale is 1 / tau.
rtmb_gmrf_scale <- function(log_sd, ln_kappa, inputs) {
  if (inputs$type == "barrier" || rtmb_areal(inputs)) return(exp(log_sd))
  exp(log_sd + ln_kappa + log(4 * pi) / 2)
}

rtmb_areal <- function(inputs) inputs$type %in% c("sar", "car")

# Anisotropy matrix H from its two free parameters, as in C++ `MakeH()`.
rtmb_aniso_H <- function(x) {
  "[<-" <- RTMB::ADoverload("[<-")
  H <- matrix(0, 2L, 2L)
  H[1L, 1L] <- exp(x[[1L]])
  H[1L, 2L] <- H[2L, 1L] <- x[[2L]]
  H[2L, 2L] <- (1 + x[[2L]]^2) / exp(x[[1L]])
  H
}
