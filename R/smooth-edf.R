# Effective degrees of freedom (EDF) of penalized smoothers, matching mgcv.
#
# The penalized coefficients of each smoother block k are `b_smooth ~ N(0,
# sigma_k^2)`, so the penalty Hessian is diag(1 / sigma_k^2) on those rows. The
# mgcv EDF of a coefficient is 1 - [H^-1]_jj / sigma_k^2, where H is the Hessian
# of the joint negative log likelihood with respect to the random effects *and*
# the linear (fixed effect) coefficients, conditional on the variance
# parameters. Each unpenalized null-space column of a smoother adds exactly 1.
#
# H^-1 restricted to the random effects is obtained from the inner-problem
# Cholesky factor that TMB has already analyzed; the linear coefficients enter
# through a Schur complement (cross terms by differencing the joint gradient),
# so no new tape or fresh sparse factorization is needed.
#
# Returns a (smoother block x linear predictor) matrix aligned with
# `ln_smooth_sigma`, or NULL if the EDF can't be computed.
calc_smooth_edf <- function(object) {
  sm <- object$smoothers
  if (!isTRUE(sm$has_smooths)) return(NULL)
  if (is.null(sm$zs_basis)) return(NULL) # fits from before this was stored
  tryCatch(smooth_edf_impl(object, sm), error = function(e) NULL)
}

smooth_edf_impl <- function(object, sm) {
  env <- object$tmb_obj$env
  par <- env$last.par.best
  nm <- names(par)
  r <- env$random
  s_idx <- which(nm[r] == "b_smooth")
  if (!length(s_idx)) return(NULL)

  # which (block, linear predictor) each estimated b_smooth element belongs to
  dims <- sm$sm_dims
  lss <- env$parList(par = par)$ln_smooth_sigma
  n_m <- ncol(lss)
  block <- rep(rep(seq_along(dims), dims), n_m)
  lp <- rep(seq_len(n_m), each = sum(dims))
  map <- object$tmb_map$b_smooth
  if (!is.null(map)) {
    keep <- !is.na(map)
    block <- block[keep]
    lp <- lp[keep]
  }
  if (length(block) != length(s_idx)) return(NULL)

  A <- env$spHess(par, random = TRUE)
  L <- env$L.created.by.newton
  L <- if (is.null(L)) {
    Matrix::Cholesky(A, perm = TRUE, LDL = FALSE)
  } else {
    Matrix::update(L, A)
  }
  E <- Matrix::sparseMatrix(
    i = s_idx, j = seq_along(s_idx), x = 1,
    dims = c(nrow(A), length(s_idx))
  )
  Ainv_s <- Matrix::solve(L, E, system = "A")
  hdiag <- as.numeric(Ainv_s[cbind(s_idx, seq_along(s_idx))])

  # linear coefficients that aren't already integrated (e.g., with `reml = TRUE`)
  f <- setdiff(which(nm %in% c("b_j", "b_j2", "bs")), r)
  if (length(f)) {
    Hf <- linear_re_hessian(env, par, r, f)
    B <- Hf[r, , drop = FALSE]
    C <- Hf[f, , drop = FALSE]
    C <- (C + t(C)) / 2
    W <- as.matrix(Matrix::solve(L, B, system = "A"))
    S <- C - crossprod(B, W)
    Ws <- W[s_idx, , drop = FALSE]
    hdiag <- hdiag + rowSums((Ws %*% solve(S)) * Ws)
  }

  sig2 <- exp(2 * lss[cbind(block, lp)])
  pen <- tapply(1 - hdiag / sig2, list(factor(block, seq_along(dims)), factor(lp, seq_len(n_m))), sum)
  pen <- matrix(pen, nrow = length(dims), ncol = n_m)

  # add each basis's null space to its first penalized block so that blocks of
  # a multi-penalty smoother (e.g., t2()) sum to the mgcv term EDF
  null <- numeric(length(dims))
  first <- !duplicated(sm$zs_basis)
  null[first] <- sm$xf_dims[sm$zs_basis[first]]
  pen + null
}

# Columns `f` of the Hessian of the joint negative log likelihood (all
# parameters), by central differences of the joint gradient. TMB's sparse
# Hessian only covers the random effects, and the joint NLL is quadratic in
# the linear coefficients for Gaussian models, so this is cheap and accurate.
linear_re_hessian <- function(env, par, r, f) {
  grad <- function(p) as.numeric(env$f(p, order = 1))
  H <- vapply(f, function(k) {
    h <- 1e-4 * max(1, abs(par[k]))
    up <- dn <- par
    up[k] <- up[k] + h
    dn[k] <- dn[k] - h
    (grad(up) - grad(dn)) / (2 * h)
  }, numeric(length(par)))
  env$f(par, order = 0) # leave the tape at the mode
  H
}
