# Effective degrees of freedom (EDF) of penalized smoothers, matching mgcv.
#
# The penalized coefficients of each smoother block k are `b_smooth ~ N(0,
# sigma_k^2)`, so the penalty Hessian is diag(1 / sigma_k^2) on those rows. The
# mgcv EDF of a coefficient is 1 - [H^-1]_jj / sigma_k^2, where H is the Hessian
# of the joint negative log likelihood with respect to the random effects *and*
# the linear (fixed effect) coefficients, conditional on the variance
# parameters. Each unpenalized null-space column of a smoother adds exactly 1.
#
# H^-1 restricted to the random effects comes from the inner-problem Cholesky
# factor, numerically refactorized at the mode but reusing TMB's symbolic
# analysis. The linear coefficients enter through a Schur complement, with
# cross terms from differencing the joint gradient, so no new tape is needed.
#
# Returns a (smoother block x linear predictor) matrix aligned with
# `ln_smooth_sigma`, or NULL if the EDF can't be computed.
calc_smooth_edf <- function(object) {
  sm <- object$smoothers
  if (!isTRUE(sm$has_smooths)) return(NULL)
  if (is.null(sm$zs_basis)) return(NULL) # fits from before this was stored
  tryCatch(smooth_edf_impl(object, sm), error = function(e) {
    cli_inform(c(
      "!" = "Could not compute the smoother effective degrees of freedom (EDF).",
      "x" = "{conditionMessage(e)}"
    ))
    NULL
  })
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
  if (length(block) != length(s_idx)) stop("Unexpected `b_smooth` mapping.")

  A <- env$spHess(par, random = TRUE)
  L <- env$L.created.by.newton
  L <- if (is.null(L)) {
    Matrix::Cholesky(A, perm = TRUE, LDL = FALSE)
  } else {
    Matrix::update(L, A)
  }
  n <- nrow(A)

  # diag(A^-1) at the smoother coefficients; solve in chunks of columns so
  # memory stays bounded with many random effects (e.g., spatiotemporal fields)
  hdiag <- numeric(length(s_idx))
  for (j in col_chunks(length(s_idx), n)) {
    ij <- cbind(s_idx[j], seq_along(j))
    E <- matrix(0, n, length(j))
    E[ij] <- 1
    hdiag[j] <- as.matrix(Matrix::solve(L, E, system = "A"))[ij]
  }

  # linear coefficients that aren't already integrated (e.g., with `reml = TRUE`)
  lin <- c("b_j", "b_j2", "bs")
  f <- setdiff(which(nm %in% lin), r)
  if (length(f)) {
    scale <- linear_coef_scale(object, nm, lin)
    Hf <- linear_re_hessian(env, par, f, step = 1e-4 / scale[f])
    B <- Hf[r, , drop = FALSE]
    S <- Hf[f, , drop = FALSE]
    S <- (S + t(S)) / 2
    Ws <- matrix(0, length(s_idx), length(f))
    for (j in col_chunks(length(f), n)) { # Schur complement C - B' A^-1 B
      W <- as.matrix(Matrix::solve(L, B[, j, drop = FALSE], system = "A"))
      S[, j] <- S[, j] - crossprod(B, W)
      Ws[, j] <- W[s_idx, , drop = FALSE]
    }
    hdiag <- hdiag + rowSums((Ws %*% solve(S)) * Ws)
  }

  # sum over the coefficients of each (block, linear predictor); NA for a block
  # with all coefficients mapped off
  sig2 <- exp(2 * lss[cbind(block, lp)])
  pen <- tapply(
    1 - hdiag / sig2,
    list(factor(block, seq_along(dims)), factor(lp, seq_len(n_m))),
    sum
  )
  pen <- matrix(pen, nrow = length(dims), ncol = n_m)

  # add each basis's null space to its first penalized block so that blocks of
  # a multi-penalty smoother (e.g., t2()) sum to the mgcv term EDF
  null <- numeric(length(dims))
  first <- !duplicated(sm$zs_basis)
  null[first] <- sm$xf_dims[sm$zs_basis[first]]
  pen + null
}

# Split `n_col` columns into chunks of at most `max_cells / n_row` columns.
col_chunks <- function(n_col, n_row, max_cells = 2e6) {
  size <- max(1L, floor(max_cells / n_row))
  split(seq_len(n_col), ceiling(seq_len(n_col) / size))
}

# For each element of `par` that is a linear coefficient (`lin`), the largest
# absolute value of its design column(s): how far a unit change moves the linear
# predictor. NA for other parameters. Mirrored coefficients take the maximum
# over their columns.
linear_coef_scale <- function(object, nm, lin) {
  d <- object$tmb_data
  colmax <- function(X) {
    if (!length(X)) return(numeric(0))
    apply(abs(X), 2L, max)
  }
  cols <- list(
    b_j = colmax(d$X_ij[[1L]]),
    b_j2 = if (length(d$X_ij) > 1L) colmax(d$X_ij[[2L]]) else numeric(0),
    bs = rep(colmax(d$Xs), ncol(object$tmb_params$bs))
  )
  out <- rep(NA_real_, length(nm))
  for (p in intersect(lin, nm)) {
    x <- cols[[p]]
    map <- object$tmb_map[[p]]
    if (!is.null(map)) x <- as.numeric(tapply(x, factor(map), max))
    if (length(x) != sum(nm == p)) stop("Unexpected `", p, "` dimensions.")
    out[nm == p] <- x
  }
  pmax(out, sqrt(.Machine$double.eps))
}

# Columns `f` of the Hessian of the joint negative log likelihood (all
# parameters), by central differences of the joint gradient with steps `step`.
# TMB's sparse Hessian only covers the random effects. The steps are scaled to
# the design so each one moves the linear predictor by the same small amount.
linear_re_hessian <- function(env, par, f, step) {
  grad <- function(p) as.numeric(env$f(p, order = 1))
  H <- vapply(seq_along(f), function(i) {
    up <- dn <- par
    up[f[i]] <- up[f[i]] + step[i]
    dn[f[i]] <- dn[f[i]] - step[i]
    (grad(up) - grad(dn)) / (2 * step[i])
  }, numeric(length(par)))
  env$f(par, order = 0) # leave the tape at the mode
  H
}
