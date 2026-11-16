#' Estimated marginal means with the \pkg{emmeans} package with \pkg{sdmTMB}
#'
#' @description
#' Methods for using the \pkg{emmeans} package with \pkg{sdmTMB}. The
#' \pkg{emmeans} package computes estimated marginal means for fixed effects.
#'
#' For delta/hurdle models, you can specify which component to analyze using the
#' `model` argument: `model = 1` for the binomial component (encounter
#' probability) or `model = 2` for the positive component (e.g., gamma for
#' `delta_gamma()`). By default, `model = 1`.
#'
#' Estimated marginal means are population-level: they include fixed effects
#' and smoothers but not random fields or random intercepts/slopes. They are
#' for an offset of 0 unless the \pkg{emmeans} `offset` argument is supplied.
#' Models with `breakpt()`, `logistic()`, `time_varying`, or covariate
#' diffusion terms aren't supported.
#'
#' @name emmeans.sdmTMB
#'
#' @references
#' \url{https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/}
#'
#' @examplesIf require("emmeans", quietly = TRUE)
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' fit <- sdmTMB(
#'   present ~ as.factor(year),
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit
#' emmeans::emmeans(fit, ~ year)
#' emmeans::emmeans(fit, pairwise ~ year)
#' emmeans::emmeans(fit, pairwise ~ year, type = "response")
#' emmeans::emmeans(fit, pairwise ~ year, adjust = "none")
#'
#' e <- emmeans::emmeans(fit, ~ year)
#' plot(e)
#'
#' e <- emmeans::emmeans(fit, pairwise ~ year)
#' confint(e)
#' summary(e, infer = TRUE)
#' as.data.frame(e)
#'
#' # interaction of factor with continuous predictor:
#' fit2 <- sdmTMB(
#'   present ~ depth_scaled * as.factor(year),
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit2
#' # slopes for each level:
#' emmeans::emtrends(fit2, ~ year, var = "depth_scaled")
#' # test difference in slopes:
#' emmeans::emtrends(fit2, pairwise ~ year, var = "depth_scaled")
#' emmeans::emmip(fit2, year ~ depth_scaled,
#'   at = list(depth_scaled = seq(-2.5, 2.5, length.out = 50)), CIs = TRUE)
#'
#' # delta/hurdle models:
#' fit_delta <- sdmTMB(
#'   density ~ as.factor(year),
#'   data = pcod_2011, spatial = "off",
#'   family = delta_gamma()
#' )
#' # binomial component (encounter probability):
#' emmeans::emmeans(fit_delta, ~ year, model = 1)
#' # positive component (gamma):
#' emmeans::emmeans(fit_delta, ~ year, model = 2)
#'
#' # smoothers:
#' fit_sm <- sdmTMB(
#'   present ~ as.factor(year) + s(depth_scaled),
#'   data = pcod_2011, spatial = "off",
#'   family = binomial()
#' )
#' emmeans::emmeans(fit_sm, ~ year, at = list(depth_scaled = 0))
#' emmeans::emmeans(fit_sm, ~ depth_scaled,
#'   at = list(depth_scaled = c(-1, 0, 1)))

NULL # don't document functions below

#' @exportS3Method NULL
recover_data.sdmTMB <- function(object, data = NULL, ...) {
  fcall <- stats::getCall(object)
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli_abort("Please install the emmeans package to use this function")
  }
  # use the stored data rather than re-evaluating `fcall$data`, which may not
  # resolve in the formula environment (e.g., a fit made inside a function)
  if (is.null(data)) data <- model.frame(object)
  # emmeans would otherwise look for the `offset` argument in the call;
  # EMMs are for an offset of 0 unless emmeans' own `offset` is supplied
  fcall$offset <- NULL
  emmeans::recover_data(
    fcall, emm_terms(object),
    attr(model.frame(object), "na.action"), data = data, ...
  )
}

# Terms spanning the predictors of all components, with each smoother's
# covariates (and any `by` variable) as plain terms so they enter the
# reference grid. Keeps the fitted `predvars` (e.g., `poly()` coefficients).
emm_terms <- function(object) {
  trms <- lapply(object$terms, stats::delete.response)
  labels <- unlist(lapply(trms, attr, which = "term.labels"))
  if (isTRUE(object$smoothers$has_smooths)) {
    sm_vars <- lapply(unlist(object$smoothers$basis_out, recursive = FALSE),
      function(b) c(b$term, if (b$by != "NA") b$by))
    labels <- c(labels, unlist(sm_vars))
  }
  labels <- unique(labels)
  if (!length(labels)) labels <- "1"
  f <- stats::reformulate(labels, env = environment(object$terms[[1]]))
  out <- stats::delete.response(stats::terms(f))
  fitted_vars <- unlist(lapply(trms, function(x) as.list(attr(x, "variables"))[-1]))
  fitted_predvars <- unlist(lapply(trms, function(x) {
    pv <- attr(x, "predvars")
    if (is.null(pv)) pv <- attr(x, "variables")
    as.list(pv)[-1]
  }))
  vars <- as.list(attr(out, "variables"))[-1]
  i <- match(vapply(vars, deparse1, ""), vapply(fitted_vars, deparse1, ""))
  vars[!is.na(i)] <- fitted_predvars[i[!is.na(i)]]
  attr(out, "predvars") <- as.call(c(quote(list), vars))
  out
}

# Model terms that aren't linear in fixed and smoother coefficients
emm_check_terms <- function(object) {
  bad <- c(
    "`breakpt()` or `logistic()`" = isTRUE(object$threshold_function > 0L),
    "`time_varying`" = !is.null(object$time_varying),
    "covariate diffusion" = !is.null(object$nonlocal_parsed)
  )
  if (any(bad)) {
    cli_abort(c(
      "`emmeans()` doesn't support models with {names(bad)[bad]} terms.",
      "i" = "Consider `predict()` on a grid of covariate values instead."
    ))
  }
}

# with help from emm_basis.glmmTMB
#' @exportS3Method NULL
emm_basis.sdmTMB <- function(object, trms, xlev, grid, ...) {
  .check_family_capability(object, "emmeans", caller = "`emmeans()`")
  emm_check_terms(object)
  dots <- list(...)
  model <- if ("model" %in% names(dots)) dots$model else 1L
  is_delta <- .object_is_delta(object, caller = "`emmeans()`")
  if (is_delta) {
    model <- as.integer(model)
    if (!model %in% c(1L, 2L)) {
      cli_abort("`model` must be 1 (binomial component) or 2 (positive component) for delta models.")
    }
    fam <- .object_family_spec(object, caller = "`emmeans()`")$family[[model]]
  } else {
    model <- 1L
    fam <- family(object)
  }
  misc <- emmeans::.std.link.labels(fam, list())

  # parametric part:
  m_trms <- stats::delete.response(object$terms[[model]])
  mf <- model.frame(m_trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- model.matrix(m_trms, mf, contrasts.arg = object$contrasts[[model]])
  b_j <- if (model == 2L) "b_j2" else "b_j"
  pos <- list(match(colnames(X), colnames(object$tmb_data$X_ij[[model]])))
  names(pos) <- b_j

  sm <- object$smoothers
  if (isTRUE(sm$has_smooths) && any(sm$block_active[, model])) {
    sm_new <- parse_smoothers(sm$formula_no_bars, data = object$data,
      newdata = grid, basis_prev = sm$basis_out)
    cols <- sm$col_active[, model]
    blocks <- sm$block_active[, model]
    Z <- do.call(cbind, sm_new$Zs[blocks])
    X_sm <- cbind(sm_new$Xs[, cols, drop = FALSE], Z)
    colnames(X_sm) <- paste0(".smooth", seq_len(ncol(X_sm)))
    X <- cbind(X, X_sm)
    pos$bs <- which(cols) + (model - 1L) * length(cols)
    coef_active <- rep(blocks, sm$sm_dims)
    pos$b_smooth <- which(coef_active) + (model - 1L) * length(coef_active)
  }

  est <- emm_coefs(object, pos)
  bhat <- stats::setNames(est$bhat, colnames(X))
  V <- est$V
  dimnames(V) <- list(colnames(X), colnames(X))

  dfargs <- list(df = df.residual(object))
  dffun <- function(k, dfargs) dfargs$df
  nbasis <- estimability::all.estble
  named_list(X, bhat, nbasis, V, dffun, dfargs, misc)
}

# Estimates and joint covariance of the elements `pos` (a list of positions
# within each named parameter array). Mapped-off elements get zero variance.
emm_coefs <- function(object, pos) {
  env <- object$tmb_obj$env
  par <- env$last.par.best
  par_list <- env$parList(par = par)
  bhat <- unlist(Map(function(p, i) as.vector(par_list[[p]])[i], names(pos), pos))
  idx <- unlist(Map(function(p, i) {
    map <- object$tmb_map[[p]]
    k <- if (is.null(map)) i else as.integer(map)[i]
    match(p, names(par)) - 1L + k
  }, names(pos), pos))
  ok <- !is.na(idx)
  V <- matrix(0, length(idx), length(idx))
  fixed <- match(idx, setdiff(seq_along(par), env$random))
  if (all(!is.na(fixed[ok]))) {
    V[ok, ok] <- object$sd_report$cov.fixed[fixed[ok], fixed[ok], drop = FALSE]
  } else {
    Q <- sd_report_joint(object)$jointPrecision
    E <- Matrix::sparseMatrix(i = idx[ok], j = seq_len(sum(ok)),
      x = 1, dims = c(nrow(Q), sum(ok)))
    V[ok, ok] <- as.matrix(Matrix::solve(Q, E))[idx[ok], , drop = FALSE]
  }
  list(bhat = unname(bhat), V = V)
}
