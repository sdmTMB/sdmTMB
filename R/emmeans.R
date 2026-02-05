#' Estimated marginal means with the \pkg{emmeans} package with \pkg{sdmTMB}
#'
#' @description
#' Methods for using the \pkg{emmeans} package with \pkg{sdmTMB}. The
#' \pkg{emmeans} package computes estimated marginal means for the fixed
#' effects.
#'
#' For delta/hurdle models, you can specify which component to analyze using the
#' `model` argument: `model = 1` for the binomial component (encounter
#' probability) or `model = 2` for the positive component (e.g., gamma for
#' `delta_gamma()`). By default, `model = 1`.
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

NULL # don't document functions below

recover_data.sdmTMB <- function(object, ...) {
  fcall <- stats::getCall(object)
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli_abort("Please install the emmeans package to use this function")
  }
  emmeans::recover_data(
    fcall,
    stats::delete.response(terms(object)),
    attr(model.frame(object), "na.action"), ...
  )
}

# with help from emm_basis.glmmTMB
emm_basis.sdmTMB <- function(object, trms, xlev, grid, ...) {
  # Extract model argument if present (defaults to 1)
  dots <- list(...)
  model <- if ("model" %in% names(dots)) dots$model else 1L
  regular_delta <- .is_regular_delta_model(object)

  # For delta models, we need to work with one component at a time
  if (regular_delta) {
    model <- as.integer(model)
    if (!model %in% c(1L, 2L)) {
      cli_abort("`model` must be 1 (binomial component) or 2 (positive component) for delta models.")
    }
    # Get the appropriate formula for the delta component
    if (!identical(object$formula[[1]], object$formula[[2]])) {
      trms <- terms(object$formula[[model]])
    }
  }

  V <- vcov(object, model = model)
  misc <- list()

  # Get the appropriate family for delta models
  if (regular_delta) {
    fam <- object$family[[model]]
  } else {
    fam <- family(object)
  }
  misc <- emmeans::.std.link.labels(fam, misc)

  contrasts <- attr(model.matrix(object), "contrasts")
  contrasts <- contrasts[names(contrasts) %in% all.vars(terms(object))]
  m <- model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- model.matrix(trms, m, contrasts.arg = contrasts)

  # Get coefficients using existing fixef method and filter to parametric terms only
  # This fixes the issue with smoothers by only including parametric terms
  if (regular_delta) {
    all_bhat <- fixef(object, model = model)
  } else {
    all_bhat <- fixef(object)
  }

  # Only keep coefficients that correspond to columns in the design matrix X
  bhat <- all_bhat[names(all_bhat) %in% colnames(X)]
  # Ensure coefficient order matches design matrix column order
  bhat <- bhat[colnames(X)]

  if (length(bhat) < ncol(X)) {
    kept <- match(names(bhat), dimnames(X)[[2]])
    full_bhat <- NA * X[1, ]
    full_bhat[kept] <- bhat
    bhat <- full_bhat
    modmat <- model.matrix(
      trms, model.frame(object),
      contrasts.arg = contrasts
    )
    nbasis <- estimability::nonest.basis(modmat)
  } else {
    nbasis <- estimability::all.estble
  }
  dfargs <- list(df = df.residual(object))
  dffun <- function(k, dfargs) dfargs$df
  named_list(X, bhat, nbasis, V, dffun, dfargs, misc)
}
