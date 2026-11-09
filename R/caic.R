#' Calculate conditional AIC
#'
#' Calculates the conditional Akaike Information criterion (cAIC).
#'
#' @param object Output from [sdmTMB()] or [profile_breakpt()].
#' @param ... Optionally, more fitted models from [sdmTMB()] or
#'   [profile_breakpt()] to compare, as with [AIC()].
#' @param what Whether to return the cAIC or the effective degrees of freedom
#'   (EDF) for each group of random effects. Must be named. `"EDF"` is only
#'   available for a single model.
#'
#' @details cAIC is designed to optimize the expected out-of-sample predictive
#' performance for new data that share the same random effects as the in-sample
#' (fitted) data, e.g., spatial interpolation.  In this sense, it should be a
#' fast approximation to optimizing the model structure based on k-fold
#' cross-validation.
#'
#' By contrast, [AIC()] calculates the marginal Akaike Information Criterion,
#' which is designed to optimize expected predictive performance for new data
#' that have new random effects, e.g., extrapolation, or inference about
#' generative parameters.
#'
#' cAIC also calculates the effective degrees of freedom (EDF) as a byproduct.
#' This is the number of fixed effects that would have an equivalent impact on
#' model flexibility as a given random effect.
#'
#' Both cAIC and EDF are calculated using Eq. 6 of Zheng, Cadigan, and Thorson
#' (2024).
#'
#' For models that include profiled fixed effects, these profiles are turned
#' off.
#'
#' When comparing several models, models without random effects contribute
#' their marginal AIC, which is what cAIC reduces to with no random effects.
#'
#' For output `p` from [profile_breakpt()], cAIC uses the final fit and counts
#' the profiled cutpoint as one more fixed effect, as `AIC(p)` does. Use
#' `cAIC(p)`, not `cAIC(p$fit)`: `p$fit` holds the cutpoint fixed, so it isn't
#' counted. The EDF (`what = "EDF"`) are those of `p$fit`.
#'
#' @return
#' For a single model, either the cAIC or the effective degrees of freedom
#' (EDF) by group of random effects depending on the argument `what`. For
#' several models, a data frame with a row per model and columns `df` (the
#' effective degrees of freedom: fixed effects plus total EDF) and `cAIC`, as
#' with [AIC()].
#'
#' @references
#' **Deriving the general approximation to cAIC used here:**
#'
#' Zheng, N., Cadigan, N., & Thorson, J. T. (2024).
#' A note on numerical evaluation of conditional Akaike information for
#' nonlinear mixed-effects models (arXiv:2411.14185). arXiv.
#' \doi{10.48550/arXiv.2411.14185}
#'
#' **The utility of EDF to diagnose hierarchical model behaviour:**
#'
#' Thorson, J. T. (2024). Measuring complexity for hierarchical
#' models using effective degrees of freedom. Ecology,
#' 105(7), e4327 \doi{10.1002/ecy.4327}
#'
#' @examples
#' mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 15)
#' fit <- sdmTMB(catch_weight ~ s(log(depth)),
#'   time_varying = ~1,
#'   time_varying_type = "ar1",
#'   time = "year",
#'   spatiotemporal = "off",
#'   mesh = mesh,
#'   family = tweedie(),
#'   data = dogfish,
#'   offset = log(dogfish$area_swept)
#' )
#' cAIC(fit)
#' cAIC(fit, what = "EDF")
#' AIC(fit)
#'
#' fit2 <- update(fit, catch_weight ~ 1)
#' cAIC(fit, fit2)
#' @export
cAIC <- function(object, ..., what = c("cAIC", "EDF")) {
  UseMethod("cAIC", object)
}

#' @exportS3Method
cAIC.sdmTMB <- function(object, ..., what = c("cAIC", "EDF")) {
  dots <- list(...)
  ## Short-lived check: `what` used to be the second argument
  ## TODO: remove in a future release
  if (length(dots) && is.character(dots[[1L]])) {
    cli_abort(c(
      "`what` must now be named.",
      i = "For example, use `cAIC(fit, what = \"{dots[[1L]][1L]}\")`."
    ))
  }
  what <- tolower(what)
  what <- match.arg(what, choices = c("caic", "edf"))

  if (!length(dots)) {
    return(caic_one(object, what = what))
  }

  if (what == "edf") {
    cli_abort("`what = \"EDF\"` is only available for a single model.")
  }
  objects <- c(list(object), dots)
  is_fit <- vapply(objects, inherits, logical(1L), what = c("sdmTMB", "breakpt_profile"))
  if (!all(is_fit)) {
    cli_abort("All objects passed to `cAIC()` must be fitted `sdmTMB()` models or `profile_breakpt()` output.")
  }
  nobs_all <- vapply(objects, stats::nobs, numeric(1L))
  if (length(unique(nobs_all)) > 1L) {
    cli_warn("Models are not all fitted to the same number of observations.")
  }
  out <- lapply(objects, caic_one, what = "caic", table = TRUE)
  out <- data.frame(
    df = vapply(out, `[[`, numeric(1L), "df"),
    cAIC = vapply(out, `[[`, numeric(1L), "cAIC")
  )
  Call <- match.call()
  Call$what <- NULL
  row.names(out) <- vapply(as.list(Call)[-1L], deparse1, character(1L))
  out
}

# Same code: match.call() in a method sees the generic's call, so row names
# still come from the caller's arguments
#' @exportS3Method cAIC breakpt_profile
cAIC.breakpt_profile <- cAIC.sdmTMB

caic_one <- function(object, what = c("caic", "edf"), table = FALSE) {
  what <- match.arg(what)
  if (inherits(object, "breakpt_profile")) {
    if (is.null(object$fit)) {
      cli_abort("`cAIC()` needs the final fit: use `profile_breakpt(final = TRUE)`.")
    }
    out <- caic_one(object$fit, what = what, table = table)
    if (what == "edf" || is.null(out)) return(out)
    # the final fit holds the profiled cutpoint fixed; count it, as AIC() does
    if (table) return(list(cAIC = out$cAIC + 2, df = out$df + 1))
    return(out + 2)
  }
  check_fit_version(object)
  .check_family_capability(object, "cAIC")

  if ("edf" %in% names(object) && what == "edf") {
    return(object$edf)
  }

  if (is.null(object$tmb_random)) {
    if (table) {
      ## with no random effects, cAIC reduces to the marginal AIC
      return(list(cAIC = stats::AIC(object), df = attr(stats::logLik(object), "df")))
    }
    cli_inform(c("This model has no random effects.", "cAIC and EDF only apply to models with random effects."))
    return(invisible(NULL))
  }

  tmb_data <- object$tmb_data

  ## Ensure profile = NULL
  if (is.null(object$control$profile)) {
    obj <- object$tmb_obj
  } else {
    obj <- make_sdmTMB_adfun(
      data = tmb_data,
      parameters = object$parlist,
      map = object$tmb_map,
      random = object$tmb_random,
          profile = NULL #<
    )
  }

  ## Make obj_new
  tmb_data$weights_i[] <- 0
  obj_new <- make_sdmTMB_adfun(
    data = tmb_data,
    parameters = object$parlist,
    map = object$tmb_map,
    random = object$tmb_random,
      profile = NULL
  )

  # the mode, not `last.par`: that's wherever the objective was last evaluated
  parDataMode <- object$last.par.best
  q <- sum(obj$env$lrandom())
  # fixed effects of `obj`, which includes any profiled ones
  p <- length(obj$par)

  ## Sparse Hessians of the random effects:
  ## joint (Hess) and with data weights set to zero (Hess_new; prior precision of REs)
  Hess <- obj$env$spHess(parDataMode, random = TRUE)
  Hess_new <- obj_new$env$spHess(parDataMode, random = TRUE)

  ## diag(Hess^-1 %*% Hess_new) only needs elements of Hess^-1 on the sparsity
  ## pattern of Hess, which solveSubset() computes from the sparse Cholesky
  solveSubset <- utils::getFromNamespace("solveSubset", "TMB")
  negEDF <- Matrix::rowSums(solveSubset(Hess) * Hess_new)

  if (what == "caic") {
    jnll <- obj$env$f(parDataMode)
    cnll <- jnll - obj_new$env$f(parDataMode)
    df <- p + q - sum(negEDF)
    cAIC_out <- 2 * cnll + 2 * df
    if (table) return(list(cAIC = cAIC_out, df = df))
    return(cAIC_out)
  } else if (what == "edf") {
    ## Figure out group for each random-effect coefficient
    group <- names(object$last.par.best[obj$env$random])

    convert_bsmooth2names <- function(object, model = 1) {
      sn <- row.names(print_smooth_effects(object, m = model, silent = TRUE)$smooth_sds)
      sn <- gsub("^sd", "", sn)
      dms <- object$smoothers$sm_dims
      unlist(lapply(seq_along(dms), \(i) rep(sn[i], dms[i])))

    }
    s_groups <- convert_bsmooth2names(object)
    # smoothers always shared in delta models
    if (.object_has_two_components(object, caller = "`cAIC()`")) s_groups <- c(paste0("1LP-", s_groups), paste0("2LP-", s_groups))
    group[group == "b_smooth"] <- s_groups
    group <- factor(group)

    ## Calculate total EDF by group
    EDF <- tapply(negEDF, INDEX = group, FUN = length) - tapply(negEDF, INDEX = group, FUN = sum)
    return(EDF)
  } else {
    cli_abort("Option not implemented")
  }
}
