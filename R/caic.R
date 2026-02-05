#' Calculate conditional AIC
#'
#' Calculates the conditional Akaike Information criterion (cAIC).
#'
#' @param object Output from [sdmTMB()].
#' @param what Whether to return the cAIC or the effective degrees of freedom
#'   (EDF) for each group of random effects.
#' @param ... Other arguments for specific methods. Not used.
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
#' @return
#' Either the cAIC or the effective degrees of freedom (EDF) by group
#' of random effects depending on the argument `what`.
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
#' @export
cAIC <- function(object, what = c("cAIC", "EDF"), ...) {
  UseMethod("cAIC", object)
}

#' @exportS3Method
cAIC.sdmTMB <- function(object, what = c("cAIC", "EDF"), ...) {

  what <- tolower(what)
  what <- match.arg(what, choices = c("caic", "edf"))

  if ("edf" %in% names(object) && what == "edf") {
    return(object$edf)
  }

  tmb_data <- object$tmb_data

  ## Ensure profile = NULL
  if (is.null(object$control$profile)) {
    obj <- object$tmb_obj
  } else {
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = object$parlist,
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      profile = NULL #<
    )
  }

  ## Make obj_new
  tmb_data$weights_i[] <- 0
  obj_new <- TMB::MakeADFun(
    data = tmb_data,
    parameters = object$parlist,
    map = object$tmb_map,
    random = object$tmb_random,
    DLL = "sdmTMB",
    profile = NULL
  )

  par <- obj$env$parList()
  parDataMode <- obj$env$last.par
  indx <- obj$env$lrandom()
  q <- sum(indx)
  p <- length(object$model$par)

  ## use '-' for Hess because model returns negative loglikelihood
  if (is.null(object$tmb_random)) {
    cli_inform(c("This model has no random effects.", "cAIC and EDF only apply to models with random effects."))
    return(invisible(NULL))
  }
  Hess_new <- -Matrix::Matrix(obj_new$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess_new <- Hess_new[indx, indx] ## marginal precision matrix of REs

  ## Joint hessian etc
  Hess <- -Matrix::Matrix(obj$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess <- Hess[indx, indx]
  negEDF <- Matrix::diag(Matrix::solve(Hess, Hess_new, sparse = FALSE))

  if (what == "caic") {
    jnll <- obj$env$f(parDataMode)
    cnll <- jnll - obj_new$env$f(parDataMode)
    cAIC_out <- 2 * cnll + 2 * (p + q) - 2 * sum(negEDF)
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
    if (.is_regular_delta_model(object)) s_groups <- c(paste0("1LP-", s_groups), paste0("2LP-", s_groups))
    group[group == "b_smooth"] <- s_groups
    group <- factor(group)

    ## Calculate total EDF by group
    EDF <- tapply(negEDF, INDEX = group, FUN = length) - tapply(negEDF, INDEX = group, FUN = sum)
    return(EDF)
  } else {
    cli_abort("Option not implemented")
  }
}
