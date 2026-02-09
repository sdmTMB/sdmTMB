#' @import methods
#' @export
summary.sdmTMB <- function(object, ..., digits) {
  print(object, ...)
}

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

#' Extract the number of observations of an sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.sdmTMB <- function(object, ...) {
    sum(!is.na(object$data[all.vars(object$formula[[1]])[1]]))
}

#' Get fitted values from an sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats predict
#' @export
#' @noRd
fitted.sdmTMB <- function(object, ...) {

  if (!"offset" %in% names(object))
    cli_abort("It looks like this was fit with an older version of sdmTMB. Try sdmTMB:::update_version(fit).")
  if (.is_multi_family_model(object)) {
    pred <- predict(object, newdata = object$data, type = "response", offset = object$offset)
    return(pred$est)
  }
  if (.is_regular_delta_model(object)) {
    inv1 <- object$family[[1]]$linkinv
    p <- predict(object, type = "link", offset = object$offset)
    p1 <- inv1(p$est1)
    inv2 <- object$family[[2]]$linkinv
    p2 <- inv2(p$est2)
    p1 * p2
  } else {
    inv <- object$family$linkinv
    inv(predict(object, type = "link", offset = object$offset)$est)
  }
}

#' Get fixed-effect coefficients
#'
#' @param object The fitted sdmTMB model object
#' @param complete Currently ignored
#' @param model Linear predictor for delta models. Defaults to the first
#'   linear predictor.
#' @param ... Currently ignored
#' @importFrom stats coef
#' @export
coef.sdmTMB <- function(object, complete = FALSE, model = 1, ...) {
  regular_delta <- .is_regular_delta_model(object)
  if (regular_delta) {
    assert_that(length(model) == 1L)
    model <- as.integer(model)
    assert_that(model %in% c(1L, 2L))
    msg <- paste0("Returning coefficients from linear predictor ", model, " based on the `model` argument.")
    cli_inform(msg)
  }
  x <- tidy(object, model = model)
  out <- x$estimate
  names(out) <- x$term
  out
}

#' Get variance-covariance matrix
#'
#' @param object The fitted sdmTMB model object
#' @param complete Currently ignored
#' @param model Linear predictor for delta models. Defaults to the first
#'   linear predictor.
#' @param ... Currently ignored
#' @importFrom stats vcov
#' @export
#' @noRd
vcov.sdmTMB <- function(object, complete = FALSE, model = 1, ...) {
  regular_delta <- .is_regular_delta_model(object)
  if (regular_delta) {
    assert_that(length(model) == 1L)
    model <- as.integer(model)
    assert_that(model %in% c(1L, 2L))
  }

  sdr <- object$sd_report
  v <- sdr$cov.fixed
  fe <- tidy(object, model = model, silent = TRUE)$term
  nm <- colnames(v)

  # For delta models, identify which b_j to use
  if (regular_delta) {
    if (model == 1L) {
      i <- grepl("^b_j$", nm)
    } else {
      i <- grepl("^b_j2$", nm)
    }
  } else {
    i <- grepl("^b_j$", nm)
  }

  if (sum(i)) {
    if (sum(i) == length(fe)) { # should always be true
      nm[i] <- fe
    }
  }
  colnames(v) <- nm
  rownames(v) <- nm
  if (isTRUE(complete)) {
    return(v)
  } else {
    return(v[i,i,drop=FALSE])
  }
}

#' Get CIs
#'
#' @param object The fitted sdmTMB model object
#' @param parm Parameters to return CIs
#' @param level CI level
#' @param ... Ignored
#' @importFrom stats confint
#' @export
#' @noRd
confint.sdmTMB <- function(object, parm, level = 0.95, ...) {
  td <- tidy(object, conf.int = TRUE, conf.level = level)
  x <- matrix(nrow = nrow(td), ncol = 3L)
  x[,3L] <- td$estimate
  x[,2L] <- td$conf.high
  x[,1L] <- td$conf.low
  p <- ((1 - level) / 2) * 100
  pn <- paste(c(p, 100 - p), "%")
  colnames(x) <- c(pn, "Estimate")
  rownames(x) <- td$term
  x
}

#' Extract the log likelihood of a sdmTMB model
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.sdmTMB <- function(object, ...) {
  val <- -object$model$objective
  nobs <- nobs.sdmTMB(object)
  lpb <- names(object$tmb_obj$env$last.par.best)
  ran <- c("omega_s", "epsilon_st", "zeta_s", "b_rw_t", "epsilon_re", "RE", "b_smooth", "re_b_pars")
  df <- sum(!lpb %in% ran)
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a sdmTMB model
#'
#' @param fit The fitted sdmTMB model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.sdmTMB <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(c(edf, c(-2 * L + k * edf)))
}

#' @importFrom stats family
#' @export
family.sdmTMB <- function (object, ...) {
  if (.has_delta_attr(object)) {
    which_model <- attr(object, "delta_model_predict")
    if (is.na(which_model)) which_model <- 2L # combined; for link
    return(object$family[[which_model]])
  }
  if ("visreg_model" %in% names(object)) {
    return(object$family[[object$visreg_model]])
  } else {
    return(object$family)
  }
}

#' @importFrom nlme fixef
#' @method fixef sdmTMB
#' @export
fixef.sdmTMB <- function(object, model = 1, ...) {
  .t <- tidy(object, model = model, silent = TRUE)
  bhat <- .t$estimate
  names(bhat) <- .t$term
  bhat
}

##' @importFrom nlme ranef
#' @method ranef sdmTMB
#' @export
ranef.sdmTMB <- function(object, ...) {
  .t <- tidy(object, "ran_vals", conf.int = FALSE, silent = TRUE)
  model_list <- list()

  # For non-delta models, there is no model column
  if (!"model" %in% names(.t)) {
    .t$model <- 1
  }

  for (i in seq_len(max(.t$model))) { # loop through models
    .t_sub <- .t[which(.t$model == i), ]
    groups <- unique(.t_sub$group_name) # names of groups for this model
    group_list <- vector("list", length = length(groups)) # create empty named list
    names(group_list) <- groups
    for (j in 1:length(groups)) {
      sub <- .t_sub[which(.t_sub$group_name == groups[j]), ]
      level_ids <- unique(sub$level_ids)
      sub <- sub[, c("group_name", "term", "estimate")]
      if (nrow(sub) > 0) {
        # convert long to wide, storing just estimates
        split_data <- split(sub$estimate, sub$term)
        wide_df <- as.data.frame(split_data) # Convert to wide format
        names(wide_df) <- unique(sub$term) # rename, fix .X issue
        rownames(wide_df) <- level_ids # add rownames, like lmer does
        # Create a list with the dataframe as an element named 'Dog'
        group_list[[j]] <- wide_df
        # names(group_list[[j]]) <- sub$group_name[1]
      } # end if
    } # end for j
    model_list[[i]] <- group_list
  }
  model_list
}

#' @importFrom stats deviance
#' @method deviance sdmTMB
#' @importFrom stats residuals
#' @export
deviance.sdmTMB <- function(object, ...) {
  regular_delta <- .is_regular_delta_model(object)
  implemented <- c("poisson", "Gamma", "binomial",
    "gaussian", "lognormal", "tweedie", "nbinom1", "nbinom2")
  if (!regular_delta) {
    if (!object$family$family %in% implemented) {
      cli_abort("Deviance not implemented for the fitted family")
    }
  } else {
    if (!object$family$family[[2]] %in% implemented) {
      cli_abort("Deviance not implemented for the fitted family")
    }
  }
  if (regular_delta) {
    r1 <- residuals(object, type = "deviance", model = 1)
    r2 <- residuals(object, type = "deviance", model = 2)
    r <- sum(r1^2 + r2^2)
  } else {
    r <- residuals(object, type = "deviance")
    r <- sum(r^2)
  }
  r
}

#' @importFrom stats df.residual
#' @method df.residual sdmTMB
#' @export
df.residual.sdmTMB <- function(object, ...) {
  nobs(object) - length(object$model$par)
}

.has_delta_attr <- function(x) {
  "delta_model_predict" %in% names(attributes(x))
}

#' @export
formula.sdmTMB <- function (x, ...) {
  if (.has_delta_attr(x)) {
    which_model <- attr(x, "delta_model_predict")
    if (!identical(x$formula[[1]], x$formula[[2]]) && is.na(which_model)) {
      cli_abort("Delta component formulas are not the same but ggeffects::ggpredict() is trying to predict on the combined model. For now, predict on one or the other component, or keep the formulas the same, or write your own prediction and plot code.")
    }
    if (is.na(which_model)) which_model <- 1L # combined take 1!?
    return(x$formula[[which_model]])
  }
  if (length(x$formula) > 1L) {
    if ("visreg_model" %in% names(x)) {
      return(x$formula[[x$visreg_model]])
    } else {
      return(x$formula)
    }
  } else {
    return(x$formula[[1]])
  }
}

#' @importFrom stats terms
#' @export
terms.sdmTMB <- function(x, ...) {
  # DELTA FIXME: hardcoded to model 1!
  # Get the base terms object (without smoothers)
  class(x) <- "glm" # fake
  out <- stats::terms(x)
  out <- out[[1]]

  # If model has smoothers, add the underlying variables to term.labels
  # This ensures ggeffects can properly create prediction grids
  if (!is.null(x$smoothers) && isTRUE(x$smoothers$has_smooths)) {
    # Extract variable names from smoother labels
    # e.g., "s(depth)" -> "depth", "t2(x,y)" -> c("x", "y")
    smooth_vars <- character(0)
    for (label in x$smoothers$labels) {
      # Remove s( or t2( prefix and ) suffix
      # Handle both single variable smooths like s(x) and multi-variable like t2(x,y)
      inner <- gsub("^[st]2?\\((.*)\\)$", "\\1", label)
      # Split by comma for multi-variable smooths
      vars <- strsplit(inner, ",")[[1]]
      # Trim whitespace
      vars <- trimws(vars)
      smooth_vars <- c(smooth_vars, vars)
    }

    # Add unique smooth variables to term.labels
    smooth_vars <- unique(smooth_vars)
    existing_labels <- attr(out, "term.labels")
    # Only add variables that aren't already in term.labels
    new_labels <- setdiff(smooth_vars, existing_labels)
    attr(out, "term.labels") <- c(existing_labels, new_labels)
  }

  out
}

#' Calculate effects
#'
#' Used by effects package
#'
#' @inheritParams effects::Effect
#'
#' @importFrom stats formula poisson
#'
#' @return
#' Output from [effects::effect()]. Can then be plotted with with associated
#' `plot()` method.
#'
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(effects::Effect, sdmTMB)
#' } else {
#'   export(Effect.sdmTMB)
#' }
#' @examplesIf require("effects", quietly = TRUE)
#' fit <- sdmTMB(present ~ depth_scaled, data = pcod_2011, family = binomial(),
#'   spatial = "off")
#' effects::effect("depth_scaled", fit)
#' plot(effects::effect("depth_scaled", fit))
Effect.sdmTMB <- function(focal.predictors, mod, ...) {
  if (!requireNamespace("effects", quietly = TRUE)) {
    cli_abort("Please install the effects package")
  }

  if (.is_delta_like_model(mod)) {
    msg <- paste0("Effect() and ggeffects::ggeffect() do not yet work with ",
      "sdmTMB delta/hurdle models. Please use ggeffects::ggpredict() instead.")
    cli_abort(msg)
  }

  vc <- vcov(mod)
  b <- tidy(mod, silent = TRUE)

  dummyfuns <- list(
    variance = function(mu) mu,
    initialize = expression(mustart = y + 0.1),
    dev.resids = function(...) stats::poisson()$dev.res(...)
  )

  fam <- family(mod)

  # from glmmTMB:
  dummyfuns <- list(
    variance = function(mu) mu,
    initialize = expression(mustart <- y + 0.1),
    dev.resids = function(...) poisson()$dev.res(...)
  )
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
  }

  coefs <- b$estimate
  names(coefs) <- b$term
  args <- list(
    call = mod$call,
    coefficients = coefs,
    vcov = vc,
    family = fam,
    # formula = formula(mod) # includes random intercepts...
    formula = remove_s_and_t2(mod$split_formula[[1]]$form_no_bars)
  )

  e <- utils::getS3method("Effect", "default", envir = asNamespace("effects"))
  e(focal.predictors, mod, ..., sources = args)
}

# get_term_names <- function(model) {
#   .names <- gsub(" ", "", labels(terms(model)))
#   .names
# }

#' @export
model.frame.sdmTMB <- function(formula, ...) {
  as.data.frame(formula$data) # no tibbles!
}

#' Update an sdmTMB model
#'
#' This method updates an sdmTMB model with new arguments, automatically
#' handling the mesh object to avoid environment issues when loading
#' models from saved files.
#'
#' @param object An sdmTMB model object
#' @param formula. Optional updated formula
#' @param ... Other arguments to update in the model call
#' @param evaluate If `TRUE` (default), the updated call is evaluated;
#'   if `FALSE`, the call is returned unevaluated
#'
#' @return An updated sdmTMB model object (if `evaluate = TRUE`) or
#'   an unevaluated call (if `evaluate = FALSE`)
#'
#' @examples
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log"))
#' fit2 <- update(fit, family = delta_gamma())
#' @export
#' @importFrom stats update
update.sdmTMB <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  new_args <- list(...)

  # handle formula update if provided
  if (!missing(formula.)) {
    if (is.null(call$formula)) {
      call$formula <- formula.
    } else {
      call$formula <- stats::update.formula(call$formula, formula.)
    }
  }

  # update other arguments
  if (length(new_args) > 0) {
    for (arg_name in names(new_args)) {
      call[[arg_name]] <- new_args[[arg_name]]
    }
  }

  # if mesh is not provided in new arguments,
  # use the mesh from the fitted object to avoid environment issues
  if (!"mesh" %in% names(new_args)) {
    call$mesh <- object$spde
  }

  # if data is not provided, use the original data
  if (!"data" %in% names(new_args)) {
    call$data <- object$data
  }

  # evaluate the updated call if requested
  if (evaluate) {
    # create an environment with the necessary objects
    eval_env <- new.env(parent = parent.frame())

    # add the objects from the fitted model to the evaluation environment
    eval_env$mesh <- object$spde
    eval_env$data <- object$data

    # evaluate the call in this environment
    eval(call, envir = eval_env)
  } else {
    call
  }
}

#' Extract residual standard deviation or dispersion parameter
#'
#' @param object The fitted sdmTMB model object
#' @param ... Currently ignored
#' @importFrom stats sigma
#' @method sigma sdmTMB
#' @export
sigma.sdmTMB <- function(object, ...) {
  if (isTRUE(object$has_dispformula)) {
    cli_abort(c(
      "`sigma()` is not available when `dispformula` is used because dispersion varies by observation.",
      "i" = "Use `tidy(object, effects = 'dispersion')` to inspect dispersion coefficients, or `dharma_residuals()` for residual diagnostics."
    ))
  }

  regular_delta <- .is_regular_delta_model(object)

  # Get family
  fam <- if (regular_delta) {
    # For delta models, use the positive model
    object$family[[2]]
  } else {
    object$family
  }

  family_name <- fam$family

  # Gaussian/lognormal have dispersion/scale parameter
  if (family_name %in% c("gaussian", "lognormal")) {
    # phi = dispersion parameter from TMB report
    tmb_obj <- object$tmb_obj

    phi_idx <- which(names(tmb_obj$env$last.par.best) == "ln_phi")
    if (length(phi_idx) > 0) {
      ln_phi <- as.numeric(tmb_obj$env$last.par.best[phi_idx])
      return(exp(ln_phi))
    }
  }

  if (family_name == "Gamma") {
    # For Gamma, phi is shape, sigma = 1/sqrt(shape)
    tmb_obj <- object$tmb_obj
    phi_idx <- which(names(tmb_obj$env$last.par.best) == "ln_phi")
    if (length(phi_idx) > 0) {
      ln_phi <- as.numeric(tmb_obj$env$last.par.best[phi_idx])
      shape <- exp(ln_phi)
      return(1 / sqrt(shape))
    }
  }

  if (family_name %in% c("nbinom1", "nbinom2", "tweedie")) {
    # For negative binomial, return phi (dispersion parameter)
    tmb_obj <- object$tmb_obj
    phi_idx <- which(names(tmb_obj$env$last.par.best) == "ln_phi")
    if (length(phi_idx) > 0) {
      ln_phi <- as.numeric(tmb_obj$env$last.par.best[phi_idx])
      return(exp(ln_phi))
    }
  }

  # Poisson, binomial don't have dispersion, return 1 (as in glmmTMB)
  if (family_name %in% c("poisson", "binomial")) {
    return(1)
  }

  cli_warn("sigma() not implemented for family: {family_name}")
  return(NA_real_)
}
