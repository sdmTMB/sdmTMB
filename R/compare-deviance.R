#' Compare deviance among models
#'
#' Calculates the proportion of deviance explained by one or more models
#' relative to a reduced model (a pseudo-\eqn{R^2}), and the gain from each
#' added term or component, such as a spatial random field.
#'
#' @param object Output from [sdmTMB()]. The most complex model. Its shape
#'   parameters (see Details) are used for all models.
#' @param ... Other [sdmTMB()] models fit to the same response with the same
#'   family, ordered from most to least complex. The last is the reduced model,
#'   usually intercept-only with no random effects or random fields. Name the
#'   arguments to label the output rows. Names are arbitrary and only to
#'   help reading the output. These models can be set up with
#'   `do_fit = FALSE` and are then fit here. Alternatively, they can be fit
#'   first, in which case they are refit (if needed, see below) starting
#'   from the previously estimated parameter values for speed.
#'
#' @details
#' Deviance explained, \eqn{1 - D / D_0}, is a common pseudo-\eqn{R^2} for
#' generalized linear models (Cameron and Windmeijer 1997) and is what
#' \pkg{mgcv} reports. It reduces to \eqn{R^2} for a Gaussian linear model.
#'
#' **Sequences of models.** `gain` is the increase in deviance explained over
#' the next simpler model. Like sequential (type I) sums of squares, it depends
#' on order: a spatial field added after covariates usually gains less than one
#' added to an intercept-only model, since both can explain the same spatial
#' pattern.
#'
#' **Shape parameters.** For some families, the deviance depends on parameters
#' other than the mean: `phi` for negative binomial, [compois()], and
#' [lognormal()]; `p` for [tweedie()]; and `phi` and `Q` for [gengamma()].
#' Deviances computed with different values of these are not comparable, so
#' the other models are refit with them fixed at their values in `object`, as
#' \pkg{mgcv} does for its null deviance. Previously fit models are refit
#' starting from their estimates for speed. Reordering the models changes the
#' answer because it changes the comparisons.
#'
#' **Random effects and random fields.** The deviance is conditional on the
#' predicted random effects, so deviance explained measures in-sample fit of
#' the fixed and random effects together. Because it does not penalize
#' complexity, more flexible fields almost always increase it. Do not use it to
#' choose among random effect structures or meshes; use [cAIC()], [AIC()], or
#' [sdmTMB_cv()] instead. For families with shape parameters, fixing `phi` can
#' shift variation between random effects and observation error in refit
#' models, so gains from random effects are approximate.
#'
#' For delta models, the deviance is summed over both components.
#'
#' @return A data frame with one row per model with columns `model`,
#'   `deviance`, `deviance_explained` (relative to the last model), and `gain`
#'   (over the next simpler model). Any fixed shape parameter values (on their
#'   internal scale) are in attribute `"fixed"`.
#'
#' @references
#' Cameron, A. C., and Windmeijer, F. A. G. (1997). An R-squared measure of
#' goodness of fit for some common nonlinear regression models. Journal of
#' Econometrics, 77(2), 329--342. \doi{10.1016/S0304-4076(96)01818-0}
#'
#' @seealso [stats::deviance()], [cAIC()]
#'
#' @examples
#' fit <- sdmTMB(density ~ depth_scaled + depth_scaled2,
#'   data = pcod_2011, spatial = "off", family = tweedie())
#' fit0 <- sdmTMB(density ~ 1,
#'   data = pcod_2011, spatial = "off", family = tweedie())
#' compare_deviance(fit, fit0)
#'
#' # Not comparable, since each model estimates its own Tweedie p:
#' 1 - deviance(fit) / deviance(fit0)
#'
#' \donttest{
#' # Gain from adding a spatial random field after depth:
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' fit_sp <- update(fit, spatial = "on", mesh = mesh)
#' compare_deviance(fit_sp, depth = fit, intercept = fit0)
#'
#' # Models after the first can be left unfitted
#' # with `do_fit = FALSE`, since they are fit here with
#' # Tweedie p fixed anyway:
#' fit_sp_only <- update(fit0, spatial = "on", mesh = mesh, do_fit = FALSE)
#' compare_deviance(fit_sp, spatial = fit_sp_only, intercept = fit0)
#' }
#' @export
compare_deviance <- function(object, ...) {
  models <- c(list(object), list(...))
  labels <- vapply(as.list(substitute(list(object, ...)))[-1L], deparse1, "")
  if (!is.null(names(models))) {
    named <- names(models) != ""
    labels[named] <- names(models)[named]
  }
  if (length(models) < 2L) {
    cli_abort("Supply at least two models, ending with the reduced model.")
  }
  lapply(models, function(m) assert_that(inherits(m, "sdmTMB")))
  unfitted <- vapply(models, function(m) isFALSE(m$call$do_fit), logical(1L))
  if (unfitted[1L]) cli_abort("`object` must be a fitted model.")
  for (m in models[-1L]) .check_deviance_comparable(object, m)
  deviance(object) # checks the family has a deviance

  fam <- .object_family_spec(object)$components$family_name
  pars <- unique(unlist(.deviance_shape_pars[fam]))
  fixed <- object$parlist[pars]
  if ("ln_phi" %in% pars &&
    any(vapply(models, function(m) isTRUE(m$has_dispformula), logical(1L)))) {
    cli_abort(c(
      "The deviance for this family depends on `phi`, which cannot be held
      fixed when it is modelled with `dispformula`.",
      "i" = "Refit all models with `dispformula = ~ 1`."
    ))
  }
  if (length(fixed)) {
    par_labels <- .deviance_par_labels[pars]
    cli_inform(paste("Refitting the other models with {par_labels} fixed at",
      "the value{?s} from `object`, since the deviance depends on {?it/them}."))
  }
  refit <- unfitted | (length(fixed) > 0L & seq_along(models) > 1L)
  models[refit] <- lapply(models[refit], .refit_for_deviance, fixed = fixed)

  dev <- vapply(models, deviance, numeric(1L))
  explained <- 1 - dev / dev[length(dev)]
  out <- data.frame(
    model = labels,
    deviance = dev,
    deviance_explained = explained,
    gain = c(explained[-length(dev)] - explained[-1L], NA),
    row.names = NULL
  )
  attr(out, "fixed") <- fixed
  out
}

# Fit `m` with the parameters in `fixed` held at those values. A fitted `m` is
# warm started from its estimates. The deviance doesn't need an sdreport.
.refit_for_deviance <- function(m, fixed) {
  control <- m$control
  if (!isFALSE(m$call$do_fit)) {
    control$start <- m$parlist
    control$multiphase <- FALSE
  }
  control$start[names(fixed)] <- fixed
  control$map[names(fixed)] <- lapply(fixed,
    function(x) factor(rep(NA, length(x))))
  control$getsd <- FALSE
  suppressMessages(update(m, control = control, do_fit = TRUE))
}

.check_deviance_comparable <- function(object, m) {
  if (!identical(.object_family_spec(object)$components$family_name,
    .object_family_spec(m)$components$family_name)) {
    cli_abort("All models must use the same family.")
  }
  if (!isTRUE(all.equal(object$response, m$response,
    check.attributes = FALSE))) {
    cli_abort("All models must be fit to the same response data.")
  }
  if (!isTRUE(all.equal(object$tmb_data$weights_i, m$tmb_data$weights_i))) {
    cli_abort("All models must use the same likelihood `weights`.")
  }
  if (!isTRUE(all.equal(object$tmb_data$upr, m$tmb_data$upr))) {
    cli_abort("All models must use the same `censored_upper` bounds.")
  }
  invisible(TRUE)
}

# Parameters other than the mean that each family's deviance depends on.
.deviance_shape_pars <- list(
  nbinom1 = "ln_phi",
  nbinom2 = "ln_phi",
  compois = "ln_phi",
  lognormal = "ln_phi",
  tweedie = "thetaf",
  gengamma = c("ln_phi", "gengamma_Q")
)
.deviance_par_labels <- c(ln_phi = "phi", thetaf = "Tweedie p",
  gengamma_Q = "Q")
