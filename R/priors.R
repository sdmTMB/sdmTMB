#' Prior distributions
#'
#' @description
#' Optional priors/penalties on model parameters. This results in penalized
#' likelihood within TMB or can be used as priors if the model is passed to
#' \pkg{tmbstan} (see the Bayesian vignette).
#'
#' **Note that Jacobian adjustments are only made if `bayesian = TRUE`** when the
#' [sdmTMB()] model is fit. In other words, if the final model will be fit with
#' \pkg{tmbstan} and priors are specified, then `bayesian` should be set to
#' `TRUE`. Otherwise, leave `bayesian = FALSE`.
#'
#' @details
#' Pass these objects to the `priors` argument in [sdmTMB()].
#'
#' @details
#' `normal()` and `halfnormal()` define normal and half-normal priors that, for
#' now, must have a location (mean) parameter of 0. `halfnormal()` is the
#' same as `normal()` but can be used to make the syntax clearer. It is intended
#' to be used for parameters that have support `> 0`.
#'
#' @details
#' See \url{https://arxiv.org/abs/1503.00256} for a description of the
#' PC prior for Gaussian random fields. Quoting the discussion (and substituting
#' the argument names in `pc_matern()`):
#' "In the simulation study we observe good coverage of the equal-tailed 95%
#' credible intervals when the prior satisfies `P(sigma > sigma_lt) = 0.05` and
#' `P(range < range_gt) = 0.05`, where `sigma_lt` is between 2.5 to 40 times
#' the true marginal standard deviation and `range_gt` is between 1/10 and 1/2.5
#' of the true range."
#'
#' @details
#' Keep in mind that the range is dependent on the units and scale of the
#' coordinate system. In practice, you may choose to try fitting the model
#' without a PC prior and then constraining the model from there. A better
#' option would be to simulate from a model with a given range and sigma to
#' choose reasonable values for the system or base the prior on knowledge from a
#' model fit to a similar system but with more spatial information in the data.
#'
#' @references
#' Fuglstad, G.-A., Simpson, D., Lindgren, F., and Rue, H. (2016) Constructing
#' Priors that Penalize the Complexity of Gaussian Random Fields.
#' arXiv:1503.00256
#'
#' Simpson, D., Rue, H., Martins, T., Riebler, A., and Sørbye, S. (2015)
#' Penalising model component complexity: A principled, practical approach to
#' constructing priors. arXiv:1403.4630
#'
#' @param matern_s A PC (Penalized Complexity) prior (`pc_matern()`) on the
#'   spatial random field Matérn parameters.
#' @param matern_st Same as `matern_s` but for the spatiotemporal random field.
#'   Note that you will likely want to set `share_range = FALSE` if you choose
#'   to set both a spatial and spatiotemporal Matérn PC prior since they both
#'   include a prior on the spatial range parameter. A shared range (see
#'   `share_range` and `range_groups` in [sdmTMB()]) gets the range part of
#'   the prior once, from the first field sharing it that is on and has a
#'   prior: spatial before spatiotemporal, and the first delta component
#'   before the second, then spatially varying coefficients (`matern_svc`).
#'   If fields sharing a range have priors with different range parts
#'   (`range_gt` or `range_prob`), a warning is issued, since only the first
#'   is applied; different sigma parts are fine. Priors are skipped for fields
#'   that are off, except that with only spatially varying coefficients
#'   (`spatial = "off"`), `matern_s` sets the prior on their shared range
#'   (unless `matern_svc` is set).
#' @param matern_svc Same as `matern_s` but for the spatially varying
#'   coefficient fields (`spatial_varying` in [sdmTMB()]). One prior applies
#'   to every coefficient field in every model component. The sigma part is on
#'   the SD of each coefficient field, which is on the scale of the linear
#'   predictor per unit of its covariate. A single `sigma_lt` is therefore a
#'   modeling assumption that the coefficient fields have similar SDs.
#'   Standardizing covariates puts them on similar scales, but in a delta
#'   model the components are on different link scales (e.g., log-odds of
#'   encounter versus log positive density), so standardization alone doesn't
#'   make one `sigma_lt` suitable for both components. The range part doesn't
#'   depend on covariate scaling but, as for the other fields, is in the units
#'   of the spatial coordinates. A range shared with the spatial or
#'   spatiotemporal field (the default, or see `range_groups` in [sdmTMB()])
#'   gets the range part from those fields first. When set, `matern_svc`
#'   rather than `matern_s` sets the range prior for spatially varying
#'   coefficients with `spatial = "off"`. Requires the RTMB backend.
#' @param phi A `halfnormal()` prior for the dispersion parameter in the
#'   observation distribution.
#' @param ar1_rho A `normal()` prior for the AR1 random field parameter. Note
#'   the parameter has support `-1 < ar1_rho < 1`.
#' @param tweedie_p A `normal()` prior for the Tweedie power parameter. Note the
#'   parameter has support `1 < tweedie_p < 2` so choose a mean appropriately.
#' @param b `normal()` priors for the main population-level 'beta' effects.
#' @param sigma_V `gamma_cv()` or `lognormal_prior()` priors for any
#'   time-varying parameter SDs. Supply a single prior to apply it to all
#'   time-varying coefficients or a vector with one element per coefficient
#'   (`NA` for no prior). In delta models, the priors apply to both components.
#' @param threshold_breakpt_slope A `normal()` prior for the slope of the
#'   linear (hockey stick) function.
#' @param threshold_breakpt_cut A `normal()` prior for the cutoff of the
#'   linear (hockey stick) function.
#' @param threshold_logistic_s50 A `normal()` prior for the parameter at which
#'   f(x) = 0.5.
#' @param threshold_logistic_s95 A `normal()` prior for the parameter at which
#'   f(x) = 0.95.
#' @param threshold_logistic_smax A `normal()` prior for the parameter at which
#'   f(x) is maximized.
#' @param custom Optional function of `(par, theta)` returning one or more log
#'   density contributions (RTMB backend only; see **Custom priors** below).
#' @param custom_log_jacobian Optional function of `(par, theta)` returning the
#'   log absolute Jacobian for `custom`. Only used if `bayesian = TRUE`. Requires
#'   `custom`.
#'
#' @section Custom priors:
#' With the RTMB backend (the default), `custom` adds arbitrary log densities
#' to the joint objective. The function is called as `custom(par, theta)`, where
#' `par` is the list of raw (internal) parameters and `theta` is a list of
#' natural-scale versions of them. See [get_prior_parameters()] for their
#' names, shapes, and labels. It must return a numeric scalar, vector, or array
#' of log densities (names are optional), which are summed and subtracted from the
#' negative log likelihood. Custom terms are evaluated inside the objective,
#' so they are included in gradients, Hessians, the Laplace approximation, and
#' standard errors, and may involve random effects. They supplement, rather than
#' replace, the built-in priors and random effect distributions. Use
#' [get_prior_densities()] to see the contributions at the estimated parameters.
#'
#' The function must:
#'
#' * use operations and densities that \pkg{RTMB} can differentiate, such as
#'   `RTMB::dnorm(x, mean, sd, log = TRUE)`; `stats::dnorm()` can't take
#'   AD values;
#' * return *log* densities (remember `log = TRUE`, which can't be checked);
#' * not branch on parameter values (e.g., `if (par$b_j[1] > 0)`);
#' * be deterministic and free of side effects;
#' * also work on ordinary numeric inputs.
#'
#' Fixed hyperparameters can be defined in the function's enclosing environment.
#' The functions are saved with the fitted model, so keep their environments
#' small (e.g., avoid defining them inside a function holding large data).
#'
#' **Jacobians.** With `bayesian = FALSE`, custom priors are penalties: the
#' optimum is the posterior mode on the scale the prior is written on, and no
#' Jacobian is needed. With `bayesian = TRUE` (e.g., for sampling with
#' \pkg{tmbstan}), the result of `custom_log_jacobian(par, theta)` is also
#' added. No Jacobian adjustment is ever inferred: without
#' `custom_log_jacobian`, a prior on a transformed parameter (e.g.,
#' `theta$phi`) is not a valid prior density on that parameter's natural scale
#' for sampling. For a prior on `phi = exp(ln_phi)`, the log Jacobian is
#' `par$ln_phi`.
#'
#' **Matérn range and SDs.** The Matérn `range` and the field SDs (`sigma_O`,
#' `sigma_E`) are more complex than most parameters because they depend on each
#' other: each SD is a function of both `ln_tau_*` and `ln_kappa`. A Jacobian
#' for one of them alone isn't well defined; it needs the joint transformation
#' from (`ln_kappa`, `ln_tau_*`) to (`range`, `sigma_*`), whose log Jacobian is
#' `log(range) + log(sigma_*)`. With a shared range (the default), add
#' `log(range)` once plus `log(sigma_*)` for each field. For priors on these
#' parameters we suggest the built-in PC priors, [pc_matern()] via the
#' `matern_s` and `matern_st` arguments, which apply the Jacobian for you if
#' `bayesian = TRUE`.
#'
#' **Limitations.** Custom terms are not evaluated in the first phase of
#' `sdmTMBcontrol(multiphase = TRUE)`, which only finds starting values with
#' random fields turned off. Custom priors do not define a sampler:
#' [simulate.sdmTMB()] and other simulation methods ignore any change that
#' custom terms make to the distribution of random effects. Held-out
#' log likelihoods in [sdmTMB_cv()] exclude all priors. As with the built-in
#' priors, custom terms are part of the objective, so they are included in
#' [stats::logLik()] and [stats::AIC()].
#'
#' @rdname priors
#'
#' @return
#' A named list with values for the specified priors.
#'
#' @export
sdmTMBpriors <- function(
  matern_s = pc_matern(range_gt = NA, sigma_lt = NA),
  matern_st = pc_matern(range_gt = NA, sigma_lt = NA),
  phi = halfnormal(NA, NA),
  ar1_rho = normal(NA, NA),
  tweedie_p = normal(NA, NA),
  b = normal(NA, NA),
  sigma_V = gamma_cv(NA, NA),
  threshold_breakpt_slope = normal(NA, NA),
  threshold_breakpt_cut = normal(NA, NA),
  threshold_logistic_s50 = normal(NA, NA),
  threshold_logistic_s95 = normal(NA, NA),
  threshold_logistic_smax = normal(NA, NA),
  matern_svc = pc_matern(range_gt = NA, sigma_lt = NA),
  custom = NULL,
  custom_log_jacobian = NULL
) {
  assert_that(attr(matern_s, "dist") == "pc_matern")
  assert_that(attr(matern_st, "dist") == "pc_matern")
  assert_that(attr(matern_svc, "dist") == "pc_matern")
  assert_that(attr(phi, "dist") == "normal")
  assert_that(attr(sigma_V, "dist") %in% c("gamma", "lognormal"))
  assert_that(attr(tweedie_p, "dist") == "normal")
  assert_that(attr(b, "dist") %in% c("normal", "mvnormal"))
  assert_that(attr(threshold_breakpt_slope, "dist") == "normal")
  assert_that(attr(threshold_breakpt_cut, "dist") == "normal")
  assert_that(attr(threshold_logistic_s50, "dist") == "normal")
  assert_that(attr(threshold_logistic_s95, "dist") == "normal")
  assert_that(attr(threshold_logistic_smax, "dist") == "normal")
  if (!is.null(custom) && !is.function(custom)) {
    cli_abort("`custom` must be a function of `(par, theta)` or `NULL`.")
  }
  if (!is.null(custom_log_jacobian)) {
    if (!is.function(custom_log_jacobian)) {
      cli_abort("`custom_log_jacobian` must be a function of `(par, theta)` or `NULL`.")
    }
    if (is.null(custom)) cli_abort("`custom_log_jacobian` requires `custom`.")
  }
  list(
    matern_s = matern_s,
    matern_st = matern_st,
    phi = phi,
    ar1_rho = ar1_rho,
    tweedie_p = tweedie_p,
    b = b,
    sigma_V = sigma_V,
    threshold_breakpt_slope = threshold_breakpt_slope,
    threshold_breakpt_cut = threshold_breakpt_cut,
    threshold_logistic_s50 = threshold_logistic_s50,
    threshold_logistic_s95 = threshold_logistic_s95,
    threshold_logistic_smax = threshold_logistic_smax,
    # last numeric prior: C++ reads the others by position
    matern_svc = matern_svc,
    custom = custom,
    custom_log_jacobian = custom_log_jacobian
  )
}

#' @param location Location parameter(s). Typically the mean.
#' @param scale Scale parameter. For `normal()`/`halfnormal()`: standard
#'   deviation(s). For `mvnormal()`: variance-covariance matrix.
#' @export
#' @rdname priors
#' @examples
#' normal(0, 1)
normal <- function(location = 0, scale = 1) {
  assert_that(all(scale[!is.na(scale)] > 0))
  assert_that(length(location) == length(scale))
  assert_that(sum(is.na(location)) == sum(is.na(scale)))
  x <- matrix(c(location, scale), ncol = 2L)
  `attr<-`(x, "dist", "normal")
}

#' @export
#' @rdname priors
#' @examples
#' halfnormal(0, 1)
halfnormal <- function(location = 0, scale = 1) {
  normal(location, scale)
}

#' @export
#' @rdname priors
#' @param cv Coefficient of variation (SD/mean).
#' @examples
#' gamma_cv(0.5, 0.2)
gamma_cv <- function(location, cv) {
  assert_that(all(cv[!is.na(cv)] > 0))
  assert_that(length(location) == length(cv))
  assert_that(sum(is.na(location)) == sum(is.na(cv)))
  ## dgamma(shape, scale) in TMB
  ## check:
  # mu <- 1.8;cv <- 0.3
  # x <- rgamma(1e6, shape = 1/cv^2, scale = cv^2*mu)
  # mean(x);sd(x) / mean(x)
  x <- matrix(c(1/cv^2, cv^2*location), ncol = 2L)
  `attr<-`(x, "dist", "gamma")
}

#' @export
#' @rdname priors
#' @param meanlog Mean of the distribution on the log scale.
#' @param sdlog Standard deviation of the distribution on the log scale.
#' @details
#' `lognormal_prior()` defines a lognormal prior with the same parameterization
#' as [stats::dlnorm()]. The median is `exp(meanlog)`.
#' @examples
#' lognormal_prior(log(0.2), 0.5)
lognormal_prior <- function(meanlog, sdlog) {
  assert_that(all(sdlog[!is.na(sdlog)] > 0))
  assert_that(length(meanlog) == length(sdlog))
  assert_that(sum(is.na(meanlog)) == sum(is.na(sdlog)))
  x <- matrix(c(meanlog, sdlog), ncol = 2L)
  `attr<-`(x, "dist", "lognormal")
}

#' @export
#' @rdname priors
#' @examples
#' mvnormal(c(0, 0))
mvnormal <- function(location = 0, scale = diag(length(location))) {
  assert_that(length(location) == dim(scale)[1])
  # return single matrix, where first col = locations, rest = Sigma
  x <- cbind(as.matrix(location,ncol=1), scale)
  `attr<-`(x, "dist", "mvnormal")
}

#' @param range_gt A value one expects the spatial or spatiotemporal range is
#'   **g**reater **t**han with `1 - range_prob` probability.
#' @param sigma_lt A value one expects the spatial or spatiotemporal marginal
#'   standard deviation (`sigma_O` or `sigma_E` internally) is **l**ess **t**han
#'   with `1 - sigma_prob` probability.
#' @param range_prob Probability. See description for `range_gt`.
#' @param sigma_prob Probability. See description for `sigma_lt`.
#' @export
#' @rdname priors
#'
#' @seealso
#' [plot_pc_matern()]
#'
#' @description
#' `pc_matern()` is the Penalized Complexity prior for the Matérn
#' covariance function.
#'
#' @examples
#' pc_matern(range_gt = 5, sigma_lt = 1)
#' plot_pc_matern(range_gt = 5, sigma_lt = 1)
#'
#' \donttest{
#' d <- subset(pcod, year > 2011)
#' pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
#'
#' # - no priors on population-level effects (`b`)
#' # - halfnormal(0, 10) prior on dispersion parameter `phi`
#' # - Matern PC priors on spatial `matern_s` and spatiotemporal
#' #   `matern_st` random field parameters
#' m <- sdmTMB(density ~ s(depth, k = 3),
#'   data = d, mesh = pcod_spde, family = tweedie(),
#'   share_range = FALSE, time = "year",
#'   priors = sdmTMBpriors(
#'     phi = halfnormal(0, 10),
#'     matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
#'     matern_st = pc_matern(range_gt = 5, sigma_lt = 1)
#'   )
#' )
#'
#' # - no prior on intercept
#' # - normal(0, 1) prior on depth coefficient
#' # - no prior on the dispersion parameter `phi`
#' # - Matern PC prior
#' m <- sdmTMB(density ~ depth_scaled,
#'   data = d, mesh = pcod_spde, family = tweedie(),
#'   spatiotemporal = "off",
#'   priors = sdmTMBpriors(
#'     b = normal(c(NA, 0), c(NA, 1)),
#'     matern_s = pc_matern(range_gt = 5, sigma_lt = 1)
#'   )
#' )
#'
#' # You get a prior, you get a prior, you get a prior!
#' # (except on the annual means; see the `NA`s)
#' m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
#'   data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"),
#'   share_range = FALSE, spatiotemporal = "AR1",
#'   priors = sdmTMBpriors(
#'     b = normal(c(0, 0, NA, NA, NA), c(2, 2, NA, NA, NA)),
#'     phi = halfnormal(0, 10),
#'     # tweedie_p = normal(1.5, 2),
#'     ar1_rho = normal(0, 1),
#'     matern_s = pc_matern(range_gt = 5, sigma_lt = 1),
#'     matern_st = pc_matern(range_gt = 5, sigma_lt = 1))
#' )
#' }
pc_matern <- function(range_gt, sigma_lt, range_prob = 0.05, sigma_prob = 0.05) {
  assert_that(range_prob > 0 && range_prob < 1)
  assert_that(sigma_prob > 0 && sigma_prob < 1)
  if (!is.na(range_gt)) assert_that(range_gt > 0)
  if (!is.na(sigma_lt)) assert_that(sigma_lt > 0)
  if (!is.na(range_gt)) assert_that(!is.na(sigma_lt))
  if (!is.na(sigma_lt)) assert_that(!is.na(range_gt))
  x <- c(range_gt, sigma_lt, range_prob, sigma_prob)
  `attr<-`(x, "dist", "pc_matern")
}
