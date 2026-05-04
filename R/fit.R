#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal Gaussian random field generalized linear
#' mixed effects model (GLMM) with the TMB (Template Model Builder) R package and
#' the SPDE (stochastic partial differential equation) approach. This can be
#' useful for (dynamic) species distribution models and relative abundance index
#' standardization among many other uses.
#'
#' @param formula Model formula. IID random intercepts are possible using
#'   \pkg{lme4} syntax, e.g., `+ (1 | g)` where `g` is a column of class
#'   character or factor representing groups. Penalized splines are possible via
#'   \pkg{mgcv} with `s()`. Optionally a list for delta (hurdle) models.  See
#'   examples and details below.
#' @param data A data frame.
#' @param mesh An object from [make_mesh()].
#' @param time An optional time column name (as character). Can be left as
#'   `NULL` for a model with only spatial random fields; however, if the data
#'   are actually spatiotemporal and you wish to use [get_index()] or [get_cog()]
#'   downstream, supply the time argument.
#' @param family The family and link. Supports [gaussian()], [Gamma()],
#'   [binomial()], [poisson()], \code{\link[sdmTMB:families]{Beta()}},
#'   \code{\link[sdmTMB:families]{nbinom2()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom2()}},
#'   \code{\link[sdmTMB:families]{nbinom1()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom1()}},
#'   \code{\link[sdmTMB:families]{censored_poisson()}},
#'   \code{\link[sdmTMB:families]{gamma_mix()}},
#'   \code{\link[sdmTMB:families]{lognormal_mix()}},
#'   \code{\link[sdmTMB:families]{student()}}, and
#'   \code{\link[sdmTMB:families]{tweedie()}}. Supports the delta/hurdle models:
#'   \code{\link[sdmTMB:families]{delta_beta()}},
#'   \code{\link[sdmTMB:families]{delta_gamma()}},
#'   \code{\link[sdmTMB:families]{delta_gamma_mix()}},
#'   \code{\link[sdmTMB:families]{delta_lognormal_mix()}},
#'   \code{\link[sdmTMB:families]{delta_lognormal()}}, and
#'   \code{\link[sdmTMB:families]{delta_truncated_nbinom2()}},
#'   For binomial family options, see 'Binomial families' in the Details
#'   section below.
#' @param spatial Estimate spatial random fields? Options are `'on'` / `'off'`
#'   or `TRUE` / `FALSE`. Optionally, a list for delta models, e.g. `list('on',
#'   'off')`.
#' @param spatiotemporal Estimate the spatiotemporal random fields as `'iid'`
#'   (independent and identically distributed; default), stationary `'ar1'`
#'   (first-order autoregressive), as a random walk (`'rw'`), or as fixed at 0
#'   `'off'`. Will be set to `'off'` if `time = NULL`. If a delta model, can be
#'   a list. E.g., `list('off', 'ar1')`. Note that the spatiotemporal standard
#'   deviation represents the marginal steady-state standard deviation of the
#'   process in the case of the AR1. I.e., it is scaled according to the
#'   correlation. See the [TMB
#'   documentation](https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html). If the AR1
#'   correlation coefficient (rho) is estimated close to 1, say > 0.99, then you
#'   may wish to switch to the random walk `'rw'`. Capitalization is ignored.
#'   `TRUE` gets converted to `'iid'` and `FALSE` gets converted to `'off'`.
#' @param share_range Logical: estimate a shared spatial and spatiotemporal
#'   range parameter (`TRUE`, default) or independent range parameters (`FALSE`). If a
#'   delta model, can be a list. E.g., `list(TRUE, FALSE)`.
#' @param time_varying An optional one-sided formula describing covariates that
#'   should be modelled as a random walk through time. Be careful not to include
#'   covariates (including the intercept) in both the main and time-varying
#'   formula since the first time step is estimated independently. I.e., at
#'   least one should have `~ 0` or `~ -1`. Structure must currently be shared
#'   in delta models.
#' @param time_varying_type Type of time-varying process to apply to
#'   `time_varying` formula.
#' @param spatial_varying An optional one-sided formula of coefficients that
#'   should vary in space as random fields. Note that you likely want to include
#'   a fixed effect for the same variable to improve interpretability since the
#'   random field is assumed to have a mean of 0. If a (scaled) time column is
#'   used, it will represent a local-time-trend model. See
#'   \doi{10.1111/ecog.05176} and the [spatial trends
#'   vignette](https://pbs-assess.github.io/sdmTMB/articles/spatial-trend-models.html).
#'   Note this predictor should be centered to have mean zero and have a
#'   standard deviation of approximately 1 and should likely also be included as
#'   a main effect. Structure must currently be shared in delta models.
#' @param weights A numeric vector representing optional likelihood weights for
#'   the conditional model. Implemented as in \pkg{glmmTMB}: weights do not have
#'   to sum to one and are not internally modified. Can also be used for trials
#'   with the binomial family; the `weights` argument needs to be a vector and not
#'   a name of the variable in the data frame. See the Details section below.
#' @param offset A numeric vector representing the model offset. In delta/hurdle
#'   models, this applies only to the positive component. Usually a log
#'   transformed variable. The `offset` argument needs to be a vector and not a name of the variable in the data frame.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See the Details
#'   section below.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation
#'   rather than maximum likelihood? Internally, this adds the fixed effects to
#'   the list of random effects to integrate over.
#' @param silent Silent or include optimization details? Helpful to set to
#'   `FALSE` for models that take a while to fit.
#' @param anisotropy Logical: allow for anisotropy (spatial correlation that is
#'   directionally dependent)? See [plot_anisotropy()].
#'   Must be shared across delta models.
#' @param control Optimization control options via [sdmTMBcontrol()].
#' @param priors Optional penalties/priors via [sdmTMBpriors()]. Must currently
#'   be shared across delta models.
#' @param knots Optional named list containing knot values to be used for basis
#'   construction of smoothing terms. See [mgcv::gam()] and [mgcv::gamm()].
#'   E.g., `s(x, bs = 'cc', k = 4), knots = list(x = c(1, 2, 3, 4))`
#' @param previous_fit A previously fitted sdmTMB model to initialize the
#'   optimization with. Can greatly speed up fitting. Note that the model must
#'   be set up *exactly* the same way. However, the data and `weights` arguments
#'   can change, which can be useful for cross-validation.
#' @param do_fit Fit the model (`TRUE`) or return the processed data without
#'   fitting (`FALSE`)?
#' @param do_index Do index standardization calculations while fitting? Saves
#'   memory and time when working with large datasets or projection grids since
#'   the TMB object doesn't have to be rebuilt with [predict.sdmTMB()] and
#'   [get_index()]. If `TRUE`, then `predict_args` must have a `newdata` element
#'   supplied and `area` can be supplied to `index_args`.
#' @param predict_args A list of arguments to pass to [predict.sdmTMB()] if
#'   `do_index = TRUE`.
#' @param index_args A list of arguments to pass to [get_index()] if
#'   `do_index = TRUE`. Currently, only `area` is supported. Bias correction
#'   can be done when calling [get_index()] on the resulting fitted object.
#' @param bayesian Logical indicating if the model will be passed to
#'   \pkg{tmbstan}. If `TRUE`, Jacobian adjustments are applied to account for
#'   parameter transformations when priors are applied.
#' @param experimental A named list for esoteric or in-development options. Here
#'   be dragons.
#   (Experimental) A column name (as character) of a predictor of a
#   linear trend (in log space) of the spatiotemporal standard deviation. By
#   default, this is `NULL` and fits a model with a constant spatiotemporal
#   variance. However, this argument can also be a character name in the
#   original data frame (a covariate that ideally has been standardized to have
#   mean 0 and standard deviation = 1). Because the spatiotemporal field varies
#   by time step, the standardization should be done by time. If the name of a
#   predictor is included, a log-linear model is fit where the predictor is
#   used to model effects on the standard deviation, e.g. `log(sd(i)) = B0 + B1
#   * epsilon_predictor(i)`. The 'epsilon_model' argument may also be
#   specified. This is the name of the model to use to modeling time-varying
#   epsilon. This can be one of the following: "trend" (default, fits a linear
#   model without random effects), "re" (fits a model with random effects in
#   epsilon_st, but no trend), and "trend-re" (a model that includes both the
#   trend and random effects)
#' @importFrom methods as is
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom mgcv s t2
#' @importFrom stats gaussian model.frame model.matrix as.formula
#'   model.response terms model.offset
#' @importFrom lifecycle deprecated is_present deprecate_warn deprecate_stop
#'
#' @return
#' An object (list) of class `sdmTMB`. Useful elements include:
#'
#' * `sd_report`: output from [TMB::sdreport()]
#' * `gradients`: log likelihood gradients with respect to each fixed effect
#' * `model`: output from [stats::nlminb()]
#' * `data`: the fitted data
#' * `mesh`: the object that was supplied to the `mesh` argmument
#' * `family`: the family object, which includes the inverse link function as `family$linkinv()`
#' * `tmb_params`: The parameters list passed to [TMB::MakeADFun()]
#' * `tmb_map`: The 'map' list passed to [TMB::MakeADFun()]
#' * `tmb_data`: The data list passed to [TMB::MakeADFun()]
#' * `tmb_obj`: The TMB object created by [TMB::MakeADFun()]
#'
#' @details
#'
#' **Model description**
#'
#' See the [model description](https://pbs-assess.github.io/sdmTMB/articles/model-description.html)
#' vignette or the relevant appendix of the preprint on sdmTMB:
#' \doi{10.1101/2022.03.24.485545}
#'
#' **Binomial families**
#'
#' Following the structure of [stats::glm()] and \pkg{glmmTMB}, a binomial
#' family can be specified in one of 4 ways: (1) the response may be a factor
#' (and the model classifies the first level versus all others), (2) the
#' response may be binomial (0/1), (3) the response can be a matrix of form
#' `cbind(success, failure)`, and (4) the response may be the observed
#' proportions, and the 'weights' argument is used to specify the Binomial size
#' (N) parameter (`prob ~ ..., weights = N`).
#'
#' **Smooth terms**
#'
#' Smooth terms can be included following GAMs (generalized additive models)
#' using `+ s(x)`, which implements a smooth from [mgcv::s()]. \pkg{sdmTMB} uses
#' penalized smooths, constructed via [mgcv::smooth2random()]. This is a similar
#' approach implemented in \pkg{gamm4} and \pkg{brms}, among other packages.
#' Within these smooths, the same syntax commonly used in [mgcv::s()] can be
#' applied, e.g. 2-dimensional smooths may be constructed with `+ s(x, y)`;
#' smooths can be specific to various factor levels, `+ s(x, by = group)`; the
#' basis function dimensions may be specified, e.g. `+ s(x, k = 4)`; and various
#' types of splines may be constructed such as cyclic splines to model
#' seasonality, `+ s(month, bs = "cc", k = 12)`. Prior to version 0.0.18,
#' \pkg{sdmTMB} implemented unpenalized splines.
#'
#' **Threshold models**
#'
#' A linear break-point relationship for a covariate can be included via
#' `+ breakpt(variable)` in the formula, where `variable` is a single covariate
#' corresponding to a column in `data`. In this case, the relationship is linear
#' up to a point and then constant (hockey-stick shaped).
#'
#' Similarly, a logistic-function threshold model can be included via
#' `+ logistic(variable)`. This option models the relationship as a logistic
#' function of the 50% and 95% values. This is similar to length- or size-based
#' selectivity in fisheries, and is parameterized by the points at which f(x) =
#' 0.5 or 0.95. See the [threshold vignette](https://pbs-assess.github.io/sdmTMB/articles/threshold-models.html).
#'
#' A special metabolic-index threshold can be fit with `+ logistic(mi)`. Here,
#' `mi` is a reserved name rather than a required column in `data`. \pkg{sdmTMB}
#' estimates the temperature-sensitivity parameter `Eo` and internally computes
#' the metabolic index as `po2 * exp(Eo * invtemp)`. This requires columns named
#' `po2` for oxygen partial pressure and `invtemp` for inverse temperature in
#' both the fitting data and any prediction data. The oxygen and inverse
#' temperature columns can appear in any order in the data frame.
#'
#' Note that only a single threshold covariate can be included and the same covariate
#' is included in both components for the delta families.
#'
#' **Extra time: forecasting or interpolating**
#'
#' Extra time slices (e.g., years) can be included for interpolation or
#' forecasting with the predict function via the `extra_time` argument. The
#' predict function requires all time slices to be defined when fitting the
#' model to ensure the various time indices are set up correctly. Be careful if
#' including extra time slices that the model remains identifiable. For example,
#' including `+ as.factor(year)` in `formula` will render a model with no data
#' to inform the expected value in a missing year. [sdmTMB()] makes no attempt
#' to determine if the model makes sense for forecasting or interpolation. The
#' options `time_varying`, `spatiotemporal = "rw"`, and `spatiotemporal = "ar1"`
#' provide mechanisms to predict over missing time slices with process error.
#'
#' `extra_time` can also be used to fill in missing time steps for the purposes
#' of a random walk or AR1 spatiotemporal field if their inclusion makes the gaps
#' between time steps even.
#'
#' **Index standardization**
#'
#' For index standardization, you may wish to include `0 + as.factor(year)`
#' (or whatever the time column is called) in the formula. See a basic
#' example of index standardization in the relevant
#' [package vignette](https://pbs-assess.github.io/sdmTMB/articles/index-standardization.html).
#' You will need to specify the `time` argument. See [get_index()] and/or
#' [get_index_sims()].
#'
#' **Regularization and priors**
#'
#' You can achieve regularization via penalties (priors) on the fixed effect
#' parameters. See [sdmTMBpriors()]. You can fit the model once without
#' penalties and look at the output of `print(your_model)` or `tidy(your_model)`
#' or fit the model with `do_fit = FALSE` and inspect
#' `head(your_model$tmb_data$X_ij[[1]])` if you want to see how the formula is
#' translated to the fixed effect model matrix. Also see the
#' [Bayesian vignette](https://pbs-assess.github.io/sdmTMB/articles/bayesian.html).
#'
#' **Delta/hurdle models**
#'
#' Delta models (also known as hurdle models) can be fit as two separate models
#' or at the same time by using an appropriate delta family. E.g.:
#'   \code{\link[sdmTMB:families]{delta_gamma()}},
#'   \code{\link[sdmTMB:families]{delta_beta()}},
#'   \code{\link[sdmTMB:families]{delta_lognormal()}}, and
#'   \code{\link[sdmTMB:families]{delta_truncated_nbinom2()}}.
#' If fit with a delta family, by default the formula, spatial, and spatiotemporal
#' components are shared. Some elements can be specified independently for the two models
#' using a list format. These include `formula`, `spatial`, `spatiotemporal`,
#' and `share_range`. The first element of the list is for the binomial component
#' and the second element is for the positive component (e.g., Gamma).
#' Other elements must be shared for now (e.g., spatially varying coefficients,
#' time-varying coefficients). Furthermore, there are currently limitations if
#' specifying two formulas as a list: the two formulas cannot have smoothers,
#' threshold effects, or random intercepts. For now, these must be specified
#' through a single formula that is shared across the two models.
#'
#' The main advantage of specifying such models using a delta family (compared
#' to fitting two separate models) is (1) coding simplicity and (2) calculation
#' of uncertainty on derived quantities such as an index of abundance with
#' [get_index()] using the generalized delta method within TMB. Also, parameters
#' can be shared across the models.
#'
#' See the [delta-model vignette](https://pbs-assess.github.io/sdmTMB/articles/delta-models.html).
#'
#' @references
#'
#' *Main reference introducing the package to cite when using sdmTMB:*
#'
#' Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett. 2022. sdmTMB: an R
#' package for fast, flexible, and user-friendly generalized linear mixed effects
#' models with spatial and spatiotemporal random fields.
#' bioRxiv 2022.03.24.485545; \doi{10.1101/2022.03.24.485545}.
#'
#' *Reference for local trends:*
#'
#' Barnett, L.A.K., E.J. Ward, S.C. Anderson. Improving estimates of species
#' distribution change by incorporating local trends. Ecography. 44(3):427-439.
#' \doi{10.1111/ecog.05176}.
#'
#' *Further explanation of the model and application to calculating climate
#' velocities:*
#'
#' English, P., E.J. Ward, C.N. Rooper, R.E. Forrest, L.A. Rogers, K.L. Hunter,
#' A.M. Edwards, B.M. Connors, S.C. Anderson. 2021. Contrasting climate velocity
#' impacts in warm and cool locations show that effects of marine warming are
#' worse in already warmer temperate waters. In press at Fish and Fisheries.
#' \doi{10.1111/faf.12613}.
#'
#' *Discussion of and illustration of some decision points when fitting these
#' models:*
#'
#' Commander, C.J.C., Barnett, L.A.K., Ward, E.J., Anderson, S.C., and
#' Essington, T.E. 2022. The shadow model: how and why small choices in
#' spatially explicit species distribution models affect predictions. PeerJ 10:
#' e12783. \doi{10.7717/peerj.12783}.
#'
#' *Application and description of the threshold/break-point models:*
#'
#' Essington, T.E. S.C. Anderson, L.A.K. Barnett, H.M. Berger, S.A. Siedlecki,
#' E.J. Ward. Advancing statistical models to reveal the effect of dissolved
#' oxygen on the spatial distribution of marine taxa using thresholds and a
#' physiologically based index. In press at Ecography. \doi{10.1111/ecog.06249}.
#'
#' *Application to fish body condition:*
#'
#' Lindmark, M., S.C. Anderson, M. Gogina, M. Casini. Evaluating drivers of
#' spatiotemporal individual condition of a bottom-associated marine fish.
#' bioRxiv 2022.04.19.488709. \doi{10.1101/2022.04.19.488709}.
#'
#' *A number of sections of the original TMB model code were adapted from the
#' VAST R package:*
#'
#' Thorson, J.T., 2019. Guidance for decisions using the Vector Autoregressive
#' Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate
#' assessments. Fish. Res. 210:143–161.
#' \doi{10.1016/j.fishres.2018.10.013}.
#'
#' *Code for the `family` R-to-TMB implementation, selected parameterizations of
#' the observation likelihoods, general package structure inspiration, and the
#' idea behind the TMB prediction approach were adapted from the glmmTMB R
#' package:*
#'
#' Mollie E. Brooks, Kasper Kristensen, Koen J. van Benthem, Arni Magnusson,
#' Casper W. Berg, Anders Nielsen, Hans J. Skaug, Martin Maechler and Benjamin
#' M. Bolker (2017). glmmTMB Balances Speed and Flexibility Among Packages for
#' Zero-inflated Generalized Linear Mixed Modeling. The R Journal, 9(2):378-400.
#' \doi{10.32614/rj-2017-066}.
#'
#' *Implementation of geometric anisotropy with the SPDE and use of random
#' field GLMMs for index standardization*:
#'
#' Thorson, J.T., Shelton, A.O., Ward, E.J., and Skaug, H.J. 2015.
#' Geostatistical delta-generalized linear mixed models improve precision for
#' estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci.
#' 72(5): 1297–1310. \doi{10.1093/icesjms/fsu243}.
#'
#' @export
#'
#' @examplesIf inla_installed() && require("visreg", quietly = TRUE)
#' library(sdmTMB)
#'
#' # Build an SPDE mesh with INLA:
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' # * this example uses a fairly coarse mesh so these examples run quickly
#' # * `cutoff` is the minimum distance between mesh vertices in units of the
#' #   x and y coordinates
#' # * `cutoff = 10` or `cutoff = 15` might make more sense in applied situations
#' #   for this dataset
#' # * or build any mesh in INLA and pass it to the `mesh` argument in `make_mesh()`
#' # * not needed if you will be turning off all spatial/spatiotemporal random fields
#'
#' # Quick mesh plot:
#' plot(mesh)
#'
#' # Fit a Tweedie spatial random field GLMM with a smoother for depth:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Extract coefficients:
#' tidy(fit, conf.int = TRUE)
#' tidy(fit, effects = "ran_par", conf.int = TRUE)
#'
#' # Perform several 'sanity' checks:
#' sanity(fit)
#'
#' # Visualize depth effect: (see ?visreg_delta)
#' visreg::visreg(fit, xvar = "depth") # link space; randomized quantile residuals
#' visreg::visreg(fit, xvar = "depth", scale = "response")
#' visreg::visreg(fit, xvar = "depth", scale = "response", gg = TRUE, rug = FALSE)
#'
#' # Predict on the fitted data; see ?predict.sdmTMB
#' p <- predict(fit)
#'
#' # Predict on new data:
#' nd <- subset(qcs_grid, year == 2017)
#' p <- predict(fit, newdata = nd)
#' head(p)
#'
#' # Add spatiotemporal random fields:
#' fit <- sdmTMB(
#'   density ~ 0 + as.factor(year),
#'   time = "year", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Make the fields AR1:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "ar1", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Make the fields a random walk:
#' fit <- sdmTMB(
#'   density ~ s(depth),
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "rw", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Depth smoothers by year:
#' fit <- sdmTMB(
#'   density ~ s(depth, by = as.factor(year)), #<
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "rw",
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # 2D depth-year smoother:
#' fit <- sdmTMB(
#'   density ~ s(depth, year), #<
#'   spatial = "off",
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Turn off spatial random fields:
#' fit <- sdmTMB(
#'   present ~ poly(log(depth)),
#'   spatial = "off", #<
#'   data = pcod_2011, mesh = mesh,
#'   family = binomial()
#' )
#' fit
#'
#' # Which, matches glm():
#' fit_glm <- glm(
#'   present ~ poly(log(depth)),
#'   data = pcod_2011,
#'   family = binomial()
#' )
#' summary(fit_glm)
#' AIC(fit, fit_glm)
#'
#' # Delta/hurdle binomial-Gamma model:
#' fit_dg <- sdmTMB(
#'   density ~ poly(log(depth), 2),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = delta_gamma() #<
#' )
#' fit_dg
#'
#' # Delta model with different formulas and spatial structure:
#' fit_dg <- sdmTMB(
#'   list(density ~ depth_scaled, density ~ poly(depth_scaled, 2)), #<
#'   data = pcod_2011, mesh = mesh,
#'   spatial = list("off", "on"), #<
#'   family = delta_gamma()
#' )
#' fit_dg
#'
#' # Delta/hurdle truncated NB2:
#' pcod_2011$count <- round(pcod_2011$density)
#' fit_nb2 <- sdmTMB(
#'   count ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = delta_truncated_nbinom2() #<
#' )
#' fit_nb2
#'
#' # Regular NB2:
#' fit_nb2 <- sdmTMB(
#'   count ~ s(depth),
#'   data = pcod_2011, mesh = mesh,
#'   spatial = "off",
#'   family = nbinom2() #<
#' )
#' fit_nb2
#'
#' # IID random intercepts by year:
#' pcod_2011$fyear <- as.factor(pcod_2011$year)
#' fit <- sdmTMB(
#'   density ~ s(depth) + (1 | fyear), #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
#' fit
#'
#' # Spatially varying coefficient of year:
#' pcod_2011$year_scaled <- as.numeric(scale(pcod_2011$year))
#' fit <- sdmTMB(
#'   density ~ year_scaled,
#'   spatial_varying = ~ 0 + year_scaled, #<
#'   data = pcod_2011, mesh = mesh, family = tweedie(), time = "year"
#' )
#' fit
#'
#' # Time-varying effects of depth and depth squared:
#' fit <- sdmTMB(
#'   density ~ 0 + as.factor(year),
#'   time_varying = ~ 0 + depth_scaled + depth_scaled2, #<
#'   data = pcod_2011, time = "year", mesh = mesh,
#'   family = tweedie()
#' )
#' print(fit)
#' # Extract values:
#' est <- as.list(fit$sd_report, "Estimate")
#' se <- as.list(fit$sd_report, "Std. Error")
#' est$b_rw_t[, , 1]
#' se$b_rw_t[, , 1]
#'
#' # Linear break-point effect of depth:
#' fit <- sdmTMB(
#'   density ~ breakpt(depth_scaled), #<
#'   data = pcod_2011,
#'   mesh = mesh,
#'   family = tweedie()
#' )
#' fit

sdmTMB <- function(
  formula,
  data,
  mesh,
  time = NULL,
  family = gaussian(link = "identity"),
  spatial = c("on", "off"),
  spatiotemporal = c("iid", "ar1", "rw", "off"),
  share_range = TRUE,
  time_varying = NULL,
  time_varying_type = c("rw", "ar1"),
  spatial_varying = NULL,
  weights = NULL,
  offset = NULL,
  extra_time = NULL,
  reml = FALSE,
  silent = TRUE,
  anisotropy = FALSE,
  control = sdmTMBcontrol(),
  priors = sdmTMBpriors(),
  knots = NULL,
  bayesian = FALSE,
  previous_fit = NULL,
  do_fit = TRUE,
  do_index = FALSE,
  predict_args = NULL,
  index_args = NULL,
  experimental = NULL
  ) {


  delta <- isTRUE(family$delta)
  n_m <- if (delta) 2L else 1L

  if (!missing(spatial)) {
    if (length(spatial) > 1 && !is.list(spatial)) {
      cli_abort("`spatial` should be a single value or a list")
    }
  }
  if (!missing(spatiotemporal)) {
    if (length(spatiotemporal) > 1 && !is.list(spatiotemporal)) {
      cli_abort("`spatiotemporal` should be a single value or a list")
    }
    if (delta && !is.list(spatiotemporal)) {
      spatiotemporal <- rep(spatiotemporal[[1]], 2L)
    }
    spatiotemporal <- vapply(seq_along(spatiotemporal),
      function(i) check_spatiotemporal_arg(
        spatiotemporal, time = time, .which = i), FUN.VALUE = character(1L))
  } else {
    if (is.null(time))
      spatiotemporal <- rep("off", n_m)
    else
      spatiotemporal <- rep("iid", n_m)
  }


  if (is.null(time)) {
    spatial_only <- rep(TRUE, n_m)
  } else {
    spatial_only <- ifelse(spatiotemporal == "off", TRUE, FALSE)
    # if(all(spatiotemporal == "off")) {
    #   cli_abort("Time needs to be null if spatiotemporal fields are not included")
    # }
  }

  if (is.list(spatial)) {
    spatial <- vapply(spatial, parse_spatial_arg, FUN.VALUE = character(1L))
  } else {
    spatial <- rep(parse_spatial_arg(spatial), n_m)
  }
  include_spatial <- "on" %in% spatial

  if (!include_spatial && all(spatiotemporal == "off") || !include_spatial && all(spatial_only)) {
    # message("Both spatial and spatiotemporal fields are set to 'off'.")
    # control$map_rf <- TRUE
    no_spatial <- TRUE
    if (missing(mesh)) {
      mesh <- sdmTMB::pcod_mesh_2011 # internal data; fake!
    }
  } else {
    no_spatial <- FALSE
  }

  share_range <- unlist(share_range)
  if (length(share_range) == 1L) share_range <- rep(share_range, n_m)
  share_range[spatiotemporal == "off"] <- TRUE

  spde <- mesh
  epsilon_model <- NULL
  epsilon_predictor <- NULL
  if (!is.null(experimental)) {
    if ("epsilon_predictor" %in% names(experimental)) {
      epsilon_predictor <- experimental$epsilon_predictor
    } else {
      epsilon_predictor <- NULL
    }

    if ("epsilon_model" %in% names(experimental)) {
      epsilon_model <- experimental$epsilon_model
    } else {
      epsilon_model <- NULL
    }

    if ("lwr" %in% names(experimental) && "upr" %in% names(experimental)) {
      lwr <- experimental$lwr
      upr <- experimental$upr
    } else {
      #epsilon_predictor <- NULL
      lwr <- 0
      upr <- Inf
    }
  } else {
    #epsilon_predictor <- NULL
    lwr <- 0
    upr <- Inf
  }
  normalize <- control$normalize
  nlminb_loops <- control$nlminb_loops
  newton_loops <- control$newton_loops
  quadratic_roots <- control$quadratic_roots
  start <- control$start
  multiphase <- control$multiphase
  map <- control$map
  lower <- control$lower
  upper <- control$upper
  get_joint_precision <- control$get_joint_precision

  dot_checks <- c("lower", "upper", "profile", "parallel",
    "nlminb_loops", "newton_steps", "mgcv", "quadratic_roots", "multiphase",
    "newton_loops", "start", "map", "get_joint_precision", "normalize")
  .control <- control
  # FIXME; automate this from sdmTMcontrol args?
  for (i in dot_checks) .control[[i]] <- NULL # what's left should be for nlminb

  ar1_fields <- spatiotemporal == "ar1"
  rw_fields <- spatiotemporal == "rw"
  assert_that(
    is.logical(reml), is.logical(anisotropy), is.logical(share_range),
    is.logical(silent), is.logical(multiphase), is.logical(normalize)
  )

  if (!is.null(spatial_varying)) assert_that(class(spatial_varying) %in% c("formula","list"))
  if (!is.null(time_varying)) assert_that(class(time_varying) %in% c("formula","list"))
  if (!is.null(previous_fit)) assert_that(identical(class(previous_fit), "sdmTMB"))
  assert_that(is.list(priors))
  assert_that(is.list(.control))
  if (!is.null(time)) assert_that(is.character(time))
  assert_that(inherits(spde, "sdmTMBmesh"))
  assert_that(class(formula) %in% c("formula", "list"))
  assert_that(inherits(data, "data.frame"))
  time_varying_type <- match.arg(time_varying_type)
  if (!is.null(map) && length(map) != length(start)) {
    cli_warn(c("`length(map) != length(start)`.",
      "You likely want to specify `start` values if you are setting the `map` argument."))
  }
  if (family$family[1] == "censored_poisson") {
    assert_that("lwr" %in% names(experimental) && "upr" %in% names(experimental),
      msg = "`lwr` and `upr` must be specified in `experimental` as elements of a named list to use the censored Poisson likelihood.")
    assert_that(length(lwr) == nrow(data) && length(upr) == nrow(data))
    assert_that(length(lwr) == length(upr))
    assert_that(mean(upr-lwr, na.rm = TRUE)>=0)
  }

  if (!is.null(time)) {
    assert_that(time %in% names(data),
      msg = "Specified `time` column is missing from `data`.")
  }
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1)
      cli_abort("There is at least one NA value in the time column. Please remove it.")
  }
  if (is.factor(data[[time]])) {
    if (length(levels(data[[time]])) > length(unique(data[[time]]))) {
      cli_abort("The time column is a factor and there are extra factor levels. ",
        "Please remove these or turn your time column into an integer.")
    }
  }

  if (!no_spatial) {
    if (!identical(nrow(spde$loc_xy), nrow(data))) {
      msg <- c(
        "Number of x-y coordinates in `mesh` does not match `nrow(data)`.",
        "Is it possible you passed a different data frame to `make_mesh()` and `sdmTMB()`?" # discussion 137
      )
      cli::cli_abort(msg)
    }
  }
  # FIXME parallel setup here?

  # thresholds not yet enabled for delta model, where formula is a list
  if (inherits(formula, "formula")) {
    original_formula <- replicate(n_m, list(formula))
    thresh <- list(check_and_parse_thresh_params(formula, data))
    if (delta) {
      formula <- list(thresh[[1]]$formula, thresh[[1]]$formula)
    } else {
      formula <- list(thresh[[1]]$formula)
    }
  } else {
    original_formula <- formula
    thresh <- list(check_and_parse_thresh_params(formula[[1]], data),
                   check_and_parse_thresh_params(formula[[2]], data))
    formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
  }

  if (!is.null(extra_time)) { # for forecasting or interpolating
    data[["__sdmTMB_offset__"]] <- offset
    data <- expand_time(df = data, time_slices = extra_time, time_column = time)
    offset <- data[["__sdmTMB_offset__"]]
    data[["__sdmTMB_offset__"]] <- NULL
    weights <- data$weight_sdmTMB
    spde$loc_xy <- as.matrix(data[,spde$xy_cols,drop=FALSE])
    spde$A_st <- INLA::inla.spde.make.A(spde$mesh, loc = spde$loc_xy)
    spde$sdm_spatial_id <- seq(1, nrow(data)) # FIXME
  }
  check_irregalar_time(data, time, spatiotemporal, time_varying)

  spatial_varying_formula <- spatial_varying # save it
  if (!is.null(spatial_varying)) {
    z_i <- model.matrix(spatial_varying, data)
    .int <- grep("(Intercept)", colnames(z_i))
    if (sum(.int) > 0) z_i <- z_i[,-.int,drop=FALSE]
    spatial_varying <- colnames(z_i)
  } else {
    z_i <- matrix(0, nrow(data), 0L)
  }
  n_z <- ncol(z_i)

  if (any(grepl("offset\\(", formula)))
    cli_abort("Detected `offset()` in formula. Offsets in sdmTMB must be specified via the `offset` argument.")
  contains_offset <- check_offset(formula[[1]]) # deprecated check

  split_formula <- list() # passed to out structure, not TMB
  RE_indexes <- list() # ncols passed into TMB
  nobs_RE <- list() # ncols passed into TMB
  ln_tau_G_index<- list() # passed into TMB
  X_ij = list() # main effects, passed into TMB
  mf <- list()
  mt <- list()
  sm <- list()

  for (ii in seq_along(formula)) {
    contains_offset <- check_offset(formula[[ii]])

    # anything in a list here needs to be saved for tmb data
    split_formula[[ii]] <- glmmTMB::splitForm(formula[ii][[1]])
    RE_names <- barnames(split_formula[[ii]]$reTrmFormulas)

    fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
    RE_indexes[[ii]] <- vapply(RE_names, function(x) as.integer(data[[x]]) - 1L, rep(1L, nrow(data)))
    nobs_RE[[ii]] <- unname(apply(RE_indexes[[ii]], 2L, max)) + 1L
    if (length(nobs_RE[[ii]]) == 0L) nobs_RE[[ii]] <- 0L
    formula[[ii]] <- split_formula[[ii]]$fixedFormula
    ln_tau_G_index[[ii]] <- unlist(lapply(seq_along(nobs_RE[[ii]]), function(i) rep(i, each = nobs_RE[[ii]][i]))) - 1L

    formula_no_sm <- remove_s_and_t2(formula[[ii]])
    X_ij[[ii]] <- model.matrix(formula_no_sm, data)
    mf[[ii]] <- model.frame(formula_no_sm, data)
    # Check for random slopes:
    if (length(split_formula[[ii]]$reTrmFormulas)) {
      termsfun <- function(x) {
        # using glmTMB and lme4:::mkBlist approach
        ff <- eval(substitute(~foo, list(foo = x[[2]])))
        tt <- try(terms(ff, data =  mf[[ii]]), silent = TRUE)
        tt
      }
      reXterms <- lapply(split_formula[[ii]]$reTrmFormulas, termsfun)
      if (length(attr(reXterms[[1]], "term.labels")))
        cli_abort("This model appears to have a random slope specified (e.g., y ~ (1 + b | group)). sdmTMB currently can only do random intercepts (e.g., y ~ (1 | group)).")
    }

    mt[[ii]] <- attr(mf[[ii]], "terms")
    # parse everything mgcv + smoothers:
    sm[[ii]] <- parse_smoothers(formula = formula[[ii]], data = data, knots = knots)
  }

  if (delta) {
    if (any(unlist(lapply(nobs_RE, function(.x) .x > 0)))) {
      if (original_formula[[1]] != original_formula[[2]]) {
        msg <- paste0("For now, if delta models contain random intercepts, both ",
          "components must have the same main-effects formula.")
        cli_abort(msg)
      }
    }
    if (any(unlist(lapply(sm, `[[`, "has_smooths")))) {
      if (original_formula[[1]] != original_formula[[2]]) {
        msg <- paste0("For now, if delta models contain smoothers, both components ",
            "must have the same main-effects formula.")
        cli_abort(msg)
      }
    }
  }

  # split_formula <- split_formula[[1]] # Delete this and next 7 lines as smooths / random effects added
  RE_indexes <- RE_indexes[[1]]
  nobs_RE <- nobs_RE[[1]]
  ln_tau_G_index <- ln_tau_G_index[[1]]
  sm <- sm[[1]]

  y_i <- model.response(mf[[1]], "numeric")
  if (delta) {
    y_i2 <- model.response(mf[[2]], "numeric")
    if (!identical(y_i, y_i2))
      cli_abort("Response variable should be the same in both parts of the delta formula.")
  }

  if (family$family[1] %in% c("Gamma", "lognormal") && min(y_i) <= 0 && !delta) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }

  # This is taken from approach in glmmTMB to match how they handle binomial
  # yobs could be a factor -> treat as binary following glm
  # yobs could be cbind(success, failure)
  # yobs could be binary
  # (yobs, weights) could be (proportions, size)
  # On the C++ side 'yobs' must be the number of successes.
  size <- rep(1, nrow(X_ij[[1]])) # for non-binomial case TODO: change hard coded index
  if (identical(family$family[1], "binomial") && !delta) {
    ## call this to catch the factor / matrix cases
    y_i <- model.response(mf[[1]], type = "any")
    if (is.factor(y_i)) {
      ## following glm, ‘success’ is interpreted as the factor not
      ## having the first level (and hence usually of having the
      ## second level).
      y_i <- pmin(as.numeric(y_i) - 1, 1)
      size <- rep(1, length(y_i))
    } else {
      if (is.matrix(y_i)) { # yobs=cbind(success, failure)
        size <- y_i[, 1] + y_i[, 2]
        yobs <- y_i[, 1] # successes
        y_i <- yobs
      } else {
        if (all(y_i %in% c(0, 1))) { # binary
          size <- rep(1, length(y_i))
        } else { # proportions
          y_i <- weights * y_i
          size <- weights
          weights <- rep(1, length(y_i))
        }
      }
    }
  }

  if (identical(family$link[1], "log") && min(y_i, na.rm = TRUE) < 0 && !delta) {
    cli_abort("`link = 'log'` but the reponse data include values < 0.")
  }

  if (is.null(offset)) offset <- rep(0, length(y_i))
  assert_that(length(offset) == length(y_i), msg = "Offset doesn't match length of data")

  if (!is.null(time_varying)) {
    X_rw_ik <- model.matrix(time_varying, data)
  } else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
  }

  n_s <- nrow(spde$mesh$loc)

  barrier <- "spde_barrier" %in% names(spde)
  if (barrier && anisotropy) {
    cli_warn("Using a barrier mesh; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }
  if (any(c(!is.na(priors$matern_s[1:2]), !is.na(priors$matern_st[1:2]))) && anisotropy) {
    cli_warn("Using PC Matern priors; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }

  df <- if (family$family[1] == "student" && "df" %in% names(family)) family$df else 3

  est_epsilon_model <- 0L
  epsilon_covariate <- rep(0, length(unique(data[[time]])))
  if (!is.null(epsilon_predictor) & !is.null(epsilon_model)) {
    if(epsilon_model %in% c("trend","trend-re")) {
      # covariate vector dimensioned by number of time steps
      time_steps <- unique(data[[time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- data[data[[time]] == time_steps[i],
                                     epsilon_predictor, drop = TRUE][[1]]
      }
      est_epsilon_model <- 1L
    }
  }
  # flags for turning off the trend and random effects
  est_epsilon_slope <- 0
  if(!is.null(epsilon_model)) {
    if(epsilon_model %in% c("trend","trend-re")) {
      est_epsilon_slope <- 1L
      est_epsilon_model <- 1L
    }
  }
  est_epsilon_re <- 0
  if(!is.null(epsilon_model)) {
    if(epsilon_model[1] %in% c("re","trend-re")) {
      est_epsilon_re <- 1L
      est_epsilon_model <- 1L
    }
  }


  priors_b <- priors$b
  .priors <- priors
  .priors$b <- NULL # removes this in the list, so not passed in as data
  if (nrow(priors_b) == 1L && ncol(X_ij[[1]]) > 1L) { # TODO change hard coded index on X_ij
    if (!is.na(priors_b[[1]])) {
      message("Expanding `b` priors to match model matrix.")
    }
    # creates matrix that is 2 columns of NAs, rows = number of unique bs
    # Instead of passing in a 2-column matrix of NAs, pass in a matrix that
    # has means in first col and the remainder is Var-cov matrix
    priors_b <- mvnormal(rep(NA, ncol(X_ij[[1]])))# TODO change hard coded index on X_ij
  }
  # ncol(X_ij) may occur if time varying model, no intercept
  if (ncol(X_ij[[1]]) > 0 & !identical(nrow(priors_b), ncol(X_ij[[1]])))# TODO change hard coded index on X_ij
    cli_abort("The number of 'b' priors does not match the model matrix.")
  if (ncol(priors_b) == 2 && attributes(priors_b)$dist == "normal") {
    # normal priors passed in by user; change to MVN diagonal matrix
    # as.numeric() here prevents diag(NA) taking on factor
    if (length(priors_b[, 2]) == 1L) {
      if (is.na(priors_b[, 2])) priors_b[, 2] <- 1
    }
    priors_b <- mvnormal(location = priors_b[, 1], scale = diag(as.numeric(priors_b[, 2]), ncol = nrow(priors_b)))
  }
  # in some cases, priors_b will be a mix of NAs (no prior) and numeric values
  # easiest way to deal with this is to subset the Sigma matrix
  not_na <- which(!is.na(priors_b[,1]))
  if (length(not_na) == 0L) {
    priors_b_Sigma <- diag(2) # TMB can't handle 0-dim matrix
  } else {
    Sigma <- as.matrix(priors_b[, -1])
    priors_b_Sigma <- as.matrix(Sigma[not_na, not_na])
  }

  .priors$threshold <- NULL # passed in separately
  priors_threshold <- priors$threshold
  if (all(is.na(priors_threshold)))
    priors_threshold <- matrix(NA, nrow = thresh[[1]]$n_par, ncol = 2L)
  msg <- paste0(
    "Number of threshold priors does not match the number of threshold parameters (",
    thresh[[1]]$n_par, ")")
  assert_that(nrow(priors_threshold) == thresh[[1]]$n_par, msg = msg)

  if (!"A_st" %in% names(spde)) cli_abort("`mesh` was created with an old version of `make_mesh()`.")
  if (delta) y_i <- cbind(ifelse(y_i > 0, 1, 0), ifelse(y_i > 0, y_i, NA_real_))
  if (!delta) y_i <- matrix(y_i, ncol = 1L)

  # TODO: make this cleaner
  X_ij_list <- list()
  #X_ij_array <- array(data = NA, dim = c(nrow(X_ij[[1]]), ncol(X_ij[[1]]), n_m))
  for (i in seq_len(n_m)) X_ij_list[[i]] <- X_ij[[i]]

  n_t <- length(unique(data[[time]]))
  tmb_data <- list(
    y_i        = y_i,
    n_t        = n_t,
    z_i        = z_i,
    offset_i   = offset,
    proj_offset_i = 0,
    A_st       = spde$A_st,
    sim_re     = if ("sim_re" %in% names(experimental)) as.integer(experimental$sim_re) else rep(0L, 6),
    A_spatial_index = spde$sdm_spatial_id - 1L,
    year_i     = make_year_i(data[[time]]),
    ar1_fields = ar1_fields,
    simulate_t = rep(1L, n_t),
    rw_fields =  rw_fields,
    X_ij       = X_ij_list,
    X_rw_ik    = X_rw_ik,
    Zs         = sm$Zs, # optional smoother basis function matrices
    Xs         = sm$Xs, # optional smoother linear effect matrix
    proj_Zs    = list(),
    proj_Xs    = matrix(nrow = 0L, ncol = 0L),
    b_smooth_start = sm$b_smooth_start,
    proj_lon   = 0,
    proj_lat   = 0,
    do_predict = 0L,
    calc_se    = 0L,
    pop_pred   = 0L,
    short_newdata = 0L,
    exclude_RE = rep(0L, ncol(RE_indexes)),
    weights_i  = if (!is.null(weights)) weights else rep(1, length(y_i)),
    area_i     = rep(1, length(y_i)),
    normalize_in_r = 0L, # not used first time
    flag = 1L, # part of TMB::normalize()
    calc_index_totals = 0L,
    calc_cog = 0L,
    random_walk = as.integer(!is.null(time_varying) && time_varying_type == "rw"),
    ar1_time = as.integer(!is.null(time_varying) && time_varying_type == "ar1"),
    priors_threshold = priors_threshold,
    priors_b_n = length(not_na),
    priors_b_index = not_na - 1L,
    priors_b_mean = priors_b[not_na,1],
    priors_b_Sigma = priors_b_Sigma,
    priors = as.numeric(unlist(.priors)),
    share_range = as.integer(share_range),
    include_spatial = as.integer(include_spatial),
    proj_mesh  = Matrix::Matrix(c(0,0,2:0), 3, 5), # dummy
    proj_X_ij  = list(matrix(0, ncol = 1, nrow = 1)), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year  = 0, # dummy
    proj_spatial_index = 0, # dummy
    proj_z_i = matrix(0, nrow = 1, ncol = n_m), # dummy
    spde_aniso = make_anisotropy_spde(spde, anisotropy),
    spde       = spde$spde$param.inla[c("M0","M1","M2")],
    barrier = as.integer(barrier),
    spde_barrier = make_barrier_spde(spde),
    barrier_scaling = if (barrier) spde$barrier_scaling else c(1, 1),
    anisotropy = as.integer(anisotropy),
    family     = .valid_family[family$family],
    size = c(size),
    link       = .valid_link[family$link],
    df         = df,
    spatial_only = as.integer(spatial_only),
    spatial_covariate = as.integer(!is.null(spatial_varying)),
    calc_quadratic_range = as.integer(quadratic_roots),
    X_threshold = thresh[[1]]$X_threshold, # TODO: don't hardcode index thresh[[1]]
    mi_est = as.integer(thresh[[1]]$mi_est),
    proj_X_threshold = 0, # dummy
    threshold_func = thresh[[1]]$threshold_func,# TODO: don't hardcode index thresh[[1]]
    RE_indexes = RE_indexes,
    proj_RE_indexes = matrix(0, ncol = 0, nrow = 1), # dummy
    nobs_RE = nobs_RE,
    ln_tau_G_index = ln_tau_G_index,
    est_epsilon_model = as.integer(est_epsilon_model),
    epsilon_predictor = epsilon_covariate,
    est_epsilon_slope = as.integer(est_epsilon_slope),
    est_epsilon_re = as.integer(est_epsilon_re),
    has_smooths = as.integer(sm$has_smooths),
    upr = upr,
    lwr = lwr,
    poisson_link_delta = as.integer(isTRUE(family$type == "poisson_link_delta")),
    stan_flag = as.integer(bayesian),
    no_spatial = no_spatial
  )

  b_thresh <- matrix(0, 2L, n_m)
  if (thresh[[1]]$threshold_func == 2L) b_thresh <- matrix(0, 3L, n_m) # logistic #TODO: change hard coding on index of thresh[[1]]

  tmb_params <- list(
    ln_H_input = matrix(0, nrow = 2L, ncol = n_m),
    b_j        = rep(0, ncol(X_ij[[1]])), # TODO: verify ok
    b_j2       = if (delta) rep(0, ncol(X_ij[[2]])) else numeric(0), # TODO: verify ok
    bs         = if (sm$has_smooths) matrix(0, nrow = ncol(sm$Xs), ncol = n_m) else array(0),
    ln_tau_O   = rep(0, n_m),
    ln_tau_Z = matrix(0, n_z, n_m),
    ln_tau_E   = rep(0, n_m),
    ln_kappa   = matrix(0, 2L, n_m),
    # ln_kappa   = rep(log(sqrt(8) / median(stats::dist(spde$mesh$loc))), 2),
    thetaf     = 0,
    logit_p_mix = 0,
    log_ratio_mix = 0,
    ln_phi     = rep(0, n_m),
    ln_tau_V   = matrix(0, ncol(X_rw_ik), n_m),
    rho_time_unscaled = matrix(0, ncol(X_rw_ik), n_m),
    ar1_phi    = rep(0, n_m),
    ln_tau_G   = matrix(0, ncol(RE_indexes), n_m),
    RE         = matrix(0, sum(nobs_RE), n_m),
    b_rw_t     = array(0, dim = c(tmb_data$n_t, ncol(X_rw_ik), n_m)),
    omega_s    = matrix(0, n_s, n_m),
    zeta_s    = array(0, dim = c(n_s, n_z, n_m)),
    epsilon_st = array(0, dim = c(n_s, tmb_data$n_t, n_m)),
    b_threshold = matrix(0, thresh[[1]]$n_par, n_m),
    b_epsilon = rep(0, n_m),
    ln_epsilon_re_sigma = rep(0, n_m),
    epsilon_re = matrix(0, tmb_data$n_t, n_m),
    b_smooth = if (sm$has_smooths) matrix(0, sum(sm$sm_dims), n_m) else array(0),
    ln_smooth_sigma = if (sm$has_smooths) matrix(0, length(sm$sm_dims), n_m) else array(0)
  )
  if (identical(family$link, "inverse") && family$family[1] %in% c("Gamma", "gaussian", "student") && !delta) {
    fam <- family
    if (family$family == "student") fam$family <- "gaussian"
    temp <- mgcv::gam(formula = formula[[1]], data = data, family = fam)
    tmb_params$b_j <- stats::coef(temp)
  }

  #if (delta && !is.null(thresh$threshold_parameter)) cli_abort("Thresholds not implemented with delta models yet.") # TODO DELTA

  # Map off parameters not needed
  tmb_map <- map_all_params(tmb_params)
  tmb_map$b_j <- NULL
  if (delta) tmb_map$b_j2 <- NULL
  if (family$family[[1]] == "tweedie") tmb_map$thetaf <- NULL
  if (family$family[[1]] %in% c("gamma_mix", "lognormal_mix")) {
    tmb_map$log_ratio_mix <- NULL
    tmb_map$logit_p_mix <- NULL
  }
  if (delta) {
    if(family$family[[2]] %in% c("gamma_mix","lognormal_mix")) {
      tmb_map$log_ratio_mix <- NULL
      tmb_map$logit_p_mix <- NULL
    }
  }
  tmb_map$ln_phi <- rep(1, n_m)
  if (family$family[[1]] %in% c("binomial", "poisson", "censored_poisson"))
    tmb_map$ln_phi[1] <- factor(NA)
  if (delta) {
    if (family$family[[2]] %in% c("binomial", "poisson", "censored_poisson"))
      tmb_map$ln_phi[2] <- factor(NA)
    else
      tmb_map$ln_phi[2] <- 2
  }
  tmb_map$ln_phi <- as.factor(tmb_map$ln_phi)
  if (!is.null(thresh[[1]]$threshold_parameter)) tmb_map$b_threshold <- NULL

  # optional models on spatiotemporal sd parameter
  # if (est_epsilon_re == 0L) {
  #   tmb_map <- c(tmb_map,
  #                list(
  #                  ln_epsilon_re_sigma = factor(rep(NA, n_m)),
  #                  epsilon_re = factor(rep(NA, tmb_data$n_t))
  #                ))
  # }
  if (est_epsilon_re == 1L) {
    tmb_map <- unmap(tmb_map, c("ln_epsilon_re_sigma","epsilon_re"))
  }
  if (est_epsilon_slope == 1L) {
     tmb_map <- unmap(tmb_map, "b_epsilon")
  }

  if (multiphase && is.null(previous_fit) && do_fit) {

    original_tmb_data <- tmb_data
    # # much faster on first phase!?
    tmb_data$no_spatial <- 1L
    tmb_data$include_spatial <- 0L
    # tmb_data$spatial_only <- rep(1L, length(tmb_data$spatial_only))

    for (i in seq_along(start)) {
      tmb_params[[names(start)[i]]] <- start[[i]]
    }
    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params,
      map = tmb_map, DLL = "sdmTMB", silent = silent)
    lim <- set_limits(tmb_obj1, lower = lower, upper = upper, silent = TRUE)

    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn,
      lower = lim$lower, upper = lim$upper,
      gradient = tmb_obj1$gr, control = .control)

    tmb_data <- original_tmb_data # restore
    # Set starting values based on phase 1:
    tmb_params <- tmb_obj1$env$parList()
    # tmb_data$no_spatial <- FALSE
    # often causes optimization problems if set from phase 1!?
    tmb_params$b_threshold <- matrix(0, thresh[[1]]$n_par, n_m)
  }

  tmb_random <- c()
  if (any(spatial == "on")) {
    tmb_random <- c(tmb_random, "omega_s")
    tmb_map <- unmap(tmb_map, c("omega_s", "ln_tau_O"))
  }
  if (!all(spatiotemporal == "off")) {
    tmb_random <- c(tmb_random, "epsilon_st")
    tmb_map <- unmap(tmb_map, c("ln_tau_E", "epsilon_st"))
  }
  if (!is.null(spatial_varying)) {
    tmb_random <- c(tmb_random, "zeta_s")
    tmb_map <- unmap(tmb_map, c("zeta_s", "ln_tau_Z"))
  }
  if (anisotropy) tmb_map <- unmap(tmb_map, "ln_H_input")
  if (!is.null(time_varying)) {
    tmb_random <- c(tmb_random, "b_rw_t")
    tmb_map <- unmap(tmb_map, c("b_rw_t", "ln_tau_V"))
    if (time_varying_type == "ar1")
      tmb_map <- unmap(tmb_map, "rho_time_unscaled")
  }
  if (est_epsilon_re) {
    tmb_random <- c(tmb_random, "epsilon_re")
    tmb_map <- unmap(tmb_map, c("epsilon_re"))
  }

  tmb_map$ar1_phi <- as.numeric(tmb_map$ar1_phi) # strip factors
  for (i in seq_along(spatiotemporal)) {
    if (spatiotemporal[i] == "ar1") tmb_map$ar1_phi[i] <- i
  }
  tmb_map$ar1_phi <- as.factor(as.integer(as.factor(tmb_map$ar1_phi)))

  if (nobs_RE[[1]] > 0) {
    tmb_random <- c(tmb_random, "RE")
    tmb_map <- unmap(tmb_map, c("ln_tau_G", "RE"))
  }
  if (reml) tmb_random <- c(tmb_random, "b_j")
  if (reml && delta) tmb_random <- c(tmb_random, "b_j2")

##  if (est_epsilon_model >= 2) {
##    # model 2 = re model, model 3 = loglinear-re
##    tmb_random <- c(tmb_random, "epsilon_rw")
##  }

  if (sm$has_smooths) {
    tmb_random <- c(tmb_random, "b_smooth") # smooth random effects
    tmb_map <- unmap(tmb_map, c("b_smooth", "ln_smooth_sigma", "bs"))
  }

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }

  tmb_data$normalize_in_r <- as.integer(normalize)

  if (!is.null(previous_fit)) tmb_map <- previous_fit$tmb_map

  # kappa mapping:
  tmb_map$ln_kappa <- matrix(seq_len(length(tmb_params$ln_kappa)),
    nrow(tmb_params$ln_kappa), ncol(tmb_params$ln_kappa))
  for (m in seq_len(n_m)) {
    if (share_range[m]) tmb_map$ln_kappa[,m] <- tmb_map$ln_kappa[1,m]
    if (spatiotemporal[m] == "off" && spatial[m] == "off") tmb_map$ln_kappa[,m] <- NA
    if (spatiotemporal[m] == "off" && spatial[m] == "on") tmb_map$ln_kappa[,m] <- m
  }
  tmb_map$ln_kappa <- as.factor(as.integer(as.factor(tmb_map$ln_kappa)))

  for (i in seq_along(start)) {
    cli_inform(c(i = paste0("Initiating `", names(start)[i],
      "` at specified starting value(s) of:"),
      paste0("  ", paste(round(start[[i]], 3), collapse = ", "))))
    tmb_params[[names(start)[i]]] <- start[[i]]
  }

  if (!is.matrix(tmb_params[["ln_kappa"]]) && "ln_kappa" %in% names(start)) {
      msg <- c("Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).",
        "It should be the same value in each row if `share_range = TRUE`.")
    cli_abort(msg)
  }
  if (nrow(tmb_params[["ln_kappa"]]) != 2L && "ln_kappa" %in% names(start)) {
      msg <- c("Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).",
        "It should be the same value in each row if `share_range = TRUE`.")
    cli_abort(msg)
  }

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL
  data$sdmTMB_X_ <- data$sdmTMB_Y_ <- NULL

  # delta spatiotemporal mapping:
  if (delta && "off" %in% spatiotemporal) {
    tmb_map$epsilon_st <- array(
      seq_len(length(tmb_params$epsilon_st)),
      dim = dim(tmb_params$epsilon_st)
    )
    tmb_map$ln_tau_E <- seq_len(length(tmb_params$ln_tau_E))
    for (i in which(spatiotemporal == "off")) {
      tmb_map$epsilon_st[,,i] <- NA
      tmb_map$ln_tau_E[i] <- NA
    }
    tmb_map$epsilon_st <- as.factor(tmb_map$epsilon_st)
    tmb_map$ln_tau_E <- as.factor(tmb_map$ln_tau_E)
  }
  # delta spatial mapping:
  if (delta && "off" %in% spatial) {
    tmb_map$omega_s <- array(
      seq_len(length(tmb_params$omega_s)),
      dim = dim(tmb_params$omega_s)
    )
    tmb_map$ln_tau_O <- seq_len(length(tmb_params$ln_tau_O))
    for (i in which(spatial == "off")) {
      tmb_map$omega_s[,i] <- NA
      tmb_map$ln_tau_O[i] <- NA
    }
    tmb_map$omega_s <- as.factor(tmb_map$omega_s)
    tmb_map$ln_tau_O <- as.factor(tmb_map$ln_tau_O)
  }

  if (anisotropy && delta && !"ln_H_input" %in% map) {
    tmb_map$ln_H_input <- factor(c(1, 2, 1, 2)) # share anistropy as in VAST
  }

  if (tmb_data$threshold_func > 0) tmb_map$b_threshold <- NULL

  if (control$profile && delta)
    cli_abort("Profile not yet working with delta models.")

  for (i in seq_along(map)) { # user supplied
    cli_inform(c(i = paste0("Fixing (mapping) `", names(map)[i],
      "` at specified starting value(s) of:"),
      paste0("  ",
        paste(round(tmb_params[[names(map)[i]]], 3), collapse = ", "))))
  }
  tmb_map <- c(map, tmb_map) # add user-supplied mapping

  prof <- c("b_j")
  if (delta) prof <- c(prof, "b_j2")



  out_structure <- structure(list(
    data       = data,
    spde       = spde,
    mesh       = spde,
    formula    = original_formula,
    split_formula = split_formula,
    time_varying = time_varying,
    threshold_parameter = thresh[[1]]$threshold_parameter,
    threshold_function = thresh[[1]]$threshold_func,
    epsilon_predictor = epsilon_predictor,
    time       = time,
    family     = family,
    smoothers = sm,
    response   = y_i,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    spatial_varying = spatial_varying,
    spatial = spatial,
    spatiotemporal = spatiotemporal,
    spatial_varying_formula = spatial_varying_formula,
    reml       = reml,
    priors     = priors,
    nlminb_control = .control,
    control  = control,
    contrasts  = lapply(X_ij, attr, which = "contrasts"),
    terms  = lapply(mf, attr, which = "terms"),
    extra_time = extra_time,
    xlevels    = lapply(seq_along(mf), function(i) stats::.getXlevels(mt[[i]], mf[[i]])),
    call       = match.call(expand.dots = TRUE),
    version    = utils::packageVersion("sdmTMB")),
    class      = "sdmTMB")

  if (do_index) {
    args <- list(object = out_structure, return_tmb_data = TRUE)
    args <- c(args, predict_args)
    tmb_data <- do.call(predict.sdmTMB, args)
    if (!"newdata" %in% names(predict_args)) {
      cli_warn("`newdata` must be supplied if `do_index = TRUE`.")
    }
    if ("bias_correct" %in% names(index_args)) {
      cli_warn("`bias_correct` must be done later with `get_index(..., bias_correct = TRUE)`.")
      index_args$bias_correct <- NULL
    }
    if (!"area" %in% names(index_args)) {
      cli_warn("`area` not supplied to `index_args` but `do_index = TRUE`. Using `area = 1`.")
      if (is.null(index_args)) index_args <- list()
      index_args[["area"]] <- 1
    }
    if (length(index_args$area) == 1L) {
      tmb_data$area_i <- rep(index_args[["area"]], nrow(predict_args[["newdata"]]))
    } else {
      if (length(index_args$area) != nrow(predict_args[["newdata"]]))
        cli_abort("`area` length does not match `nrow(newdata)`.")
      tmb_data$area_i <- index_args[["area"]]
    }
    tmb_data$calc_index_totals <- 1L
    tmb_params[["eps_index"]] <- numeric(0) # for bias correction
    out_structure$do_index <- TRUE
  } else {
    out_structure$do_index <- FALSE
  }

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    profile = if (control$profile) prof else NULL,
    random = tmb_random, DLL = "sdmTMB", silent = silent)
  lim <- set_limits(tmb_obj, lower = lower, upper = upper,
    loc = spde$mesh$loc, silent = FALSE)

  out_structure$tmb_obj <- tmb_obj
  out_structure$tmb_obj <- tmb_obj
  out_structure$tmb_data <- tmb_data
  out_structure$tmb_params <- tmb_params
  out_structure$lower <- lim$lower
  out_structure$upper <- lim$upper

  if (!do_fit) return(out_structure)

  if (normalize) tmb_obj <- TMB::normalize(tmb_obj, flag = "flag", value = 0)

  if (length(tmb_obj$par)) {
    tmb_opt <- stats::nlminb(
      start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
      lower = lim$lower, upper = lim$upper, control = .control)
  } else {
    tmb_opt <- list(par = tmb_obj$par, objective = tmb_obj$fn(tmb_obj$par))
  }

  if (nlminb_loops > 1) {
    if (!silent) cli_inform("running extra nlminb loops\n")
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- tmb_opt[c("iterations", "evaluations")]
      tmb_opt <- stats::nlminb(
        start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
        control = .control, lower = lim$lower, upper = lim$upper)
      tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
      tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  if (newton_loops > 0) {
    if (!silent) cli_inform("running newtonsteps\n")
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(tmb_obj$gr(tmb_opt$par))
      h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
      tmb_opt$par <- tmb_opt$par - solve(h, g)
      tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
    }
  }
  check_bounds(tmb_opt$par, lim$lower, lim$upper)

  sd_report <- TMB::sdreport(tmb_obj, getJointPrecision = get_joint_precision)
  conv <- get_convergence_diagnostics(sd_report)

  out_structure$tmb_obj <- tmb_obj
  out <- c(out_structure, list(
    model      = tmb_opt,
    sd_report  = sd_report,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    pos_def_hessian = sd_report$pdHess))
  `class<-`(out, "sdmTMB")
}

check_bounds <- function(.par, lower, upper) {
  for (i in seq_along(.par)) {
    if (is.finite(lower[i]) || is.finite(upper[i])) {
      if (abs(.par[i] - lower[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its lower bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE)
      }
      if (abs(.par[i] - upper[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its upper bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE)
      }
    }
  }
}

set_limits <- function(tmb_obj, lower, upper, loc = NULL, silent = TRUE) {
  .lower <- stats::setNames(rep(-Inf, length(tmb_obj$par)), names(tmb_obj$par))
  .upper <- stats::setNames(rep(Inf, length(tmb_obj$par)), names(tmb_obj$par))
  for (i_name in names(lower)) {
    if (i_name %in% names(.lower)) {
      .lower[names(.lower) %in% i_name] <- lower[[i_name]]
      if (!silent) message("Setting lower limit for ", i_name, " to ",
        lower[[i_name]], ".")
    }
    if (i_name %in% names(.upper)) {
      .upper[names(.upper) %in% i_name] <- upper[[i_name]]
      if (!silent) message("Setting upper limit for ", i_name, " to ",
        upper[[i_name]], ".")
    }
  }
  if ("ar1_phi" %in% names(tmb_obj$par) &&
      !"ar1_phi" %in% union(names(lower), names(upper))) {
    .lower["ar1_phi"] <- stats::qlogis((-0.999 + 1) / 2)
    .upper["ar1_phi"] <- stats::qlogis((0.999 + 1) / 2)
  }

  # if (!is.null(loc) && !"ln_kappa" %in% union(names(lower), names(upper))) {
  #   .dist <- stats::dist(loc)
  #   range_limits <- c(min(.dist) * 0.25, max(.dist) * 4)
  #   .upper[names(.upper) == "ln_kappa"] <- log(sqrt(8) / range_limits[1])
  #   .lower[names(.upper) == "ln_kappa"] <- log(sqrt(8) / range_limits[2])
  #   message("Setting limits on range to ",
  #     round(range_limits[1], 1), " and ",
  #     round(range_limits[2], 1), ":\n",
  #     "half the minimum and twice the maximum knot distance."
  #   )
  #   if (.upper[names(.upper) == "ln_kappa"][1] <= tmb_obj$par[["ln_kappa"]][1]) {
  #     .upper[names(.upper) == "ln_kappa"] <- tmb_obj$par[["ln_kappa"]][1] + 0.1
  #   }
  #   if (.lower[names(.lower) == "ln_kappa"][1] >= tmb_obj$par[["ln_kappa"]][1]) {
  #     .lower[names(.lower) == "ln_kappa"] <- tmb_obj$par[["ln_kappa"]][1] - 0.1
  #   }
  # }
  list(lower = .lower, upper = .upper)
}

check_spatiotemporal_arg <- function(x, time, .which = 1) {
  sp_len <- length(x)
  .spatiotemporal <- tolower(as.character(x[[.which]]))
  assert_that(.spatiotemporal %in% c("iid", "ar1", "rw", "off", "true", "false"),
    msg = "`spatiotemporal` argument value not valid")
  if (.spatiotemporal == "true") .spatiotemporal <- "iid"
  if (.spatiotemporal == "false") .spatiotemporal <- "off"
  .spatiotemporal <- match.arg(tolower(.spatiotemporal), choices = c("iid", "ar1", "rw", "off"))
  if (is.null(time) && .spatiotemporal != "off" && sp_len >= 1) {
    cli_abort("`spatiotemporal` is set but the `time` argument is missing.")
  }
  .spatiotemporal
}

map_all_params <- function(x) {
  m <- list()
  nm <- names(x)
  for (i in seq_along(x)) {
    m[[i]]<- factor(rep(NA, length(x[[i]])))
  }
 stats::setNames(m, names(x))
}

parse_spatial_arg <- function(spatial) {
  if (!is.logical(spatial[[1]])) {
    spatial <- match.arg(tolower(spatial), choices = c("on", "off"))
  }
  if (identical(spatial, "on") || isTRUE(spatial)) {
    spatial <- "on"
  } else {
    spatial <- "off"
  }
  spatial
}

check_irregalar_time <- function(data, time, spatiotemporal, time_varying) {
  if (any(spatiotemporal %in% c("ar1", "rw")) || !is.null(time_varying)) {
    ti <- sort(unique(data[[time]]))
    if (length(unique(diff(ti))) > 1L) {
      missed <- find_missing_time(data[[time]])
      msg <- c(
        "Detected irregular time spacing with an AR(1) or random walk process.",
        "Consider filling in the missing time slices with `extra_time`.",
        if (length(missed)) {
          paste0("`extra_time = c(", paste(missed, collapse = ", "), ")`")
        }
      )
      cli_inform(msg)
    }
  }
}

find_missing_time <- function(x) {
  if (!is.factor(x)) {
    ti <- sort(unique(x))
    mindiff <- min(diff(ti))
    allx <- seq(min(ti), max(ti), by = mindiff)
    setdiff(allx, ti)
  }
}

unmap <- function(x, v) {
  for (i in v) x[[i]] <- NULL
  x
}
