#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal generalized linear mixed effects model (GLMM)
#' with the TMB (Template Model Builder) R package. Spatial and spatiotemporal
#' random fields are approximated using the SPDE (stochastic partial differential
#' equation) approach, which allows for efficient modeling of data that are
#' correlated in space and/or time. See the [model description
#' vignette](https://sdmTMB.github.io/sdmTMB/articles/model-description.html) for
#' details.
#'
#' @param formula Model formula. IID random intercepts and slopes are possible using
#'   \pkg{lme4} syntax, e.g., `+ (1 | g)` or `+ (0 + depth | g)` or `+ (1 + depth | g)` where `g` is a column of class
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
#'   \code{\link[sdmTMB:families]{betabinomial()}},
#'   \code{\link[sdmTMB:families]{nbinom2()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom2()}},
#'   \code{\link[sdmTMB:families]{nbinom1()}},
#'   \code{\link[sdmTMB:families]{truncated_nbinom1()}},
#'   \code{\link[sdmTMB:families]{censored_poisson()}},
#'   \code{\link[sdmTMB:families]{gamma_mix()}},
#'   \code{\link[sdmTMB:families]{lognormal_mix()}},
#'   \code{\link[sdmTMB:families]{student()}},
#'   \code{\link[sdmTMB:families]{tweedie()}}, and
#'   \code{\link[sdmTMB:families]{gengamma()}}.
#'   Delta/hurdle models (for zero-inflated data) include:
#'   \code{\link[sdmTMB:families]{delta_beta()}},
#'   \code{\link[sdmTMB:families]{delta_gamma()}},
#'   \code{\link[sdmTMB:families]{delta_gamma_mix()}},
#'   \code{\link[sdmTMB:families]{delta_lognormal_mix()}},
#'   \code{\link[sdmTMB:families]{delta_lognormal()}}, and
#'   \code{\link[sdmTMB:families]{delta_truncated_nbinom2()}}.
#'   See the [delta-model
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html) for
#'   details. For binomial family options, see 'Binomial families' in the Details
#'   section below.
#' @param spatial Estimate spatial random fields? Options are `'on'` / `'off'`
#'   or `TRUE` / `FALSE`. Optionally, a list for delta models, e.g. `list('on',
#'   'off')`.
#' @param spatiotemporal Estimate the spatiotemporal random fields as `'iid'`
#'   (independent and identically distributed; default), stationary `'ar1'`
#'   (first-order autoregressive), a random walk (`'rw'`), or fixed at 0
#'   `'off'`. Will be set to `'off'` if `time = NULL`. If a delta model, can be
#'   a list. E.g., `list('off', 'ar1')`. Guidance: Use `'iid'` if temporal
#'   correlation is negligible or already accounted for in fixed effects; `'ar1'`
#'   if correlation between consecutive time steps decays gradually; `'rw'` if
#'   changes between time steps are cumulative (each step builds on the last). If
#'   the AR1 correlation coefficient (rho) is estimated close to 1 (say > 0.99),
#'   consider switching to `'rw'`. See the [model description
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/model-description.html) for
#'   mathematical details. Capitalization is ignored. `TRUE` gets converted to
#'   `'iid'` and `FALSE` gets converted to `'off'`.
#' @param share_range Logical: estimate a shared spatial and spatiotemporal
#'   range parameter (`TRUE`, default) or independent range parameters
#'   (`FALSE`). If a delta model, can be a list. E.g., `list(TRUE, FALSE)`.
#' @param time_varying An optional one-sided formula describing covariates
#'   that should be modelled as a time-varying process. Set the type of
#'   process with `time_varying_type`. See the help for `time_varying_type`
#'   for warnings about modelling the first time step. Structure shared in
#'   delta models.
#' @param time_varying_type Type of time-varying process to apply to
#'   `time_varying` formula. Options: `'rw'` (random walk, default), `'rw0'`
#'   (random walk with mean-zero prior on first time step), or `'ar1'`
#'   (autoregressive, for coefficients that fluctuate around a mean). For `'rw0'`
#'   and `'ar1'`, the coefficient starts at zero in the first time step. For
#'   `'rw'` (default), the first time step is estimated separately—in this case,
#'   avoid including the same covariates in both `formula` and `time_varying` to
#'   prevent non-identifiability (use `~ 0` or `~ -1` in at least one). Structure
#'   shared in delta models.
#' @param spatial_varying An optional one-sided formula of coefficients that
#'   should vary in space as random fields. Allows the effect of a covariate to
#'   differ spatially. You likely want to include the same variable as a fixed
#'   effect in `formula` to estimate the average effect—the spatial field then
#'   represents deviations from that average (since it has mean zero). For example,
#'   use `formula = y ~ depth` and `spatial_varying = ~ 0 + depth` to model an
#'   average depth effect plus spatially varying deviations. If a (scaled) time
#'   column is used, this creates a local-time-trend model. See
#'   \doi{10.1111/ecog.05176} and the [spatial trends
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/spatial-trend-models.html).
#'   Predictors should usually be centered to have mean zero and standard deviation
#'   approximately 1. **The spatial intercept is controlled by the `spatial`
#'   argument**; set `spatial = 'on'` or `'off'` to include or exclude it. For
#'   factor predictors, if `spatial_varying` excludes the intercept (`~ 0` or `~
#'   -1`), set `spatial = 'off'` to match. Structure must be shared in delta
#'   models.
#' @param weights A numeric vector representing optional likelihood weights for
#'   the conditional model. Implemented as in \pkg{glmmTMB}: weights do not have
#'   to sum to one and are not internally modified. Can also be used for trials
#'   with the binomial family; the `weights` argument needs to be a vector and not
#'   a name of the variable in the data frame. See the Details section below.
#' @param offset A numeric vector representing the model offset *or* a character
#'   value representing the column name of the offset. In delta/hurdle models,
#'   this applies only to the positive component. Usually a log transformed
#'   variable.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See the Details
#'   section below.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation
#'   rather than maximum likelihood? REML accounts for uncertainty in estimating
#'   fixed effects and can reduce bias in variance parameter estimates, but
#'   prevents likelihood-based model comparison (e.g., AIC) between models with
#'   different fixed effects. Use `TRUE` if your focus is on random effect
#'   variance parameters; use `FALSE` (default) if comparing models with different
#'   fixed effects or performing index standardization.
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
#'   Most users can ignore this option. The fitted object can be passed directly
#'   to [get_index()].
#' @param predict_args A list of arguments to pass to [predict.sdmTMB()] **if**
#'   `do_index = TRUE`. Most users can ignore this option.
#' @param index_args A list of arguments to pass to [get_index()] **if**
#'   `do_index = TRUE`. Currently, only `area` is supported. Bias correction
#'   can be done when calling [get_index()] on the resulting fitted object.
#'   Most users can ignore this option.
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
#' * `gradients`: marginal log likelihood gradients with respect to each fixed effect
#' * `model`: output from [stats::nlminb()]
#' * `data`: the fitted data
#' * `spde`: the object that was supplied to the `mesh` argument
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
#' sdmTMB fits GLMMs with spatial and/or spatiotemporal random fields, which
#' account for correlation in the data due to spatial proximity, or alternatively,
#' latent spatial and spatiotemporal effects. Spatial fields represent
#' consistent spatial patterns, while
#' spatiotemporal fields represent spatial patterns that vary over time. See the [model
#' description](https://sdmTMB.github.io/sdmTMB/articles/model-description.html)
#' vignette for mathematical details and the paper: \doi{10.18637/jss.v115.i02}
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
#' Within these smooths, the same syntax commonly used in [mgcv::s()] or
#' [mgcv::t2()] can be applied, e.g. 2-dimensional smooths may be constructed
#' with `+ s(x, y)` or `+ t2(x, y)`; smooths can be specific to various factor
#' levels, `+ s(x, by = group)`; the basis function dimensions may be specified,
#' e.g. `+ s(x, k = 4)`; and various types of splines may be constructed such as
#' cyclic splines to model seasonality (perhaps with the `knots` argument also
#' be supplied).
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
#' 0.5 or 0.95. See the
#' [threshold vignette](https://sdmTMB.github.io/sdmTMB/articles/threshold-models.html).
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
#' options `time_varying`, `spatiotemporal = "rw"`, `spatiotemporal = "ar1"`,
#' or a smoother on the time column provide mechanisms to predict over missing
#' time slices with process error.
#'
#' `extra_time` can also be used to fill in missing time steps for the purposes
#' of a random walk or AR(1) process if the gaps between time steps are uneven.
#'
#' `extra_time` can include only extra time steps or all time steps including
#' those found in the fitted data. This latter option may be simpler.
#'
#' **Regularization and priors**
#'
#' You can achieve regularization via penalties (priors) on the fixed effect
#' parameters. See [sdmTMBpriors()]. You can fit the model once without
#' penalties and look at the output of `print(your_model)` or `tidy(your_model)`
#' or fit the model with `do_fit = FALSE` and inspect
#' `head(your_model$tmb_data$X_ij[[1]])` if you want to see how the formula is
#' translated to the fixed effect model matrix. Also see the
#' [Bayesian vignette](https://sdmTMB.github.io/sdmTMB/articles/bayesian.html).
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
#' specifying two formulas as a list: the two formulas cannot have smoothers or
#' threshold effects. For now, these must be specified
#' through a single formula that is shared across the two models.
#'
#' The main advantage of specifying such models using a delta family (compared
#' to fitting two separate models) is (1) coding simplicity and (2) calculation
#' of uncertainty on derived quantities such as an index of abundance with
#' [get_index()] using the generalized delta method within TMB. Also, selected
#' parameters can be shared across the models.
#'
#' See the [delta-model vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html).
#'
#' **Index standardization**
#'
#' For index standardization, you may wish to include `0 + as.factor(year)`
#' (or whatever the time column is called) in the formula. See a basic
#' example of index standardization in the relevant
#' [package vignette](https://sdmTMB.github.io/sdmTMB/articles/index-standardization.html).
#' You will need to specify the `time` argument. See [get_index()].
#'
#' @references
#'
#' **Main reference introducing the package to cite when using sdmTMB:**
#'
#' Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett., J.T. Thorson. 2025.
#' sdmTMB: an R package for fast, flexible, and user-friendly generalized linear
#' mixed effects models with spatial and spatiotemporal random fields.
#' Journal of Statistical Software. 115(2):1--46. \doi{10.18637/jss.v115.i02}.
#'
#' *Reference for local trends:*
#'
#' Barnett, L.A.K., E.J. Ward, S.C. Anderson. 2021. Improving estimates of species
#' distribution change by incorporating local trends. Ecography. 44(3):427-439.
#' \doi{10.1111/ecog.05176}.
#'
#' *Further explanation of the model and application to calculating climate
#' velocities:*
#'
#' English, P., E.J. Ward, C.N. Rooper, R.E. Forrest, L.A. Rogers, K.L. Hunter,
#' A.M. Edwards, B.M. Connors, S.C. Anderson. 2021. Contrasting climate velocity
#' impacts in warm and cool locations show that effects of marine warming are
#' worse in already warmer temperate waters. Fish and Fisheries. 23(1) 239-255.
#' \doi{10.1111/faf.12613}.
#'
#' *Discussion of and illustration of some decision points when fitting these
#' models:*
#'
#' Commander, C.J.C., L.A.K. Barnett, E.J. Ward, S.C. Anderson, T.E.
#' Essington. 2022. The shadow model: how and why small choices in
#' spatially explicit species distribution models affect predictions. PeerJ 10:
#' e12783. \doi{10.7717/peerj.12783}.
#'
#' *Application and description of threshold/break-point models:*
#'
#' Essington, T.E., S.C. Anderson, L.A.K. Barnett, H.M. Berger, S.A. Siedlecki,
#' E.J. Ward. 2022. Advancing statistical models to reveal the effect of
#' dissolved oxygen on the spatial distribution of marine taxa using thresholds
#' and a physiologically based index. Ecography. 2022: e06249
#' \doi{10.1111/ecog.06249}.
#'
#' *Application to fish body condition:*
#'
#' Lindmark, M., S.C. Anderson, M. Gogina, M. Casini. Evaluating drivers of
#' spatiotemporal individual condition of a bottom-associated marine fish.
#' bioRxiv 2022.04.19.488709. \doi{10.1101/2022.04.19.488709}.
#'
#' *Several sections of the original TMB model code were adapted from the
#' VAST R package:*
#'
#' Thorson, J.T. 2019. Guidance for decisions using the Vector Autoregressive
#' Spatio-Temporal (VAST) package in stock, ecosystem, habitat and climate
#' assessments. Fish. Res. 210:143–161.
#' \doi{10.1016/j.fishres.2018.10.013}.
#'
#' *Code for the `family` R-to-TMB implementation, selected parameterizations of
#' the observation likelihoods, general package structure inspiration, and the
#' idea behind the TMB prediction approach were adapted from the glmmTMB R
#' package:*
#'
#' Brooks, M.E., K. Kristensen, K.J. van Benthem, A. Magnusson,
#' C.W. Berg, A. Nielsen, H.J. Skaug, M. Maechler, B.M. Bolker.
#'  2017. glmmTMB Balances Speed and Flexibility Among Packages for
#' Zero-inflated Generalized Linear Mixed Modeling. The R Journal, 9(2):378-400.
#' \doi{10.32614/rj-2017-066}.
#'
#' *Implementation of geometric anisotropy with the SPDE and use of random
#' field GLMMs for index standardization*:
#'
#' Thorson, J.T., A.O. Shelton, E.J. Ward, H.J. Skaug. 2015.
#' Geostatistical delta-generalized linear mixed models improve precision for
#' estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci.
#' 72(5): 1297–1310. \doi{10.1093/icesjms/fsu243}.
#'
#' @export
#'
#' @examplesIf require("visreg", quietly = TRUE) && require("ggeffects", quietly = TRUE)
#' library(sdmTMB)
#'
#' # Build a mesh to implement the SPDE approach:
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#'
#' # - this example uses a fairly coarse mesh so these examples run quickly
#' # - 'cutoff' is the minimum distance between mesh vertices in units of the
#' #   x and y coordinates
#' # - 'cutoff = 10' might make more sense in applied situations for this dataset
#' # - or build any mesh in 'fmesher' and pass it to the 'mesh' argument in make_mesh()`
#' # - the mesh is not needed if you will be turning off all
#' #   spatial/spatiotemporal random fields
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
#' # Predict on the fitted data; see ?predict.sdmTMB
#' p <- predict(fit)
#'
#' # Predict on new data:
#' p <- predict(fit, newdata = qcs_grid)
#' head(p)
#'
#' \donttest{
#' # Visualize the depth effect with ggeffects:
#' ggeffects::ggpredict(fit, "depth [all]") |> plot()
#'
#' # Visualize depth effect with visreg: (see ?visreg_delta)
#' visreg::visreg(fit, xvar = "depth") # link space; randomized quantile residuals
#' visreg::visreg(fit, xvar = "depth", scale = "response")
#' visreg::visreg(fit, xvar = "depth", scale = "response", gg = TRUE, rug = FALSE)
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
#' # IID random slopes and intercepts (implicit) by year:
#' fit <- sdmTMB(
#'   density ~ (depth | fyear), #<
#'   data = pcod_2011, mesh = mesh,
#'   family = tweedie(link = "log")
#' )
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
#' }
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
    time_varying_type = c("rw", "rw0", "ar1"),
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
    experimental = NULL) {
  data <- droplevels(data) # if data was subset, strips absent factors

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
      function(i) {
        check_spatiotemporal_arg(
          spatiotemporal,
          time = time, .which = i
        )
      },
      FUN.VALUE = character(1L)
    )
  } else {
    if (is.null(time)) {
      spatiotemporal <- rep("off", n_m)
    } else {
      spatiotemporal <- rep("iid", n_m)
    }
  }

  if (is.null(time)) {
    spatial_only <- rep(TRUE, n_m)
  } else {
    spatial_only <- ifelse(spatiotemporal == "off", TRUE, FALSE)
  }

  if (is.list(spatial)) {
    spatial <- vapply(spatial, parse_spatial_arg, FUN.VALUE = character(1L))
  } else {
    spatial <- rep(parse_spatial_arg(spatial), n_m)
  }
  include_spatial <- "on" %in% spatial

  if (!include_spatial && !is.null(spatial_varying)) {
    # move intercept into spatial_varying
    omit_spatial_intercept <- TRUE
    include_spatial <- TRUE
    spatial <- rep("on", length(spatial)) # checked and turned off later if needed
  } else {
    omit_spatial_intercept <- FALSE
  }

  if (!include_spatial && all(spatiotemporal == "off") || !include_spatial && all(spatial_only)) {
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
  share_range[spatial == "off"] <- TRUE

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
  }

  normalize <- control$normalize
  nlminb_loops <- control$nlminb_loops
  newton_loops <- control$newton_loops
  quadratic_roots <- control$quadratic_roots
  start <- control$start
  multiphase <- control$multiphase
  getsd <- control$getsd
  map <- control$map
  lower <- control$lower
  upper <- control$upper
  get_joint_precision <- control$get_joint_precision
  upr <- control$censored_upper
  suppress_nlminb_warnings <- control$suppress_nlminb_warnings
  collapse_spatial_variance <- control$collapse_spatial_variance
  collapse_threshold <- control$collapse_threshold

  dot_checks <- c(
    "lower", "upper", "profile", "parallel", "censored_upper", "getsd",
    "nlminb_loops", "newton_steps", "mgcv", "quadratic_roots", "multiphase",
    "newton_loops", "start", "map", "get_joint_precision", "normalize",
    "suppress_nlminb_warnings", "collapse_spatial_variance", "collapse_threshold"
  )
  .control <- control
  # FIXME; automate this from sdmTMcontrol args?
  for (i in dot_checks) .control[[i]] <- NULL # what's left should be for nlminb

  ar1_fields <- spatiotemporal == "ar1"
  rw_fields <- spatiotemporal == "rw"
  assert_that(
    is.logical(reml), is.logical(anisotropy), is.logical(share_range),
    is.logical(silent), is.logical(multiphase), is.logical(normalize)
  )

  if (!is.null(spatial_varying)) assert_that(class(spatial_varying) %in% c("formula", "list"))
  if (!is.null(time_varying)) assert_that(class(time_varying) %in% c("formula", "list"))
  if (!is.null(previous_fit)) assert_that(identical(class(previous_fit), "sdmTMB"))
  assert_that(is.list(priors))
  assert_that(is.list(.control))
  if (!is.null(time)) assert_that(is.character(time))
  assert_that(inherits(spde, "sdmTMBmesh"))
  assert_that(class(formula) %in% c("formula", "list"))
  assert_that(inherits(data, "data.frame"))
  time_varying_type <- match.arg(time_varying_type)
  # if (!is.null(map) && length(map) != length(start)) {
  #   cli_warn(c("`length(map) != length(start)`.",
  #     "You likely want to specify `start` values if you are setting the `map` argument."))
  # }

  if (!is.null(time)) {
    assert_that(time %in% names(data),
      msg = "Specified `time` column is missing from `data`."
    )
    assert_that(sum(is.na(as.numeric(data[[time]]))) == 0L,
      msg = "Specified `time` column can't be coerced to a numeric value or contains NAs. Please remove any NAs in the time column."
    )
    if (!is.null(extra_time)) { # 320 protect against factor(time), extra_time clash
      if (!is.list(formula)) {
        .form <- list(formula)
      } else {
        .form <- formula
      }
      x <- unlist(lapply(.form, \(x) attr(stats::terms(x), "term.labels")))
      xf <- x[grep("factor\\(", x)]
      if (any(c(grep(time, xf), grep(paste0("^", time, "$"), x)))) {
        cli::cli_warn("Detected potential formula-time column clash. E.g., assuming 'year' is your time column: `formula = ... + factor(year)` combined with `time = 'year'`, and 'extra_time' specified. This can produce a non-identiable model because extra factor levels for the missing time slices will be created. To avoid this, rename your factor time column used in your formula. E.g. create a new column 'year_factor' in your data and use that in the formula. See issue https://github.com/sdmTMB/sdmTMB/issues/320.")
      }
    }
  }
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  } else {
    if (sum(is.na(data[[time]])) > 1) {
      cli_abort("There is at least one NA value in the time column. Please remove it.")
    }
  }
  if (is.factor(data[[time]])) {
    if (length(levels(data[[time]])) > length(unique(data[[time]]))) {
      cli_abort(
        "The time column is a factor and there are extra factor levels. ",
        "Please remove these or turn your time column into an integer."
      )
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

  if (family$family[1] == "censored_poisson") {
    if ("lwr" %in% names(experimental) || "upr" %in% names(experimental)) {
      cli_abort("Detected `lwr` or `upr` in `experimental`. `lwr` is no longer needed and `upr` is now specified as `control = sdmTMBcontrol(censored_upper = ...)`.")
    }
    if (is.null(upr)) cli_abort("`censored_upper` must be defined in `control = sdmTMBcontrol()` to use the censored Poisson distribution.")
    assert_that(length(upr) == nrow(data))
  }
  if (is.null(upr)) upr <- Inf

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
    thresh <- list(
      check_and_parse_thresh_params(formula[[1]], data),
      check_and_parse_thresh_params(formula[[2]], data)
    )
    formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
  }

  if (is.character(offset)) {
    offset <- data[[offset]]
  }

  check_irregalar_time(data, time, spatiotemporal, time_varying, extra_time = extra_time)

  spatial_varying_formula <- spatial_varying # save it
  if (!is.null(spatial_varying)) {
    mf1 <- model.frame(spatial_varying, data)
    for (i in seq_len(ncol(mf1))) {
      if (is.character(mf1[[i]])) {
        cli_warn(paste0(
          "Detected '{colnames(mf1)[i]}' as a character term in the ",
          "'spatial_varying' formula. We suggest you make this a factor if you plan ",
          "to predict with only some factor levels. `as.factor({colnames(mf1)[i]})`."
        ))
      }
    }
    z_i <- model.matrix(spatial_varying, data)
    .int <- sum(grep("(Intercept)", colnames(z_i)) > 0)
    if (length(attr(z_i, "contrasts")) && !.int && !omit_spatial_intercept) { # factors with ~ 0 or ~ -1
      msg <- c(
        "Detected predictors with factor levels in `spatial_varying` with the intercept omitted from the `spatial_varying` formula.",
        "You likely want to set `spatial = 'off'` since the constant spatial field (`omega_s`) also represents a spatial intercept.`"
      )
      cli_inform(paste(msg, collapse = " "))
    }
    .int <- grep("(Intercept)", colnames(z_i))
    if (sum(.int) > 0) z_i <- z_i[, -.int, drop = FALSE]
    spatial_varying <- colnames(z_i)
    svc_contrasts <- attr(z_i, which = "contrasts")
  } else {
    z_i <- matrix(0, nrow(data), 0L)
    svc_contrasts <- NULL
  }
  n_z <- ncol(z_i)

  if (any(grepl("offset\\(", formula))) {
    cli_abort("Detected `offset()` in formula. Offsets in sdmTMB must be specified via the `offset` argument.")
  }
  contains_offset <- check_offset(formula[[1]]) # deprecated check

  split_formula <- list() # passed to out structure, not TMB
  X_ij <- list() # main effects, passed into TMB
  mf <- list()
  mt <- list()
  sm <- list()

  for (ii in seq_along(formula)) {
    contains_offset <- check_offset(formula[[ii]])

    # RE_names is either the names of the groups OR a character w/length 0 if none exist
    # anything in a list here needs to be saved for tmb data

    # smoothers removed, but random effects remain:
    formula_no_sm <- remove_s_and_t2(formula[ii][[1]])
    # random effects parsed into data structures:
    split_formula[[ii]] <- parse_formula(formula_no_sm, data)

    # save formula with no bars (but with smoothers)
    formula_no_bars <- reformulas::nobars(formula[ii][[1]])
    formula_no_bars_no_sm <- remove_s_and_t2(formula_no_bars)
    X_ij[[ii]] <- model.matrix(formula_no_bars_no_sm, data)
    mf[[ii]] <- model.frame(formula_no_bars_no_sm, data)
    vars <- colnames(mf[[ii]])
    for (g in seq_along(vars)) {
      if (any(is.infinite(mf[[ii]][, g, drop = TRUE]))) {
        cli_abort("Column `{vars[g]}` had an Inf/-Inf value. Please remove this before fitting the model.")
      }
    }

    mt[[ii]] <- attr(mf[[ii]], "terms")
    # parse everything mgcv + smoothers:
    sm[[ii]] <- parse_smoothers(formula = formula_no_bars, data = data, knots = knots)
    sm[[ii]]$split_formula <- split_formula
    sm[[ii]]$formula_no_sm <- formula_no_sm
    sm[[ii]]$formula_no_bars <- formula_no_bars
    sm[[ii]]$formula_no_bars_no_sm <- formula_no_bars_no_sm
  }

  # random slopes and intercepts --------------------------------------------
  # bind the elements of split_formula[[ii]] together to pass into TMB
  # add a new column to each dataframe storing the model number (1, 2, ...)
  split_formula_dfs <- add_model_index(split_formula, "re_df")
  re_cov_df <- do.call(rbind, split_formula_dfs) # bind all dataframes together

  re_cov_term_dfs <- add_model_index(split_formula, "re_cov_term_map")
  re_cov_df_map <- do.call(rbind, re_cov_term_dfs)

  re_b_dfs <- add_model_index(split_formula, "re_b_df")
  re_b_df <- do.call(rbind, re_b_dfs)

  re_b_map_dfs <- add_model_index(split_formula, "re_b_map")
  re_b_map <- do.call(rbind, re_b_map_dfs)

  n_re_groups <- rep(0, length(formula)) # integer array of number of groups per model
  Zt_list <- list() # list of sparse matrices for random effects
  n_re_betas <- rep(0, length(formula))
  n_re_cov <- rep(0, length(formula))
  for (ii in seq_len(length(formula))) {
    n_re_groups[ii] <- length(which(re_cov_df_map$model == ii))
    if (n_re_groups[ii] > 0) {
      # transpose of matrix, so cols = observation, row = REs:
      Zt_list[[ii]] <- split_formula[[ii]]$re_cov_terms$Zt
      # number of unique random effect deviations being estimated:
      n_re_betas[ii] <- max(re_b_df$end[which(re_b_df$model == ii)])
      # number of unique random effect variances being estimated:
      n_re_cov[ii] <- max(re_cov_df_map$end[which(re_cov_df_map$model == ii)])
    }
  }
  max_re_betas <- max(n_re_betas) + 1 # most number of estimated betas across all models
  max_re_cov <- max(n_re_cov) + 1

  var_indx_matrix <- matrix(0, nrow = max_re_betas, ncol = length(formula))
  for (i in seq_len(length(formula))) {
    sd_vec <- split_formula[[ii]]$var_indx_vector
    if (n_re_groups[ii] > 0) var_indx_matrix[seq_along(sd_vec), i] <- sd_vec
  }

  if (delta) {
    if (any(unlist(lapply(sm, `[[`, "has_smooths")))) {
      if (original_formula[[1]] != original_formula[[2]]) {
        msg <- paste0(
          "For now, if delta models contain smoothers, both components ",
          "must have the same main-effects formula."
        )
        cli_abort(msg)
      }
    }
  }

  # always shared; only keep track of one:
  sm <- sm[[1]]

  y_i <- model.response(mf[[1]], "numeric")

  # Filter weights and offset to match NA-filtered response
  # model.frame() removes NAs by default; external vectors need same filtering
  na_action <- attr(mf[[1]], "na.action")
  if (!is.null(na_action)) {
    # na.omit creates "omit" class with row indices that were removed
    if (!is.null(weights)) {
      weights <- weights[-na_action]
    }
    if (!is.null(offset)) {
      offset <- offset[-na_action]
    }
  }

  if (delta) {
    y_i2 <- model.response(mf[[2]], "numeric")
    if (!identical(y_i, y_i2)) {
      cli_abort("Response variable should be the same in both parts of the delta formula.")
    }
  }

  if (family$family[1] %in% c("Gamma", "lognormal") && min(y_i) <= 0 && !delta) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }
  if (family$family[1] == "censored_poisson") {
    assert_that(mean(upr - y_i, na.rm = TRUE) >= 0)
  }

  # This is taken from approach in glmmTMB to match how they handle binomial
  # yobs could be a factor -> treat as binary following glm
  # yobs could be cbind(success, failure)
  # yobs could be binary
  # (yobs, weights) could be (proportions, size)
  # On the C++ side 'yobs' must be the number of successes.
  size <- rep(1, nrow(X_ij[[1]])) # for non-binomial case TODO: change hard coded index

  # Helper function for binomial-type response processing
  process_binomial_response <- function(mf, weights) {
    ## call this to catch the factor / matrix cases
    y_i <- model.response(mf[[1]], type = "any")
    size <- rep(1, length(y_i))

    ## allow character
    if (is.character(y_i)) {
      y_i <- model.response(mf[[1]], type = "factor")
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
    }
    if (is.factor(y_i)) {
      if (nlevels(y_i) > 2) {
        cli_abort("More than 2 levels detected for response")
      }
      ## following glm, 'success' is interpreted as the factor not
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
          if (is.null(weights)) {
            cli_abort(c(
              "`weights` argument was not specified but proportions were supplied to the binomial or betabinomial family.",
              "Please supply the `weights` argument with the number of trials per event."))
          }
          y_i <- weights * y_i
          size <- weights
          weights <- rep(1, length(y_i))
        }
      }
    } # https://github.com/sdmTMB/sdmTMB/issues/172
    if (is.logical(y_i)) {
      msg <- paste0(
        "We recommend against using `TRUE`/`FALSE` ",
        "response values if you are going to use the `visreg::visreg()` ",
        "function after. Consider converting to integer with `as.integer()`."
      )
      cli_warn(msg)
    }

    list(y_i = y_i, size = size, weights = weights)
  }

  if ((identical(family$family[1], "binomial") || identical(family$family[1], "betabinomial")) && !delta) {
    result <- process_binomial_response(mf, weights)
    y_i <- result$y_i
    size <- result$size
    weights <- result$weights
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

  # Student-t df: can be NULL (estimate, default) or numeric (fix)
  student_df_fixed <- if (family$family[1] == "student" && "df" %in% names(family)) {
    family$df
  } else {
    NULL
  }
  estimate_student_df <- family$family[1] == "student" && is.null(student_df_fixed)

  est_epsilon_model <- 0L
  epsilon_covariate <- rep(0, length(unique(data[[time]])))
  if (!is.null(epsilon_predictor) & !is.null(epsilon_model)) {
    if (epsilon_model %in% c("trend", "trend-re")) {
      # covariate vector dimensioned by number of time steps
      time_steps <- unique(data[[time]])
      for (i in seq_along(time_steps)) {
        epsilon_covariate[i] <- data[data[[time]] == time_steps[i],
          epsilon_predictor,
          drop = TRUE
        ][[1]]
      }
      est_epsilon_model <- 1L
    }
  }
  # flags for turning off the trend and random effects
  est_epsilon_slope <- 0
  if (!is.null(epsilon_model)) {
    if (epsilon_model %in% c("trend", "trend-re")) {
      est_epsilon_slope <- 1L
      est_epsilon_model <- 1L
    }
  }
  est_epsilon_re <- 0
  if (!is.null(epsilon_model)) {
    if (epsilon_model[1] %in% c("re", "trend-re")) {
      est_epsilon_re <- 1L
      est_epsilon_model <- 1L
    }
  }

  priors_b <- priors$b
  priors_sigma_V <- priors$sigma_V
  .priors <- priors
  .priors$b <- NULL # removes this in the list, so not passed in as data
  .priors$sigma_V <- NULL # removes this in the list, so not passed in as data
  if (nrow(priors_b) == 1L && ncol(X_ij[[1]]) > 1L) { # TODO change hard coded index on X_ij
    if (!is.na(priors_b[[1]])) {
      message("Expanding `b` priors to match model matrix.")
    }
    # creates matrix that is 2 columns of NAs, rows = number of unique bs
    # Instead of passing in a 2-column matrix of NAs, pass in a matrix that
    # has means in first col and the remainder is Var-cov matrix
    priors_b <- mvnormal(rep(NA, ncol(X_ij[[1]]))) # TODO change hard coded index on X_ij
  }
  # ncol(X_ij) may occur if time varying model, no intercept
  if (ncol(X_ij[[1]]) > 0 & !identical(nrow(priors_b), ncol(X_ij[[1]]))) { # TODO change hard coded index on X_ij
    cli_abort("The number of 'b' priors does not match the model matrix.")
  }
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
  not_na <- which(!is.na(priors_b[, 1]))
  if (length(not_na) == 0L) {
    priors_b_Sigma <- diag(2) # TMB can't handle 0-dim matrix
  } else {
    Sigma <- as.matrix(priors_b[, -1])
    priors_b_Sigma <- as.matrix(Sigma[not_na, not_na])
  }

  # expand sigma_V priors if needed
  if (nrow(priors_sigma_V) == 1L && ncol(X_rw_ik) > 1L) {
    priors_sigma_V <- t(replicate(ncol(X_rw_ik), priors_sigma_V, simplify = "matrix"))
  }
  if (nrow(priors_sigma_V) != ncol(X_rw_ik)) {
    cli_abort("sigma_V (time-varying SD) priors do not match the fitted model.")
  }

  if (!"A_st" %in% names(spde)) cli_abort("`mesh` was created with an old version of `make_mesh()`.")
  if (delta) y_i <- cbind(ifelse(y_i > 0, 1, 0), ifelse(y_i > 0, y_i, NA_real_))
  if (!delta) y_i <- matrix(y_i, ncol = 1L)

  # TODO: make this cleaner
  X_ij_list <- list()
  for (i in seq_len(n_m)) X_ij_list[[i]] <- X_ij[[i]]

  time_df <- make_time_lu(data[[time]], full_time_vec = union(data[[time]], extra_time))
  n_t <- nrow(time_df)

  random_walk <- if (!is.null(time_varying)) {
    switch(time_varying_type,
      rw = 1L,
      rw0 = 2L,
      ar1 = 0L
    )
  } else {
    0L
  }
  tmb_data <- list(
    y_i = y_i,
    n_t = n_t,
    z_i = z_i,
    offset_i = offset,
    proj_offset_i = 0,
    A_st = spde$A_st,
    sim_re = if ("sim_re" %in% names(experimental)) as.integer(experimental$sim_re) else rep(0L, 6),
    sim_obs = 1L,
    A_spatial_index = spde$sdm_spatial_id - 1L,
    year_i = time_df$year_i[match(data[[time]], time_df$time_from_data)],
    ar1_fields = ar1_fields,
    simulate_t = rep(1L, n_t),
    rw_fields = rw_fields,
    X_ij = X_ij_list,
    X_rw_ik = X_rw_ik,
    Zs = sm$Zs, # optional smoother basis function matrices
    Xs = sm$Xs, # optional smoother linear effect matrix
    proj_Zs = list(),
    proj_Xs = matrix(nrow = 0L, ncol = 0L),
    b_smooth_start = sm$b_smooth_start,
    proj_lon = 0,
    proj_lat = 0,
    proj_vector = 0,
    do_predict = 0L,
    calc_se = 0L,
    pop_pred = 0L,
    short_newdata = 0L,
    weights_i = if (!is.null(weights)) weights else rep(1, length(y_i)),
    area_i = rep(1, length(y_i)),
    normalize_in_r = 0L, # not used first time
    flag = 1L, # part of TMB::normalize()
    calc_index_totals = 0L,
    calc_cog = 0L,
    calc_eao = 0L,
    calc_weighted_avg = 0L,
    random_walk = random_walk,
    ar1_time = as.integer(!is.null(time_varying) && time_varying_type == "ar1"),
    priors_b_n = length(not_na),
    priors_b_index = not_na - 1L,
    priors_b_mean = priors_b[not_na, 1],
    priors_b_Sigma = priors_b_Sigma,
    priors_sigma_V = priors_sigma_V,
    priors = as.numeric(unlist(.priors)),
    share_range = as.integer(if (length(share_range) == 1L) rep(share_range, 2L) else share_range),
    include_spatial = as.integer(include_spatial), # changed later
    omit_spatial_intercept = as.integer(omit_spatial_intercept),
    proj_mesh = Matrix::Matrix(c(0, 0, 2:0), 3, 5), # dummy
    proj_X_ij = list(matrix(0, ncol = 1, nrow = 1)), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year = 0, # dummy
    proj_time_include = rep(1L, n_t),
    proj_spatial_index = 0, # dummy
    proj_z_i = matrix(0, nrow = 1, ncol = n_m), # dummy
    indexes_total = numeric(0), # dummy
    n_integration = 0L, # dummy
    spde_aniso = make_anisotropy_spde(spde, anisotropy),
    spde = get_spde_matrices(spde),
    barrier = as.integer(barrier),
    spde_barrier = make_barrier_spde(spde),
    barrier_scaling = if (barrier) spde$barrier_scaling else c(1, 1),
    anisotropy = as.integer(anisotropy),
    family = .valid_family[family$family],
    size = c(size),
    link = .valid_link[family$link],
    spatial_only = as.integer(spatial_only),
    spatial_covariate = as.integer(!is.null(spatial_varying)),
    calc_quadratic_range = as.integer(quadratic_roots),
    X_threshold = thresh[[1]]$X_threshold, # TODO: don't hardcode index thresh[[1]]
    proj_X_threshold = 0, # dummy
    threshold_func = thresh[[1]]$threshold_func, # TODO: don't hardcode index thresh[[1]]
    est_epsilon_model = as.integer(est_epsilon_model),
    epsilon_predictor = epsilon_covariate,
    est_epsilon_slope = as.integer(est_epsilon_slope),
    est_epsilon_re = as.integer(est_epsilon_re),
    has_smooths = as.integer(sm$has_smooths),
    upr = upr,
    lwr = 0L, # in case we want to reintroduce this
    poisson_link_delta = as.integer(isTRUE(family$type == "poisson_link_delta")),
    stan_flag = as.integer(bayesian),
    no_spatial = no_spatial,
    re_cov_df_map = as.matrix(re_cov_df_map), # dataframe used to map parameters to cov matrices,
    re_cov_df = as.matrix(re_cov_df),
    n_re_groups = c(n_re_groups),
    re_b_df = as.matrix(re_b_df[, -1]), # data frame containing indidivual level ids. Don't pass in first col (can be char)
    re_b_map = as.matrix(re_b_map),
    var_indx_matrix = var_indx_matrix,
    Zt_list = Zt_list, # list of RE matrices
    Zt_list_proj = list(),
    exclude_RE = 0L
  )
  b_thresh <- matrix(0, 2L, n_m)
  if (thresh[[1]]$threshold_func == 2L) b_thresh <- matrix(0, 3L, n_m) # logistic #TODO: change hard coding on index of thresh[[1]]

  tmb_params <- list(
    ln_H_input = matrix(0, nrow = 2L, ncol = n_m),
    b_j = rep(0, ncol(X_ij[[1]])), # TODO: verify ok
    b_j2 = if (delta) rep(0, ncol(X_ij[[2]])) else numeric(0), # TODO: verify ok
    bs = if (sm$has_smooths) matrix(0, nrow = ncol(sm$Xs), ncol = n_m) else array(0),
    ln_tau_O = rep(0, n_m),
    ln_tau_Z = matrix(0, n_z, n_m),
    ln_tau_E = rep(0, n_m),
    ln_kappa = matrix(0, 2L, n_m),
    # ln_kappa   = rep(log(sqrt(8) / median(stats::dist(spde$mesh$loc))), 2),
    thetaf = 0,
    ln_student_df = if (family$family[1] == "student") {
      if (estimate_student_df) log(2) else log(student_df_fixed - 1)
    } else {
      0  # default for non-student families
    },
    gengamma_Q = 0.5, # Not defined at exactly 0
    logit_p_extreme = 0,
    log_ratio_mix = -1, # ratio is 1 + exp(log_ratio_mix) so 0 would start fairly high
    ln_phi = rep(0, n_m),
    ln_tau_V = matrix(0, ncol(X_rw_ik), n_m),
    rho_time_unscaled = matrix(0, ncol(X_rw_ik), n_m),
    ar1_phi = rep(0, n_m),
    re_cov_pars = matrix(0, max_re_cov, n_m), # defined above, mapping off pars as needed
    re_b_pars = matrix(0, max_re_betas, n_m), # defined above, mapping off pars as needed
    b_rw_t = array(0, dim = c(tmb_data$n_t, ncol(X_rw_ik), n_m)),
    omega_s = matrix(0, if (!omit_spatial_intercept) n_s else 0L, n_m),
    zeta_s = array(0, dim = c(n_s, n_z, n_m)),
    epsilon_st = array(0, dim = c(n_s, tmb_data$n_t, n_m)),
    b_threshold = if (thresh[[1]]$threshold_func == 2L) matrix(0, 3L, n_m) else matrix(0, 2L, n_m),
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

  # Map off parameters not needed
  tmb_map <- map_all_params(tmb_params)
  tmb_map$b_j <- NULL
  if (delta) tmb_map$b_j2 <- NULL
  if (family$family[[1]] == "tweedie") tmb_map$thetaf <- NULL
  if (family$family[[1]] == "student") {
    if (estimate_student_df) {
      tmb_map$ln_student_df <- NULL  # estimate
    } else {
      tmb_map$ln_student_df <- factor(NA)  # fix at initial value
    }
  } else {
    tmb_map$ln_student_df <- factor(NA)  # not student family, fix at 0
  }
  if ("gengamma" %in% family$family) tmb_map$gengamma_Q <- factor(NA)
  tmb_map$ln_phi <- rep(1, n_m)
  if (family$family[[1]] %in% c("binomial", "poisson", "censored_poisson")) {
    tmb_map$ln_phi[1] <- factor(NA)
  }
  if (delta) {
    if (family$family[[2]] %in% c("binomial", "poisson", "censored_poisson")) {
      tmb_map$ln_phi[2] <- factor(NA)
    } else {
      tmb_map$ln_phi[2] <- 2
    }
  }
  tmb_map$ln_phi <- as.factor(tmb_map$ln_phi)
  if (!is.null(thresh[[1]]$threshold_parameter)) tmb_map$b_threshold <- NULL

  if (est_epsilon_re == 1L) {
    tmb_map <- unmap(tmb_map, c("ln_epsilon_re_sigma", "epsilon_re"))
  }
  if (est_epsilon_slope == 1L) {
    tmb_map <- unmap(tmb_map, "b_epsilon")
  }

  if (multiphase && is.null(previous_fit) && do_fit) {
    original_tmb_data <- tmb_data
    # much faster on first phase!?
    tmb_data$no_spatial <- 1L
    # tmb_data$include_spatial <- 0L
    tmb_data$include_spatial <- rep(0L, length(spatial)) # for 1st phase
    # tmb_data$spatial_only <- rep(1L, length(tmb_data$spatial_only))

    # Poisson on first phase increases stability:
    if (family$family[[1]] == "censored_poisson") tmb_data$family <- .valid_family["poisson"]

    tmb_obj1 <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params,
      profile = control$profile,
      map = tmb_map, DLL = "sdmTMB", silent = silent
    )
    lim <- set_limits(tmb_obj1, lower = lower, upper = upper, silent = TRUE)

    tmb_opt1 <- stats::nlminb(
      start = tmb_obj1$par, objective = tmb_obj1$fn,
      lower = lim$lower, upper = lim$upper,
      gradient = tmb_obj1$gr, control = .control
    )

    tmb_data <- original_tmb_data # restore
    # Set starting values based on phase 1:
    tmb_params <- tmb_obj1$env$parList()
    # tmb_data$no_spatial <- FALSE
    # often causes optimization problems if set from phase 1!?
    tmb_params$b_threshold <- if (thresh[[1]]$threshold_func == 2L) matrix(0, 3L, n_m) else matrix(0, 2L, n_m)
  }

  tmb_random <- c()
  if (any(spatial == "on") && !omit_spatial_intercept) {
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
    if (time_varying_type == "ar1") {
      tmb_map <- unmap(tmb_map, "rho_time_unscaled")
    }
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
  if (sum(n_re_groups) > 0) {
    tmb_random <- c(tmb_random, "re_b_pars")
    tmb_map <- unmap(tmb_map, c("re_cov_pars", "re_b_pars"))

    # map off pars not being estimated
    tmb_map$re_b_pars <- matrix(NA, max_re_betas, n_m)
    tmb_map$re_cov_pars <- matrix(NA, max_re_cov, n_m)
    for (m in seq_len(n_m)) {
      tmb_map$re_b_pars[seq(1, n_re_betas[m] + 1), m] <- 0 # these are estimated
      tmb_map$re_cov_pars[seq(1, n_re_cov[m] + 1), m] <- 0 # these are estimated
    }
    # turn non-NA values into sequence of integers
    non_nas <- which(!is.na(c(tmb_map$re_b_pars)))
    tmb_map$re_b_pars[non_nas] <- seq(1, length(non_nas))
    non_nas <- which(!is.na(c(tmb_map$re_cov_pars)))
    tmb_map$re_cov_pars[non_nas] <- seq(1, length(non_nas))
    tmb_map$re_b_pars <- as.factor(tmb_map$re_b_pars)
    tmb_map$re_cov_pars <- as.factor(tmb_map$re_cov_pars)
  }
  if (reml) tmb_random <- c(tmb_random, "b_j")
  if (reml && delta) tmb_random <- c(tmb_random, "b_j2")

  if (sm$has_smooths) {
    if (reml) tmb_random <- c(tmb_random, "bs")
    tmb_random <- c(tmb_random, "b_smooth") # smooth random effects
    tmb_map <- unmap(tmb_map, c("b_smooth", "ln_smooth_sigma", "bs"))
  }

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$tmb_obj$env$parList()
  }

  tmb_data$normalize_in_r <- as.integer(normalize)
  tmb_data$include_spatial <- as.integer(spatial == "on")

  if (!is.null(previous_fit)) tmb_map <- previous_fit$tmb_map

  # this is complex; pulled it out into own function:
  tmb_map$ln_kappa <- get_kappa_map(n_m = n_m, spatial = spatial, spatiotemporal = spatiotemporal, share_range = share_range)

  for (i in seq_along(start)) {
    cli_inform(c(
      i = paste0(
        "Initiating `", names(start)[i],
        "` at specified starting value(s) of:"
      ),
      paste0("  ", paste(round(start[[i]], 3), collapse = ", "))
    ))
    tmb_params[[names(start)[i]]] <- start[[i]]
  }

  if (!is.matrix(tmb_params[["ln_kappa"]]) && "ln_kappa" %in% names(start)) {
    msg <- c(
      "Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).",
      "It should be the same value in each row if `share_range = TRUE`."
    )
    cli_abort(msg)
  }
  if (nrow(tmb_params[["ln_kappa"]]) != 2L && "ln_kappa" %in% names(start)) {
    msg <- c(
      "Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).",
      "It should be the same value in each row if `share_range = TRUE`."
    )
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
      tmb_map$epsilon_st[, , i] <- NA
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
      tmb_map$omega_s[, i] <- NA
      tmb_map$ln_tau_O[i] <- NA
    }
    tmb_map$omega_s <- as.factor(tmb_map$omega_s)
    tmb_map$ln_tau_O <- as.factor(tmb_map$ln_tau_O)
  }

  if (anisotropy && delta && !"ln_H_input" %in% names(map)) {
    tmb_map$ln_H_input <- factor(c(1, 2, 1, 2)) # share anisotropy as in VAST
  }

  # Handle mixture models
  if (family$family[[1]] %in% c("gamma_mix", "lognormal_mix", "nbinom2_mix")) {
    fixed_p_extreme <- family$p_extreme

    if (!is.null(fixed_p_extreme)) {
      # Fix logit_p_extreme at the user-specified value
      tmb_map$logit_p_extreme <- factor(NA)
      tmb_params$logit_p_extreme <- log(fixed_p_extreme/(1-fixed_p_extreme))
      tmb_map$log_ratio_mix <- NULL
    } else {
      tmb_map$log_ratio_mix <- NULL
      tmb_map$logit_p_extreme <- NULL
    }
  }
  # delta mixture models
  if (delta) {
    if (family$family[[2]] %in% c("gamma_mix", "lognormal_mix", "nbinom2_mix")) {
      fixed_p_extreme <- family[[2]]$p_extreme
      if (is.null(fixed_p_extreme)) fixed_p_extreme <- family$p_extreme

      if (!is.null(fixed_p_extreme)) {
        # Fix logit_p_extreme at the user-specified value
        tmb_map$logit_p_extreme <- factor(NA)
        tmb_params$logit_p_extreme <- log(fixed_p_extreme/(1-fixed_p_extreme))
        tmb_map$log_ratio_mix <- NULL
      } else {
        tmb_map$log_ratio_mix <- NULL
        tmb_map$logit_p_extreme <- NULL
      }
    }
  }

  if (tmb_data$threshold_func > 0) tmb_map$b_threshold <- NULL
  if ("gengamma" %in% family$family) tmb_map$gengamma_Q <- NULL

  for (i in seq_along(map)) { # user supplied
    cli_inform(c(i = paste0("Fixing or mirroring `", names(map)[i], "`")))
    # ,
    #   "` at specified starting value(s) of:"),
    #   paste0("  ",
    #     paste(round(tmb_params[[names(map)[i]]], 3), collapse = ", "))))
  }
  tmb_map <- c(map, tmb_map) # add user-supplied mapping

  if (isTRUE(control$profile)) {
    prof <- c("b_j")
    if (delta) prof <- c(prof, "b_j2")
    control$profile <- prof
  } else if (is.character(control$profile)) {
    prof <- control$profile
    control$profile <- prof
  } else {
    control$profile <- NULL
  }

  out_structure <- structure(
    list(
      data = data,
      offset = offset,
      spde = spde,
      formula = original_formula,
      split_formula = split_formula,
      time_varying = time_varying,
      threshold_parameter = thresh[[1]]$threshold_parameter,
      threshold_function = thresh[[1]]$threshold_func,
      epsilon_predictor = epsilon_predictor,
      time = time,
      time_lu = time_df,
      family = family,
      smoothers = sm,
      response = y_i,
      tmb_data = tmb_data,
      tmb_params = tmb_params,
      tmb_map = tmb_map,
      tmb_random = tmb_random,
      spatial_varying = spatial_varying,
      spatial = spatial,
      spatiotemporal = spatiotemporal,
      spatial_varying_formula = spatial_varying_formula,
      reml = reml,
      priors = priors,
      nlminb_control = .control,
      control = control,
      contrasts = lapply(X_ij, attr, which = "contrasts"),
      terms = lapply(mf, attr, which = "terms"),
      extra_time = extra_time,
      fitted_time = sort(unique(data[[time]])),
      xlevels = lapply(seq_along(mf), function(i) stats::.getXlevels(mt[[i]], mf[[i]])),
      call = match.call(expand.dots = TRUE),
      version = utils::packageVersion("sdmTMB")
    ),
    class = "sdmTMB"
  )

  if ((!is.null(predict_args) || !is.null(index_args)) && isFALSE(do_index)) {
    cli_abort(c(
      "`predict_args` and/or `index_args` were set but `do_index` was FALSE.",
      "`do_index` must be TRUE for these arguments to take effect."
    )) # 276
  }
  if (do_index) {
    args <- list(object = out_structure, return_tmb_data = TRUE)
    args <- c(args, predict_args)
    if (!"newdata" %in% names(predict_args)) {
      cli_warn("`newdata` must be supplied if `do_index = TRUE`.")
    }
    tmb_data <- do.call(predict.sdmTMB, args)
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
      tmb_data$area_i <- rep(index_args[["area"]], length(tmb_data$proj_year)) # proj_year includes padded extra_time! otherwise, crash
    } else {
      if (length(index_args$area) != nrow(predict_args[["newdata"]])) {
        cli_abort("`area` length does not match `nrow(newdata)`.")
      }
      tmb_data$area_i <- index_args[["area"]]
    }
    tmb_data$calc_index_totals <- 1L
    tmb_params[["eps_index"]] <- numeric(0) # for bias correction
    out_structure$do_index <- TRUE
    do_index_time_missing_from_nd <-
      out_structure$do_index_time_missing_from_nd <-
      setdiff(data[[time]], predict_args$newdata[[time]])
  } else {
    out_structure$do_index <- FALSE
  }

  tmb_obj <- TMB::MakeADFun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    profile = control$profile,
    random = tmb_random, DLL = "sdmTMB", silent = silent
  )
  lim <- set_limits(tmb_obj,
    lower = lower, upper = upper,
    loc = spde$mesh$loc, silent = FALSE
  )

  out_structure$tmb_obj <- tmb_obj
  out_structure$tmb_data <- tmb_data
  out_structure$tmb_params <- tmb_params
  out_structure$lower <- lim$lower
  out_structure$upper <- lim$upper


  if (!do_fit) {
    return(out_structure)
  }

  if (normalize) tmb_obj <- TMB::normalize(tmb_obj, flag = "flag", value = 0)

  if (length(tmb_obj$par)) {
    tmb_opt <- stats::nlminb(
      start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
      lower = lim$lower, upper = lim$upper, control = .control
    )
  } else {
    tmb_opt <- list(par = tmb_obj$par, objective = tmb_obj$fn(tmb_obj$par))
  }

  if (isTRUE(suppress_nlminb_warnings)) {
    maybe_suppress_warnings <- suppressWarnings
  } else {
    maybe_suppress_warnings <- I
  }

  if (nlminb_loops > 1) {
    if (!silent) cli_inform("running extra nlminb optimization\n")
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- tmb_opt[c("iterations", "evaluations")]
      tmb_opt <- maybe_suppress_warnings(stats::nlminb(
        start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
        control = .control, lower = lim$lower, upper = lim$upper
      ))
      tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
      tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  if (!is.null(control$upper) || !is.null(control$lower)) {
    if (newton_loops > 0) {
      cli_inform("Upper or lower limits were set. `stats::optimHess()` will ignore these limits. Set `control = sdmTMBcontrol(newton_loops = 0)` to avoid the `stats::optimHess()` optimization if desired.")
    }
  }

  check_bounds(tmb_opt$par, lim$lower, lim$upper)

  # Check if spatial/spatiotemporal variances are collapsing to zero
  # Do this before Newton steps to avoid Hessian issues if variances have collapsed
  if (collapse_spatial_variance && length(tmb_obj$par) > 0) {
    report_vals <- tmb_obj$report()
    collapse_result <- check_and_collapse_random_fields(
      tmb_obj = tmb_obj,
      report_vals = report_vals,
      spatial = spatial,
      spatiotemporal = spatiotemporal,
      n_m = n_m,
      n_t = n_t,
      omit_spatial_intercept = omit_spatial_intercept,
      collapse_threshold = collapse_threshold,
      delta = delta,
      silent = silent
    )
    if (collapse_result$do_refit) {
      if (!silent) {
        cli_inform(c(
          "i" = "Refitting model with collapsed random field(s) using update()..."
        ))
      }

      # Now refit with collapsed fields disabled
      # The rest of this sdmTMB() function call was just executed with update()
      updated_fit <- update(
        out_structure,
        spatial = collapse_result$spatial_arg,
        spatiotemporal = collapse_result$spatiotemporal_arg,
        do_fit = TRUE,
        silent = silent
      )
      return(updated_fit)
    }
  }

  # We only end up here if no collapse was detected
  tmb_opt <- run_newton_loops(newton_loops = newton_loops, tmb_opt, tmb_obj, silent)

  if (!silent && getsd) cli_inform("running TMB sdreport\n")
  if (getsd) {
    sd_report <- TMB::sdreport(tmb_obj, getJointPrecision = get_joint_precision)
    conv <- get_convergence_diagnostics(sd_report)
  } else {
    sd_report <- NULL
    conv <- NULL
  }

  # save params that families need to grab from environments:
  if (any(family$family %in% c("truncated_nbinom1", "truncated_nbinom2"))) {
    phi <- exp(tmb_obj$par[["ln_phi"]])
    if (delta) {
      assign(".phi", phi, environment(out_structure[["family"]][[2]][["linkinv"]]))
    } else {
      assign(".phi", phi, environment(out_structure[["family"]][["linkinv"]]))
    }
  }

  out_structure$tmb_obj <- tmb_obj
  out <- c(out_structure, list(
    model = tmb_opt,
    sd_report = sd_report,
    parlist = tmb_obj$env$parList(par = tmb_obj$env$last.par.best),
    last.par.best = tmb_obj$env$last.par.best,
    gradients = conv$final_grads,
    bad_eig = conv$bad_eig,
    pos_def_hessian = sd_report$pdHess
  ))
  out <- `class<-`(out, "sdmTMB")
  # if (out$tmb_data$has_smooths) {
  #   edf <- suppressMessages(cAIC(out, what = "EDF"))
  #   out$edf <- edf
  # }
  out
}

check_bounds <- function(.par, lower, upper) {
  for (i in seq_along(.par)) {
    if (is.finite(lower[i]) || is.finite(upper[i])) {
      if (abs(.par[i] - lower[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its lower bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE
        )
      }
      if (abs(.par[i] - upper[i]) < 0.0001) {
        warning("Parameter ", names(.par)[i],
          " is very close or equal to its upper bound.\n",
          "Consider changing your model configuration or bounds.",
          call. = FALSE
        )
      }
    }
  }
}

check_and_collapse_random_fields <- function(
    tmb_obj,
    report_vals,
    spatial,
    spatiotemporal,
    n_m,
    n_t,
    omit_spatial_intercept,
    collapse_threshold,
    delta,
    silent) {

  do_refit <- FALSE
  spatial_updated <- spatial
  spatiotemporal_updated <- spatiotemporal

  # Check spatial field
  if (any(spatial == "on") && !omit_spatial_intercept) {
    est_sigma_O <- report_vals$sigma_O

    if (length(est_sigma_O) > 0 && any(est_sigma_O < collapse_threshold)) {
      which_sigma <- which(est_sigma_O < collapse_threshold)

      if (!silent) {
        cli_inform(c(
          "!" = "Spatial variance below threshold ({collapse_threshold}) detected",
          "i" = "Affected model(s): {paste(which_sigma, collapse = ', ')}",
          ">" = "Refitting with spatial field(s) disabled"
        ))
      }

      for (i in which_sigma) {
        spatial_updated[i] <- "off"
      }
      do_refit <- TRUE
    }
  }

  # Check spatiotemporal field (epsilon_st)
  if (!all(spatiotemporal == "off")) {
    est_sigma_E <- report_vals$sigma_E

    if (length(est_sigma_E) > 0) {
      est_sigma_E_by_model <- numeric(n_m)

      for (m in seq_len(n_m)) {
        idx_start <- (m - 1) * n_t + 1
        idx_end <- m * n_t
        if (idx_end <= length(est_sigma_E)) {
          est_sigma_E_by_model[m] <- min(est_sigma_E[idx_start:idx_end])
        }
      }

      if (any(est_sigma_E_by_model < collapse_threshold, na.rm = TRUE)) {
        which_sigma <- which(est_sigma_E_by_model < collapse_threshold)

        if (!silent) {
          cli_inform(c(
            "!" = "Spatiotemporal variance below threshold ({collapse_threshold}) detected",
            "i" = "Affected model(s): {paste(which_sigma, collapse = ', ')}",
            ">" = "Refitting with spatiotemporal field(s) disabled"
          ))
        }

        for (m in which_sigma) {
          spatiotemporal_updated[m] <- "off"
        }
        do_refit <- TRUE
      }
    }
  }

  # Prepare spatial/spatiotemporal arguments for update()
  if (delta) {
    spatial_arg <- as.list(spatial_updated)
    spatiotemporal_arg <- as.list(spatiotemporal_updated)
  } else {
    spatial_arg <- spatial_updated[1]
    spatiotemporal_arg <- spatiotemporal_updated[1]
  }

  list(
    do_refit = do_refit,
    spatial_arg = spatial_arg,
    spatiotemporal_arg = spatiotemporal_arg
  )
}

set_limits <- function(tmb_obj, lower, upper, loc = NULL, silent = TRUE) {
  .lower <- stats::setNames(rep(-Inf, length(tmb_obj$par)), names(tmb_obj$par))
  .upper <- stats::setNames(rep(Inf, length(tmb_obj$par)), names(tmb_obj$par))
  for (i_name in names(lower)) {
    if (i_name %in% names(.lower)) {
      .lower[names(.lower) %in% i_name] <- lower[[i_name]]
      if (!silent) {
        message(
          "Setting lower limit for ", i_name, " to ",
          lower[[i_name]], "."
        )
      }
    }
  }
  for (i_name in names(upper)) {
    if (i_name %in% names(.upper)) {
      .upper[names(.upper) %in% i_name] <- upper[[i_name]]
      if (!silent) {
        message(
          "Setting upper limit for ", i_name, " to ",
          upper[[i_name]], "."
        )
      }
    }
  }
  if ("ar1_phi" %in% names(tmb_obj$par) &&
    !"ar1_phi" %in% union(names(lower), names(upper))) {
    .lower["ar1_phi"] <- stats::qlogis((-0.999 + 1) / 2)
    .upper["ar1_phi"] <- stats::qlogis((0.999 + 1) / 2)
  }

  list(lower = .lower, upper = .upper)
}

check_spatiotemporal_arg <- function(x, time, .which = 1) {
  sp_len <- length(x)
  .spatiotemporal <- tolower(as.character(x[[.which]]))
  assert_that(.spatiotemporal %in% c("iid", "ar1", "rw", "off", "true", "false"),
    msg = "`spatiotemporal` argument value not valid"
  )
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
    m[[i]] <- factor(rep(NA, length(x[[i]])))
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

check_irregalar_time <- function(data, time, spatiotemporal, time_varying, extra_time) {
  if (any(spatiotemporal %in% c("ar1", "rw")) || !is.null(time_varying)) {
    if (!is.numeric(data[[time]])) {
      cli_abort("Time column should be integer or numeric if using AR(1) or random walk processes.")
    }
    ti <- sort(union(unique(data[[time]]), extra_time))
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
    mindiff <- 1L
    allx <- seq(min(ti), max(ti), by = mindiff)
    setdiff(allx, ti)
  }
}

unmap <- function(x, v) {
  for (i in v) x[[i]] <- NULL
  x
}

get_spde_matrices <- function(x) {
  x <- x$spde[c("c0", "g1", "g2")]
  names(x) <- c("M0", "M1", "M2") # legacy INLA names needed!
  x
}
