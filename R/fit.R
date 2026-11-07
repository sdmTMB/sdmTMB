#' @useDynLib sdmTMB, .registration = TRUE
NULL

#' Fit a spatial or spatiotemporal GLMM with TMB
#'
#' Fit a spatial or spatiotemporal generalized linear mixed effects model (GLMM)
#' with \pkg{TMB} (Template Model Builder), by default through the \pkg{RTMB}
#' package. Spatial and spatiotemporal Gaussian random fields are approximated
#' with the SPDE (stochastic partial differential equation) approach, which
#' represents them as Gaussian Markov random fields on a mesh. This allows for
#' efficient modelling of data that are correlated in space and/or time. Areal
#' spatial and spatiotemporal models (simultaneous or conditional
#' autoregressive; SAR or CAR) are also available. See the
#' [model description vignette](https://sdmTMB.github.io/sdmTMB/articles/model-description.html)
#' for details.
#'
#' @param formula A model formula. Random intercepts and slopes use \pkg{lme4}
#'   syntax, e.g., `(1 | g)`, `(0 + depth | g)`, or `(1 + depth | g)`, where
#'   `g` is a character or factor column. As in \pkg{lme4}, intercepts and
#'   slopes within a term are correlated. Penalized smooths use \pkg{mgcv}
#'   syntax, e.g., `s(depth)`, and threshold terms use `breakpt()` or
#'   `logistic()`. For delta models, optionally a list of two formulas. See
#'   Details.
#' @param data A data frame. Rows with missing values in the response or any
#'   variable the model uses (including `weights` and `offset`) are omitted
#'   before fitting, as with `na.action = na.omit` in [stats::glm()]. The
#'   rows used are stored in the returned object as `data`.
#' @param mesh An object from [make_mesh()] for `spatial_model = "spde"` or
#'   from [make_areal_domain()] for `"sar"` or `"car"`. Can be omitted for a
#'   model without spatial or spatiotemporal random fields, in which case
#'   `spatial` and `spatiotemporal` default to `'off'`.
#' @param time An optional time column name (as character). Can be left as
#'   `NULL` for a model with only spatial random fields; however, if the data
#'   are actually spatiotemporal and you wish to calculate derived quantities
#'   downstream (e.g., [get_index()] or [get_cog()]), then supply the time argument.
#' @param family A family object specifying the response distribution and
#'   link: [gaussian()], [Gamma()], [binomial()], [poisson()], or one of the
#'   additional families listed in [Families] (e.g., [tweedie()],
#'   [nbinom2()], [delta_gamma()]). See 'Binomial families', 'Censored
#'   families', and 'Delta/hurdle models' in Details. Experimental multi-family models take a named list of
#'   families; `distribution_column` then assigns each row to one of them. See
#'   the [multi-family
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/multi-family.html) for
#'   supported combinations.
#' @param spatial Estimate spatial random fields? Options are `'on'` / `'off'`
#'   or equivalently `TRUE` / `FALSE`. Defaults to `'on'` if `mesh` is supplied
#'   and `'off'` otherwise. Optionally, a list for delta models,
#'   e.g. `list('on', 'off')`.
#' @param spatiotemporal Estimate the spatiotemporal random fields as `'iid'`
#'   (independent and identically distributed; default), stationary `'ar1'`
#'   (first-order autoregressive), a random walk (`'rw'`), or fixed at 0
#'   `'off'`. Will be set to `'off'` if `time = NULL` or, by default, if `mesh`
#'   is not supplied. If a delta model, can be
#'   a list. E.g., `list('off', 'ar1')`. Guidance: Use `'iid'` if spatiotemporal
#'   correlation is negligible or already accounted for in fixed effects; `'ar1'`
#'   if correlation between consecutive time steps decays gradually; `'rw'` if
#'   changes between time steps are cumulative (each step builds on the last). If
#'   the AR1 correlation coefficient (rho) is estimated close to 1 (say > 0.99),
#'   consider switching to `'rw'`. See the [model description
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/model-description.html) for
#'   mathematical details. Capitalization is ignored. `TRUE` gets converted to
#'   `'iid'` and `FALSE` gets converted to `'off'`.
#' @param offset A numeric vector representing the model offset *or* a character
#'   value representing the column name of the offset. In delta/hurdle models,
#'   this applies only to the positive component except for Poisson-link delta
#'   models, where it also enters the occurrence-probability calculation.
#'   Usually a log transformed variable.
#' @param weights An optional numeric vector (not a column name) of weights on
#'   each observation's contribution to the likelihood. As in \pkg{glmmTMB},
#'   weights need not sum to one and are not rescaled. For binomial-type
#'   families with a proportion response, `weights` gives the number of trials
#'   instead; see 'Binomial families' in Details.
#' @param spatial_model Spatial process model. `"spde"` uses the default
#'   continuous-space SPDE approximation. `"sar"` and `"car"` use areal spatial
#'   autoregressive models with an areal domain supplied to `mesh`.
#'   Capitalization is ignored.
#' @param share_range Logical: estimate a single range parameter shared by the
#'   spatial and spatiotemporal fields (`TRUE`, default) or separate range
#'   parameters (`FALSE`)? For delta models, can be a list, e.g.,
#'   `list(TRUE, FALSE)`. Applies to SPDE models only.
#' @param anisotropy Logical: allow for anisotropy (spatial correlation that is
#'   directionally dependent)? See [plot_anisotropy()]. Applies to SPDE models
#'   only and is shared by both components of a delta model.
#' @param extra_time Optional extra time slices (e.g., years) to include for
#'   interpolation or forecasting with the predict function. See the Details
#'   section below.
#' @param spatial_varying An optional one-sided formula of coefficients that
#'   should vary in space as random fields. Allows the effect of a covariate to
#'   differ spatially. You likely want to include the same variable as a fixed
#'   effect in `formula` to estimate the average effect—the spatial field then
#'   represents deviations from that average. For example,
#'   use `formula = y ~ depth` and `spatial_varying = ~ 0 + depth` to model an
#'   average depth effect plus spatially varying deviations. If a (scaled) time
#'   column is used, this creates a local-time-trend model. See
#'   \doi{10.1111/ecog.05176} and the [spatial trends
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/spatial-trend-models.html).
#'   Predictors should usually be centered to have mean zero and standard deviation
#'   approximately 1. **The spatial intercept is controlled by the `spatial`
#'   argument**; set `spatial = 'on'` or `'off'` to include or exclude it. For a
#'   factor, `~ 0 + f` gives every level its own field, which would duplicate
#'   the spatial intercept field, so set `spatial = 'off'` to match. Structure
#'   is shared in delta models.
#' @param time_varying An optional one-sided formula of coefficients that vary
#'   through time following the process set by `time_varying_type`. Whether
#'   the same covariates should also appear in `formula` depends on that type.
#'   Shared by both components of a delta model.
#' @param time_varying_type The process for `time_varying` coefficients:
#'   * `'rw'` (default): a random walk with the first value estimated freely.
#'     Do not also include these covariates in `formula`, or the model is not
#'     identifiable; e.g., with `time_varying = ~ 1`, use `formula = y ~ 0 +
#'     ...`.
#'   * `'rw0'`: a random walk whose first value has a mean-zero normal prior.
#'   * `'ar1'`: a stationary first-order autoregressive process with mean
#'     zero.
#'
#'   For `'rw0'` and `'ar1'`, include the same covariates in `formula`; the
#'   time-varying process then describes deviations from that average effect.
#'   Shared by both components of a delta model.
#' @param dispformula A one-sided formula for the dispersion parameter (e.g.,
#'   `phi`). The default, `~ 1`, estimates a single value. Ignored for
#'   families without a dispersion parameter (e.g., binomial or Poisson). In
#'   delta models, applies to the positive component. Not yet available for
#'   multi-family models or truncated negative binomial families.
#' @param nonlocal_formula An optional one-sided formula of non-local
#'   covariate effects: [diffusion()] for an effect of conditions in the
#'   surrounding area, and [time_lag()] for an effect of conditions in
#'   previous time steps, e.g., `~ diffusion(x) + time_lag(x)`. Using both on
#'   the same covariate gives one combined effect that spreads over both space
#'   and time; different covariates give separate effects. See [time_lag()]
#'   for how the lag starts (`start`). If `time` is `NULL`, a spatial-only
#'   covariate is held constant across time slices. See the [non-local
#'   covariates
#'   vignette](https://sdmTMB.github.io/sdmTMB/articles/nonlocal-covariates.html),
#'   which also explains the reported diffusion scales (MSDK and RMSDK).
#' @param nonlocal_data An optional data frame supplying the
#'   `nonlocal_formula` covariate(s) at a different resolution and/or
#'   coverage than `data` (e.g., a finer grid, or one spanning `extra_time`
#'   slices). Must contain the mesh `xy_cols`, the diffusion covariate
#'   columns, and the `time` column if `nonlocal_formula` is time-indexed
#'   (`time_lag()` terms, or `diffusion()` terms with `time` specified). In
#'   that case, it must cover every fitted (+ `extra_time`) time slice.
#'   Defaults to `NULL`, in which case `data` is used.
#' @param knots Optional named list containing knot values to be used for basis
#'   construction of smoothing terms. See [mgcv::gam()] and [mgcv::gamm()].
#'   E.g., `s(x, bs = 'cc', k = 4), knots = list(x = c(1, 2, 3, 4))`
#' @param censored_upper Upper bounds for the censored families (see
#'   [Families]): a numeric vector or the name of a column in `data`. Each
#'   response is a count between its observed value and this bound. A bound
#'   equal to the response is uncensored, a larger bound is interval-censored,
#'   and `Inf` is right-censored (capped at the number of trials for
#'   [censored_binomial()] and [censored_betabinomial()]). Non-integer bounds
#'   are rounded down. For left censoring (e.g., fewer than 5), use a response
#'   of 0 and a bound of 4.
#' @param distribution_column For experimental multi-family models, the name of
#'   the column in `data` mapping each row to a family in the named `family`
#'   list. See the multi-family vignette for the supported family and method
#'   combinations.
#' @param silent Silent or include optimization details? Helpful to set to
#'   `FALSE` for models that take a while to fit.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation
#'   rather than maximum likelihood? REML accounts for uncertainty in estimating
#'   fixed effects and can reduce bias in variance parameter estimates, but
#'   prevents likelihood-based model comparison (e.g., AIC) between models with
#'   different fixed effects. Use `TRUE` if your focus is on random effect
#'   variance parameters; use `FALSE` (default) if comparing models with different
#'   fixed effects or performing index standardization.
#' @param priors Optional penalties/priors via [sdmTMBpriors()]. Must currently
#'   be shared across delta models.
#' @param bayesian Logical indicating if the model will be passed to
#'   \pkg{tmbstan}. If `TRUE`, Jacobian adjustments are applied to account for
#'   parameter transformations when priors are applied.
#' @param control Optimization control options via [sdmTMBcontrol()].
#' @param previous_fit A previously fitted sdmTMB model to initialize the
#'   optimization with. Can greatly speed up fitting. Note that the model must
#'   be set up *exactly* the same way. However, the data and `weights` arguments
#'   can change, which can be useful for cross-validation.
#' @param do_fit Fit the model (`TRUE`) or return the processed data without
#'   fitting (`FALSE`)?
#' @param experimental A named list for esoteric or in-development options. Here
#'   be dragons.
#' @param do_index,predict_args,index_args `r lifecycle::badge("defunct")`
#'   Index calculations while fitting are no longer supported. Instead, fit the
#'   model and then pass it with `newdata` directly to [get_index()],
#'   [get_cog()], [get_eao()], or [get_weighted_average()].
#' @importFrom methods as is
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom mgcv s t2
#' @importFrom rlang %||%
#' @importFrom stats gaussian model.frame model.matrix as.formula
#' @importFrom stats model.response terms model.offset
#' @importFrom lifecycle deprecated is_present deprecate_warn
#'
#' @return
#' An object (list) of class `sdmTMB`. Useful elements include:
#'
#' * `model`: output from [stats::nlminb()]
#' * `sd_report`: output from [TMB::sdreport()] or [RTMB::sdreport()]
#' * `gradients`: gradients of the marginal log likelihood with respect to
#'   each fixed effect
#' * `data`: the data used in fitting (rows with missing values removed)
#' * `spde`: the object supplied to `mesh`
#' * `family`: the family object, including the inverse link function
#'   `family$linkinv()`
#' * `backend`: `"rtmb"` or `"tmb"`; see [sdmTMBcontrol()]
#' * `tmb_obj`: the objective function object from [TMB::MakeADFun()] or
#'   [RTMB::MakeADFun()]
#' * `tmb_data`, `tmb_params`, `tmb_map`: the data, parameter, and map lists
#'   passed to `MakeADFun()`
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
#' family can be specified in one of four ways: (1) the response may be a factor
#' (success is interpreted as any level other than the first level), (2) the
#' response may be binary (0/1), (3) the response can be a matrix of form
#' `cbind(success, failure)`, and (4) the response may be observed
#' proportions, and the `weights` argument is used to specify the binomial size
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
#' supplied).
#'
#' **Threshold models**
#'
#' A linear break-point relationship for a covariate can be included via
#' `+ breakpt(variable)` in the formula, where `variable` is a single covariate
#' corresponding to a column in `data`. In this case, the relationship is linear
#' up to a point and then constant (hockey-stick shaped). The likelihood is not
#' smooth in the cutpoint, so a fit can stop at the wrong cutpoint while
#' reporting convergence, and its standard error for the cutpoint is not
#' reliable. [profile_breakpt()] estimates the cutpoint and a confidence interval
#' by profile likelihood instead.
#'
#' Similarly, a logistic-function threshold model can be included via
#' `+ logistic(variable)`. This option models the relationship as a logistic
#' function of the 50% and 95% values. This is similar to length- or size-based
#' selectivity in fisheries, and is parameterized by the points at which f(x) =
#' 0.5 or 0.95. See the
#' [threshold vignette](https://sdmTMB.github.io/sdmTMB/articles/threshold-models.html).
#'
#' Note that only a single threshold covariate can be included. For delta
#' families, the threshold applies to both components, and threshold terms are
#' not available if `formula` is a list.
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
#' time slices; `time_varying` and spatiotemporal random-walk or AR(1) fields
#' include process error.
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
#' or at the same time by using an appropriate delta family (see the list under
#' `family`). Delta families with `type = "poisson-link"` use a Poisson-link
#' parameterization instead of a classic hurdle model; see the [Poisson-link
#' vignette](https://sdmTMB.github.io/sdmTMB/articles/poisson-link.html).
#' If fit with a delta family, by default the formula, spatial, and spatiotemporal
#' components are shared. Some elements can be specified independently for the two models
#' using a list format. These include `formula`, `spatial`, `spatiotemporal`,
#' and `share_range`. The first element of the list is for the binomial component
#' and the second element is for the positive component (e.g., Gamma).
#' Other elements must be shared for now (e.g., spatially varying coefficients,
#' time-varying coefficients), and `dispformula` applies to the positive
#' component only. Furthermore, there are currently limitations if
#' specifying two formulas as a list: smoothers and random effect terms must
#' be identical between the two formulas, and threshold effects must be
#' specified through a single formula that is shared across the two models.
#'
#' The main advantage of specifying such models using a delta family (compared
#' to fitting two separate models) is (1) coding simplicity and (2) calculation
#' of uncertainty on derived quantities such as an index of abundance with
#' [get_index()] using the generalized delta method within TMB. Also, selected
#' parameters can be shared across the models.
#'
#' See the [delta-model vignette](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html).
#'
#' **Censored families**
#'
#' Censored families treat some observations as known only to lie within a
#' range, e.g., catch counts on longlines where competition for hooks hides
#' the true number caught. Supply the bounds with `censored_upper`. All but
#' [censored_poisson()] need the RTMB backend (the default). See [Families]
#' and the [hook competition
#' vignette](https://sdmTMB.github.io/sdmTMB/articles/hook-competition.html).
#'
#' **Areal models**
#'
#' For data on areal units (e.g., grid cells or management areas) rather than
#' point locations, build a domain with [make_areal_domain()], pass it to
#' `mesh`, and set `spatial_model = "sar"` or `"car"`. `share_range` and
#' `anisotropy` do not apply. See the [areal model
#' vignette](https://sdmTMB.github.io/sdmTMB/articles/areal-sar-car-spde.html).
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
#' Anderson, S.C., E.J. Ward, P.A. English, L.A.K. Barnett, J.T. Thorson. 2025.
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
#' A.M. Edwards, B.M. Connors, S.C. Anderson. 2022. Contrasting climate velocity
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
#' Lindmark, M., S.C. Anderson, M. Gogina, M. Casini. 2023. Evaluating drivers
#' of spatiotemporal variability in individual condition of a bottom-associated
#' marine fish, Atlantic cod (*Gadus morhua*). ICES Journal of Marine Science.
#' 80(5): 1539--1550.
#'
#' *Non-local covariates:*
#'
#' Lindmark, M., Anderson, S.C., and Thorson, J.T. 2026. Estimating scale-dependent
#' covariate responses using two-dimensional diffusion derived from the stochastic
#' partial differential equation method. Methods in Ecology and Evolution 17:
#' 207–218. \doi{10.1111/2041-210X.70177}.
#'
#' Thorson, J.T., Anderson, S.C., and Lindmark, M. 2026.
#' Temperature carryover effect revealed for marine fishes using spatio-temporal
#' distributed lag models. EcoEvoRxiv. \doi{10.32942/X2W95P}.
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
#' # - or build any mesh in 'fmesher' and pass it to the 'mesh' argument in make_mesh()
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
#' visreg::visreg(fit, xvar = "depth", scale = "response", rug = FALSE)
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
#' # Which matches glm():
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
#' # Correlated random intercepts and slopes by year:
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
    offset = NULL,
    weights = NULL,
    spatial_model = c("spde", "sar", "car"),
    share_range = TRUE,
    anisotropy = FALSE,
    extra_time = NULL,
    spatial_varying = NULL,
    time_varying = NULL,
    time_varying_type = c("rw", "rw0", "ar1"),
    dispformula = ~ 1,
    nonlocal_formula = NULL,
    nonlocal_data = NULL,
    knots = NULL,
    censored_upper = NULL,
    distribution_column = NULL,
    silent = TRUE,
    reml = FALSE,
    priors = sdmTMBpriors(),
    bayesian = FALSE,
    control = sdmTMBcontrol(),
    previous_fit = NULL,
    do_fit = TRUE,
    experimental = NULL,
    do_index = FALSE,
    predict_args = NULL,
    index_args = NULL) {
  if (!isFALSE(do_index) || !is.null(predict_args) || !is.null(index_args)) {
    lifecycle::deprecate_stop(
      "1.2.0",
      "sdmTMB(do_index = )",
      details = paste(
        "Fit the model and then pass it with `newdata` directly to `get_index()`,",
        "`get_cog()`, `get_eao()`, or `get_weighted_average()`."
      )
    )
  }

  assert_that(class(formula) %in% c("formula", "list"))
  assert_that(inherits(data, "data.frame"))
  if (!is.null(time)) assert_that(is.character(time))
  nonlocal_data_arg <- nonlocal_data
  mesh_missing <- missing(mesh) || is.null(mesh)
  spatial_model <- match.arg(tolower(spatial_model[1L]), c("spde", "sar", "car"))
  if (mesh_missing && spatial_model %in% c("sar", "car")) {
    cli_abort("`spatial_model = \"{spatial_model}\"` requires an areal domain supplied to `mesh`.")
  }
  if (!inherits(dispformula, "formula") || length(dispformula) != 2L) {
    cli_abort("`dispformula` must be a one-sided formula such as `~ 1`.")
  }
  censored_upper <- .censored_upper_arg(censored_upper, data)
  # The unfiltered row inputs, for the refit after collapsing spatial fields
  refit_inputs <- list(data = data, mesh = if (!mesh_missing) mesh, weights = weights,
    offset = offset, censored_upper = censored_upper, experimental = experimental,
    control = control)
  # Omit rows with missing values in any variable the model uses, as with
  # `na.action = na.omit` in glm(). Doing this before anything is built from
  # `data` keeps the stored data, mesh rows, and post-fit methods aligned.
  complete <- .complete_model_rows(
    data,
    formulas = c(
      .formula_list(formula), .formula_list(spatial_varying),
      .formula_list(time_varying), list(dispformula),
      if (is.null(nonlocal_data)) .formula_list(nonlocal_formula)
    ),
    row_vectors = list(weights, if (is.character(offset)) data[[offset]] else offset)
  )
  if (!any(complete)) cli_abort("No rows remain after omitting rows with missing values.")
  if (!all(complete)) {
    rows <- which(complete)
    n <- nrow(data)
    if (!silent) cli_inform("Omitting {n - length(rows)} row{?s} with missing values.")
    data <- data[rows, , drop = FALSE]
    weights <- .subset_rows(weights, rows, n)
    offset <- .subset_rows(offset, rows, n)
    censored_upper <- .subset_rows(censored_upper, rows, n)
    if (!is.null(experimental$.cv_fold_weights)) {
      experimental$.cv_fold_weights <- .subset_rows(experimental$.cv_fold_weights, rows, n)
    }
    if (!mesh_missing) mesh <- .subset_mesh_rows(mesh, rows, n)
  }
  data <- droplevels(data) # if data was subset, strips absent factors
  disp_terms <- stats::terms(dispformula)
  has_dispformula <- length(attr(disp_terms, "term.labels")) > 0L || attr(disp_terms, "intercept") != 1L

  family_spec <- .compile_family_spec(
    family = family,
    data = data,
    distribution_column = distribution_column
  )
  is_multi_family <- .family_spec_is_multi_family(family_spec)
  family <- family_spec$family
  has_two_components <- .family_spec_has_two_components(family_spec)
  n_m <- family_spec$n_m
  check_component_list(formula, n_m, "formula")
  check_component_list(if (!missing(spatial)) spatial, n_m, "spatial")
  check_component_list(if (!missing(spatiotemporal)) spatiotemporal, n_m, "spatiotemporal")
  check_component_list(share_range, n_m, "share_range", vector_ok = TRUE)
  if (is_multi_family && has_dispformula) {
    cli_abort("`dispformula` is not supported for multi-family models yet. Use `dispformula = ~ 1`.")
  }
  if (!is_multi_family) {
    disp_family <- if (has_two_components) family$family[[2]] else family$family[[1]]
    if (has_dispformula && !.family_metadata(disp_family)$uses_phi) {
      cli_warn("`dispformula` ignored because this family has no dispersion parameter.")
      dispformula <- ~ 1
      has_dispformula <- FALSE
    }
    if (has_dispformula && disp_family %in% c("truncated_nbinom1", "truncated_nbinom2")) {
      cli_abort("`dispformula` is not supported for truncated negative-binomial families yet.")
    }
  }

  if (inherits(formula, "formula")) {
    formula <- normalize_bar_only_formula(formula)
  } else {
    formula <- lapply(formula, normalize_bar_only_formula)
  }

  if (!missing(spatiotemporal)) {
    if (has_two_components && !is.list(spatiotemporal)) {
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
    if (is.null(time) || mesh_missing) {
      spatiotemporal <- rep("off", n_m)
    } else {
      spatiotemporal <- rep("iid", n_m)
    }
  }

  if (missing(spatial) && mesh_missing) spatial <- "off"
  if (is.list(spatial)) {
    spatial <- vapply(spatial, parse_spatial_arg, FUN.VALUE = character(1L))
  } else {
    spatial <- rep(parse_spatial_arg(spatial), n_m)
  }
  include_spatial <- "on" %in% spatial
  spatial_user <- spatial  # preserve user's specification before internal modification

  if (!include_spatial && !is.null(spatial_varying)) {
    # spatial = "off" + SVC: keep SPDE infrastructure alive so SVC fields can
    # estimate range and marginal SD, but omega_s is omitted from the linear
    # predictor downstream (see the TMB template).
    omit_spatial_intercept <- TRUE
    include_spatial <- TRUE
    spatial <- rep("on", length(spatial)) # internal only; spatial_user retains "off"
  } else {
    omit_spatial_intercept <- FALSE
  }

  no_spatial <- !include_spatial && all(spatiotemporal == "off")
  if (no_spatial && mesh_missing) mesh <- sdmTMB::pcod_mesh_2011 # internal data; fake!
  if (!no_spatial && mesh_missing) {
    cli_abort(c(
      "`mesh` must be supplied for spatial, spatiotemporal, or spatially varying fields.",
      "i" = "Create one with `make_mesh()` or set `spatial = \"off\"` and `spatiotemporal = \"off\"`."
    ))
  }

  share_range <- unlist(share_range)
  share_range_user <- share_range
  if (length(share_range) == 1L) share_range <- rep(share_range, n_m)
  share_range[spatiotemporal == "off"] <- TRUE
  share_range[spatial == "off"] <- TRUE

  spde <- mesh
  if (any(c("epsilon_model", "epsilon_predictor") %in% names(experimental))) {
    cli_abort("`experimental$epsilon_model` and `experimental$epsilon_predictor` have been removed.")
  }

  backend <- control$backend
  if (backend == "rtmb" && isTRUE(control$normalize)) {
    cli_abort("`normalize = TRUE` is unavailable with the RTMB backend; its densities are fully normalized.")
  }
  if (isTRUE(control$profile)) {
    control$profile <- if (has_two_components) c("b_j", "b_j2") else "b_j"
  } else if (!is.character(control$profile)) {
    control$profile <- NULL
  }
  # options that aren't sdmTMBcontrol() arguments go to nlminb(); "newton_steps"
  # and "mgcv" are in controls saved by old versions
  sdmTMB_only <- c(setdiff(names(formals(sdmTMBcontrol)), c("eval.max", "iter.max", "...")),
    "newton_steps", "mgcv")
  .control <- control[setdiff(names(control), sdmTMB_only)]

  assert_that(
    is.logical(reml), is.logical(anisotropy), is.logical(share_range),
    is.logical(silent), is.logical(control$multiphase), is.logical(control$normalize)
  )

  if (!is.null(spatial_varying)) assert_that(class(spatial_varying) %in% c("formula", "list"))
  if (!is.null(time_varying)) assert_that(class(time_varying) %in% c("formula", "list"))
  if (!is.null(previous_fit)) assert_that(identical(class(previous_fit), "sdmTMB"))
  assert_that(is.list(priors))
  priors_custom <- custom_prior_spec(priors, backend)
  assert_that(is.list(.control))
  time_varying_type <- match.arg(time_varying_type)

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

  user_time_supplied <- !is.null(time)
  nonlocal_formula_parsed <- .parse_nonlocal_formula(nonlocal_formula)
  nonlocal_formula_parsed <- .validate_nonlocal_terms(
    nonlocal_formula_parsed,
    data = data,
    time = time,
    multi_family = is_multi_family
  )
  if (!is.null(nonlocal_formula_parsed) && mesh_missing) {
    cli_abort("`mesh` must be supplied when using `nonlocal_formula`.")
  }
  if (!is.null(nonlocal_data_arg) && is.null(nonlocal_formula_parsed)) {
    cli_abort("`nonlocal_data` was supplied but `nonlocal_formula` was not.")
  }
  nonlocal_time_indexed <- .nonlocal_time_indexed(
    nonlocal_formula_parsed,
    time_supplied = user_time_supplied
  )
  if (!is.null(nonlocal_formula_parsed)) {
    nonlocal_formula_parsed$time_indexed <- nonlocal_time_indexed
  }

  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }

  time_df <- make_time_lu(data[[time]], full_time_vec = union(data[[time]], extra_time))
  n_t <- nrow(time_df)
  year_i_data <- time_df$year_i[match(data[[time]], time_df$time_from_data)]

  domain <- prepare_spatial_domain(
    mesh = spde,
    data = data,
    mesh_missing = mesh_missing,
    share_range = share_range,
    anisotropy = anisotropy,
    nonlocal_formula = nonlocal_formula_parsed,
    priors = priors,
    normalize = control$normalize,
    share_range_user = share_range_user,
    spatial_model = spatial_model,
    sar_weight_style = control$sar_weight_style
  )
  is_areal <- identical(domain$type, "areal")

  nonlocal_grid_supplied <- !is.null(nonlocal_data_arg)
  if (nonlocal_grid_supplied) {
    nonlocal_grid_inputs <- .prepare_nonlocal_grid_inputs(
      grid = nonlocal_data_arg,
      nonlocal_formula = nonlocal_formula_parsed,
      mesh = spde$mesh,
      xy_cols = spde$xy_cols,
      time = time,
      time_df = time_df,
      full_time_vec = time_df$time_from_data,
      time_indexed = nonlocal_time_indexed
    )
  }

  # Covariate diffusion also maps observations to the mesh by row
  uses_mesh_rows <- !no_spatial || !is.null(nonlocal_formula_parsed)
  if (uses_mesh_rows && !is_areal) {
    if (!identical(nrow(spde$loc_xy), nrow(data))) {
      msg <- c(
        "Number of x-y coordinates in `mesh` does not match `nrow(data)`.",
        "Is it possible you passed a different data frame to `make_mesh()` and `sdmTMB()`?" # discussion 137
      )
      cli::cli_abort(msg)
    }
  }

  rtmb_only <- intersect(family_spec$components$family_name, c("censored_nbinom1",
    "censored_nbinom2", "censored_binomial", "censored_betabinomial", "compois"))
  if (length(rtmb_only) && backend != "rtmb") {
    cli_abort("`{rtmb_only[[1]]}()` needs `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  uses_censored <- any(family_spec$components$family_name %in%
    c("censored_poisson", "censored_nbinom1", "censored_nbinom2",
      "censored_binomial", "censored_betabinomial"))
  if (uses_censored) {
    if (is.null(censored_upper)) cli_abort("`censored_upper` must be supplied to use a censored family.")
  }
  if (is.null(censored_upper)) censored_upper <- Inf

  original_formula <- if (inherits(formula, "formula")) replicate(n_m, list(formula)) else formula
  thresh <- parse_thresholds(original_formula, data)
  formula <- thresh$formulas

  if (is.character(offset)) {
    offset <- data[[offset]]
  }

  check_irregalar_time(
    data,
    time,
    spatiotemporal,
    time_varying,
    extra_time = extra_time,
    nonlocal_temporal = !is.null(nonlocal_formula_parsed) &&
      isTRUE(nonlocal_formula_parsed$time_indexed) &&
      !nonlocal_grid_supplied
  )

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
    svc_design <- .make_design(spatial_varying, data, "`spatial_varying`")
    z_i <- svc_design$X
    .int <- grep("(Intercept)", colnames(z_i))
    has_intercept <- length(.int) > 0L
    svc_omega_is_intercept <- has_intercept && !omit_spatial_intercept
    if (svc_omega_is_intercept) {
      # spatial = "on" + intercept in `spatial_varying`: alias the ordinary
      # spatial field (omega_s) as the SVC intercept/reference-level field and
      # drop the duplicate intercept column from the SVC design matrix.
      if (!silent) {
        cli_inform(c(
          "i" = "Using the ordinary spatial field (`omega_s`) as the SVC intercept/reference-level field.",
          " " = "The `(Intercept)` column was removed from the `spatial_varying` design matrix to avoid a duplicate intercept field."
        ))
      }
      z_i <- z_i[, -.int, drop = FALSE]
    }
    if (has_intercept && omit_spatial_intercept &&
        length(attr(z_i, "contrasts"))) {
      cli_inform(c(
        "The behaviour of `spatial = \"off\"` with `~ 1 + factor` in `spatial_varying` has changed.",
        "It now fits a genuine SVC intercept field plus `K - 1` deviation fields (previously the intercept was dropped, leaving only the deviations with no reference field).",
        "To use the ordinary spatial field as the reference-level field instead (same resulting model), set `spatial = \"on\"`."
      ))
    }
    spatial_varying <- colnames(z_i)
    svc_design$spec$columns <- colnames(z_i)
  } else {
    z_i <- matrix(0, nrow(data), 0L)
    svc_design <- NULL
    svc_omega_is_intercept <- FALSE
  }
  n_z <- ncol(z_i)

  if (any(grepl("offset\\(", formula))) {
    cli_abort("Detected `offset()` in formula. Offsets in sdmTMB must be specified via the `offset` argument.")
  }

  # the default `~ 1` is created in this frame; don't capture (and serialize) it
  if (identical(environment(dispformula), environment())) {
    environment(dispformula) <- globalenv()
  }
  disp_design <- .make_design(dispformula, data, "`dispformula`")
  Xdisp_ij <- disp_design$X

  split_formula <- list() # passed to out structure, not TMB
  X_ij <- list() # main effects, passed into TMB
  mf <- list()

  for (ii in seq_along(formula)) {
    check_offset(formula[[ii]])

    # random effects parsed into data structures (smoothers removed):
    split_formula[[ii]] <- parse_formula(remove_s_and_t2(formula[[ii]]), data)

    mf[[ii]] <- model.frame(remove_s_and_t2(reformulas::nobars(formula[[ii]])), data)
    X_ij[[ii]] <- model.matrix(attr(mf[[ii]], "terms"), mf[[ii]])
    is_inf <- vapply(mf[[ii]], function(x) any(is.infinite(x)), logical(1))
    if (any(is_inf)) {
      cli_abort("Column `{names(mf[[ii]])[is_inf][1]}` had an Inf/-Inf value. Please remove this before fitting the model.")
    }
  }

  if (has_two_components) {
    random_effects <- lapply(split_formula, function(x) {
      vapply(x$bars, safe_deparse, character(1))
    })
    if (!identical(random_effects[[1L]], random_effects[[2L]])) {
      cli_abort(c(
        "Random-effect terms must be identical in both components of a delta model.",
        "i" = "Use the same random intercepts, slopes, and grouping variables in both formulas."
      ))
    }
  }

  re <- re_tmb_data(split_formula) # random intercepts and slopes

  has_smooths <- vapply(formula, function(f) {
    length(get_smooth_terms(all_terms(reformulas::nobars(f)))) > 0L
  }, logical(1))
  if (has_two_components && any(has_smooths) &&
      original_formula[[1]] != original_formula[[2]]) {
    cli_abort(paste0(
      "For now, if delta models contain smoothers, both components ",
      "must have the same main-effects formula."
    ))
  }

  # smoothers are shared by all components, so parse them once:
  formula_no_bars <- reformulas::nobars(formula[[1]])
  sm <- parse_smoothers(formula = formula_no_bars, data = data, knots = knots)
  sm$formula_no_sm <- remove_s_and_t2(formula[[1]])
  sm$formula_no_bars <- formula_no_bars

  y_i <- model.response(mf[[1]], "any")

  # Keep the CV inclusion mask separate from user weights. In particular,
  # binomial user weights can represent trial sizes rather than likelihood
  # weights and must be processed before the mask is applied.
  cv_fold_weights <- experimental[[".cv_fold_weights"]]
  if (!is.null(cv_fold_weights) && length(cv_fold_weights) != nrow(data)) {
    cli_abort("Internal error: CV fold weights do not match the number of data rows.")
  }

  if (has_two_components) {
    y_i2 <- model.response(mf[[2]], "any")
    if (!identical(y_i, y_i2)) {
      cli_abort("Response variable should be the same in both parts of the delta formula.")
    }
  }

  response <- .prepare_family_response(y_i, weights, family_spec, censored_upper)
  y_i <- response$y_i
  size <- response$size
  weights <- response$weights
  if (!is.null(response$upr)) censored_upper <- response$upr

  likelihood_weights <- if (!is.null(weights)) weights else rep(1, NROW(y_i))
  if (!is.null(cv_fold_weights)) {
    weights <- likelihood_weights * cv_fold_weights
  }

  if (identical(family$link[1], "log") && min(y_i, na.rm = TRUE) < 0 && !has_two_components) {
    cli_abort("`link = 'log'` but the reponse data include values < 0.")
  }
  if (!has_two_components && family$family[1] %in% c("truncated_nbinom1", "truncated_nbinom2") &&
    any(y_i < 1, na.rm = TRUE)) {
    cli_abort(c(
      "`{family$family[1]}()` requires a response > 0.",
      "i" = "Consider `delta_{family$family[1]}()` for data with zeros."
    ))
  }
  if (is.null(offset)) offset <- rep(0, length(y_i))
  assert_that(length(offset) == length(y_i), msg = "Offset doesn't match length of data")
  if (nrow(Xdisp_ij) != length(y_i)) {
    cli_abort("Internal error: `dispformula` model matrix rows do not match response length.")
  }

  if (!is.null(time_varying)) {
    tv_design <- .make_design(time_varying, data, "`time_varying`")
    X_rw_ik <- tv_design$X
    if (time_varying_type == "rw") {
      X_fixed <- if (sm$has_smooths) lapply(X_ij, cbind, sm$Xs) else X_ij
      confounded <- vapply(X_fixed, rw_confounded, logical(1), X_rw = X_rw_ik)
      if (any(confounded)) {
        cli_warn(c(
          "The `time_varying` terms overlap with the fixed effects in `formula`, so the model is likely not identifiable with `time_varying_type = 'rw'`.",
          "i" = "Use `time_varying_type = 'rw0'` (or 'ar1') and keep the terms in `formula`, or remove the overlapping terms from `formula`.",
          "i" = "This includes implicit intercepts, e.g., `formula = y ~ 0 + factor(x)` with `time_varying = ~ 1`."
        ))
      }
    }
  } else {
    X_rw_ik <- matrix(0, nrow = nrow(data), ncol = 1)
    tv_design <- NULL
  }

  n_s <- domain$n_s

  # Student-t df: can be NULL (estimate, default) or numeric (fix)
  has_student_family <- any(!is.na(family_spec$param_slot$ln_student_df))
  student_df_fixed <- if (family_spec$n_f == 1L && family$family[1] == "student" && "df" %in% names(family)) {
    family$df
  } else {
    NULL
  }
  estimate_student_df <- has_student_family && is.null(student_df_fixed)

  .priors <- priors
  .priors$b <- .priors$sigma_V <- NULL # passed in as their own data
  .priors$custom <- .priors$custom_log_jacobian <- NULL # passed separately
  # Both templates apply the `b` priors to `b_j` and `b_j2` alike:
  if (has_two_components && !all(is.na(priors$b[, 1])) &&
    !identical(colnames(X_ij[[1]]), colnames(X_ij[[2]]))) {
    cli_abort(c(
      "`priors$b` applies to both delta components, so they must have the same fixed-effect columns.",
      i = "Use the same main-effects formula in both components, or drop the `b` priors.",
      i = "For priors on one component's coefficients, use `sdmTMBpriors(custom = )` with `par$b_j` or `par$b_j2` (RTMB backend)."
    ))
  }
  priors_b_tmb <- b_priors_tmb(priors$b, ncol(X_ij[[1]]))

  A_st <- domain$A_st
  A_spatial_index <- domain$A_spatial_index
  family_tmb <- .as_tmb_family_data(family_spec)
  fit_poisson_link_delta <- identical(family_spec$families$combine_kind[[1]], "poisson_link_delta") &&
    family_spec$n_f == 1L

  nonlocal_parsed <- .build_nonlocal_tmb_data(
    nonlocal_formula = nonlocal_formula_parsed,
    data = if (nonlocal_grid_supplied) nonlocal_grid_inputs$data else data,
    A_st = if (nonlocal_grid_supplied) nonlocal_grid_inputs$A_st else A_st,
    A_spatial_index = if (nonlocal_grid_supplied) nonlocal_grid_inputs$A_spatial_index else A_spatial_index,
    year_i = if (nonlocal_grid_supplied) nonlocal_grid_inputs$year_i else year_i_data,
    n_t = n_t,
    time_values = time_df$time_from_data
  )

  if (!is.null(nonlocal_parsed)) {
    for (m in seq_len(n_m)) {
      X_ij[[m]] <- .append_nonlocal_coef_columns(
        X = X_ij[[m]],
        coef_names = nonlocal_parsed$term_coef_name
      )
    }
  }
  nonlocal_tmb <- .nonlocal_tmb_slots(nonlocal_parsed)

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
    y_i = response$response,
    n_t = n_t,
    z_i = z_i,
    offset_i = offset,
    proj_offset_i = 0,
    A_st = A_st,
    spatial_model = switch(domain$spatial_model, spde = 0L, sar = 1L, car = 2L),
    W_ss = domain$W_ss,
    sar_rho_transform = 0L,
    sim_re = if ("sim_re" %in% names(experimental)) as.integer(experimental$sim_re) else rep(0L, 6),
    sim_obs = 1L,
    A_spatial_index = A_spatial_index,
    year_i = year_i_data,
    covariate_diffusion = nonlocal_tmb,
    ar1_fields = spatiotemporal == "ar1",
    simulate_t = rep(1L, n_t),
    rw_fields = spatiotemporal == "rw",
    X_ij = X_ij,
    Xdisp_ij = Xdisp_ij,
    X_rw_ik = X_rw_ik,
    Zs = sm$Zs, # optional smoother basis function matrices
    Xs = sm$Xs, # optional smoother linear effect matrix
    proj_Zs = list(),
    proj_Xs = matrix(nrow = 0L, ncol = 0L),
    b_smooth_start = sm$b_smooth_start,
    proj_lon = 0,
    proj_lat = 0,
    proj_vector = matrix(0),
    do_predict = 0L,
    do_rsr = as.integer(isTRUE(control$get_rsr)),
    calc_se = 0L,
    pop_pred = 0L,
    short_newdata = 0L,
    weights_i = if (!is.null(weights)) weights else rep(1, NROW(y_i)),
    area_i = rep(1, NROW(y_i)),
    normalize_in_r = domain$normalize_in_r,
    flag = 1L, # part of TMB::normalize()
    calc_index_totals = 0L,
    calc_cog = 0L,
    calc_eao = 0L,
    calc_weighted_avg = 0L,
    random_walk = random_walk,
    ar1_time = as.integer(!is.null(time_varying) && time_varying_type == "ar1"),
    priors_sigma_V = sigma_V_priors_tmb(priors$sigma_V, ncol(X_rw_ik)),
    priors = as.numeric(unlist(.priors)),
    share_range = as.integer(if (length(share_range) == 1L) rep(share_range, 2L) else share_range),
    include_spatial = as.integer(spatial == "on"),
    omit_spatial_intercept = as.integer(omit_spatial_intercept),
    proj_mesh = if (is_areal) dummy_sparse_1x1() else Matrix::Matrix(c(0, 0, 2:0), 3, 5), # dummy
    proj_X_ij = list(matrix(0, ncol = 1, nrow = 1)), # dummy
    proj_X_rw_ik = matrix(0, ncol = 1, nrow = 1), # dummy
    proj_year = 0, # dummy
    proj_time_include = rep(1L, n_t),
    proj_family_id = 0L, # dummy
    proj_spatial_index = 0L, # dummy
    proj_z_i = matrix(0, nrow = 1, ncol = n_m), # dummy
    indexes_total = numeric(0), # dummy
    n_integration = 0L, # dummy
    spde_aniso = domain$spde_aniso_struct,
    spde = domain$spde_struct,
    barrier = domain$barrier,
    spde_barrier = domain$spde_barrier_struct,
    barrier_scaling = domain$barrier_scaling,
    anisotropy = domain$anisotropy,
    # Used only for derived prediction quantities (e.g., indices), so that
    # `derived_link` does not alter the likelihood link.
    link_pred = family_tmb$link_code[1L, 1L],
    size = c(size),
    spatial_only = as.integer(spatiotemporal == "off"),
    spatial_covariate = as.integer(!is.null(spatial_varying)),
    calc_quadratic_range = as.integer(control$quadratic_roots),
    X_threshold = thresh$X_threshold,
    proj_X_threshold = 0, # dummy
    threshold_func = thresh$threshold_func,
    has_smooths = as.integer(sm$has_smooths),
    has_dispersion_model = as.integer(has_dispformula),
    upr = censored_upper,
    cens_direct = .censored_direct_init(weights, NROW(y_i), control$censored_method),
    lwr = 0L, # in case we want to reintroduce this
    stan_flag = as.integer(bayesian),
    no_spatial = no_spatial,
    Zt_list_proj = list(),
    exclude_RE = 0L
  )
  tmb_data <- c(tmb_data, re, priors_b_tmb)
  # R functions, so RTMB only; kept in the data list so every rebuilt
  # objective (prediction, index, simulation, saved fits) includes them
  tmb_data$priors_custom <- priors_custom
  tmb_data <- c(tmb_data, family_tmb)
  tmb_data$poisson_link_delta <- as.integer(fit_poisson_link_delta)
  family_params <- .family_parameter_values(
    family_spec,
    estimate_student_df = estimate_student_df,
    student_df_fixed = student_df_fixed
  )

  tmb_params <- list(
    ln_H_input = matrix(0, nrow = 2L, ncol = n_m),
    b_j = rep(0, ncol(X_ij[[1]])),
    b_j2 = if (has_two_components) rep(0, ncol(X_ij[[2]])) else numeric(0),
    b_disp_k = rep(0, ncol(Xdisp_ij)),
    bs = if (sm$has_smooths) matrix(0, nrow = ncol(sm$Xs), ncol = n_m) else array(0),
    ln_tau_O = rep(0, n_m),
    ln_tau_Z = matrix(0, n_z, n_m),
    ln_tau_E = rep(0, n_m),
    ln_kappa = matrix(0, 2L, n_m),
    log_kappaS_nl = .nonlocal_log_kappaS_start(spde$loc_xy, nonlocal_tmb$n_covariates),
    log_kappaT_nl = numeric(nonlocal_tmb$n_covariates),
    thetaf = family_params$thetaf,
    ln_student_df = family_params$ln_student_df,
    gengamma_Q = family_params$gengamma_Q,
    psi = family_params$psi,
    logit_p_extreme = 0,
    log_ratio_mix = -1, # ratio is 1 + exp(log_ratio_mix) so 0 would start fairly high
    ln_phi = family_params$ln_phi,
    ln_tau_V = matrix(0, ncol(X_rw_ik), n_m),
    rho_time_unscaled = matrix(0, ncol(X_rw_ik), n_m),
    ar1_phi = rep(0, n_m),
    logit_rho_sar = if (is_areal) rep(0, n_m) else numeric(0),
    re_cov_pars = matrix(0, max(c(0, re$re_cov_df_map[, "end"])) + 1, n_m),
    re_b_pars = matrix(0, nrow(re$var_indx_matrix), n_m),
    b_rw_t = array(0, dim = c(tmb_data$n_t, ncol(X_rw_ik), n_m)),
    omega_s = matrix(0, if (!omit_spatial_intercept) n_s else 0L, n_m),
    zeta_s = array(0, dim = c(n_s, n_z, n_m)),
    epsilon_st = array(0, dim = c(n_s, tmb_data$n_t, n_m)),
    b_threshold = if (thresh$threshold_func == 2L) matrix(0, 3L, n_m) else matrix(0, 2L, n_m),
    b_smooth = if (sm$has_smooths) matrix(0, sum(sm$sm_dims), n_m) else array(0),
    ln_smooth_sigma = if (sm$has_smooths) matrix(0, length(sm$sm_dims), n_m) else array(0)
  )
  if (family_spec$n_f == 1L && identical(family$link, "inverse") && family$family[1] %in% c("Gamma", "gaussian", "student") && !has_two_components) {
    fam <- family
    if (family$family == "student") fam <- stats::gaussian(link = "inverse")
    # fit the design matrix already built; NA coefficients (e.g., the all-zero
    # nonlocal columns) start at 0
    cf <- suppressWarnings(stats::glm.fit(X_ij[[1]], c(tmb_data$y_i),
      weights = c(tmb_data$weights_i), offset = c(tmb_data$offset_i), family = fam)$coefficients)
    tmb_params$b_j <- ifelse(is.na(cf), 0, unname(cf))
  }

  # Map off parameters not needed
  tmb_map <- map_all_params(tmb_params)
  tmb_map$b_j <- NULL
  if (!is_multi_family && has_dispformula) tmb_map <- unmap(tmb_map, "b_disp_k")
  if (has_two_components) tmb_map$b_j2 <- NULL
  tmb_map <- .map_family_parameters(
    tmb_map, tmb_params,
    has_dispformula = has_dispformula,
    estimate_student_df = estimate_student_df
  )

  tmb_map$log_kappaS_nl <- .make_nonlocal_kappa_map(nonlocal_parsed$covariate_has_spatial)
  tmb_map$log_kappaT_nl <- .make_nonlocal_kappa_map(nonlocal_parsed$covariate_has_temporal)

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
  if (isTRUE(domain$anisotropy == 1L)) tmb_map <- unmap(tmb_map, "ln_H_input")
  if (!is.null(time_varying)) {
    tmb_random <- c(tmb_random, "b_rw_t")
    tmb_map <- unmap(tmb_map, c("b_rw_t", "ln_tau_V"))
    if (time_varying_type == "ar1") {
      tmb_map <- unmap(tmb_map, "rho_time_unscaled")
    }
  }

  is_ar1 <- spatiotemporal == "ar1"
  tmb_map$ar1_phi <- factor(ifelse(is_ar1, cumsum(is_ar1), NA))
  if (sum(re$n_re_groups) > 0) {
    tmb_random <- c(tmb_random, "re_b_pars")
    # components share their random effects, so every row is estimated
    tmb_map$re_b_pars <- factor(seq_along(tmb_params$re_b_pars))
    tmb_map$re_cov_pars <- factor(seq_along(tmb_params$re_cov_pars))
  }
  # fixed effects that REML integrates out; logLik() still counts them
  reml_random <- NULL
  if (reml) reml_random <- c(reml_random, "b_j")
  if (reml && has_two_components) reml_random <- c(reml_random, "b_j2")
  if (reml && sm$has_smooths) reml_random <- c(reml_random, "bs")
  tmb_random <- c(tmb_random, reml_random)

  if (sm$has_smooths) {
    tmb_random <- c(tmb_random, "b_smooth") # smooth random effects
    tmb_map <- unmap(tmb_map, c("b_smooth", "ln_smooth_sigma", "bs"))
  }

  if (!is.null(previous_fit)) {
    tmb_params <- previous_fit$parlist %||% previous_fit$tmb_params
    tmb_map <- previous_fit$tmb_map
  }

  # SVC fields use the spatial range of every component:
  kappa_spatial <- if (is.null(spatial_varying)) spatial else rep("on", n_m)
  tmb_map$ln_kappa <- get_kappa_map(n_m = n_m, spatial = kappa_spatial, spatiotemporal = spatiotemporal, share_range = share_range)
  if (is_areal) {
    tmb_map$ln_kappa <- factor(rep(NA_integer_, length(tmb_params$ln_kappa)))
    areal_field_active <- any(spatial == "on" & !omit_spatial_intercept) ||
      !all(spatiotemporal == "off") ||
      !is.null(spatial_varying)
    if (areal_field_active) {
      tmb_map <- unmap(tmb_map, "logit_rho_sar")
    }
  }

  .validate_nonlocal_control(control, nonlocal_parsed)
  for (nm in names(control$start)) {
    cli_inform(c(
      i = paste0("Initiating `", nm, "` at specified starting value(s) of:"),
      paste0("  ", paste(round(control$start[[nm]], 3), collapse = ", "))
    ))
  }
  tmb_params[names(control$start)] <- control$start

  if ("ln_kappa" %in% names(control$start) &&
      (!is.matrix(tmb_params[["ln_kappa"]]) || nrow(tmb_params[["ln_kappa"]]) != 2L)) {
    msg <- c(
      "Note that `ln_kappa` must be a matrix of nrow 2 and ncol models (regular=1, delta=2).",
      "It should be the same value in each row if `share_range = TRUE`."
    )
    cli_abort(msg)
  }

  data$sdm_x <- data$sdm_y <- data$sdm_orig_id <- data$sdm_spatial_id <- NULL
  data$sdmTMB_X_ <- data$sdmTMB_Y_ <- NULL

  if (has_two_components && "off" %in% spatiotemporal) {
    tmb_map$epsilon_st <- map_off_components(tmb_params$epsilon_st, spatiotemporal == "off")
    tmb_map$ln_tau_E <- map_off_components(tmb_params$ln_tau_E, spatiotemporal == "off")
  }
  if (has_two_components && "off" %in% spatial) {
    tmb_map$omega_s <- map_off_components(tmb_params$omega_s, spatial == "off")
    tmb_map$ln_tau_O <- map_off_components(tmb_params$ln_tau_O, spatial == "off")
  }

  if (isTRUE(domain$anisotropy == 1L) && has_two_components) {
    tmb_map$ln_H_input <- factor(c(1, 2, 1, 2)) # share anisotropy as in VAST
  }

  if (family_spec$n_f == 1L && family$family[[n_m]] %in% rtmb_mixture_families) {
    p_extreme <- (if (has_two_components) family[[2]]$p_extreme) %||% family$p_extreme
    tmb_map$log_ratio_mix <- NULL
    if (is.null(p_extreme)) {
      tmb_map$logit_p_extreme <- NULL
    } else {
      tmb_map$logit_p_extreme <- factor(NA)
      tmb_params$logit_p_extreme <- stats::qlogis(p_extreme)
    }
  }

  if (tmb_data$threshold_func > 0) tmb_map$b_threshold <- NULL

  for (nm in names(control$map)) cli_inform(c(i = paste0("Fixing or mirroring `", nm, "`")))
  tmb_map[names(control$map)] <- control$map

  # profile_breakpt() maps b_threshold to fix the cutpoint, so it doesn't trigger this
  if (do_fit && identical(as.integer(tmb_data$threshold_func), 1L) && is.null(control$map$b_threshold)) {
    cli_inform(class = "sdmTMB_breakpt_message", c(
      i = paste0("The likelihood of a `breakpt()` model isn't smooth in the cutpoint, ",
        "so the fit can stop at the wrong cutpoint while reporting convergence, ",
        "and the cutpoint's standard error isn't reliable."),
      i = breakpt_advice(has_two_components)
    ))
  }

  # Phase 1 runs on the full map and starting values, so user `map` and
  # `start` are respected; the user's `start` values always win
  if (control$multiphase && is.null(previous_fit) && do_fit) {
    tmb_params <- fit_first_phase(tmb_data, tmb_params, tmb_map,
      profile = control$profile, backend = backend, lower = control$lower, upper = control$upper,
      mesh = if (is_areal) NULL else spde$mesh, nlminb_control = .control,
      silent = silent, suppress_warnings = isTRUE(control$suppress_nlminb_warnings))
    tmb_params[names(control$start)] <- control$start
  }

  out_structure <- structure(
    list(
      data = data,
      offset = offset,
      spde = spde,
      formula = original_formula,
      dispformula = dispformula,
      has_dispformula = has_dispformula,
      split_formula = split_formula,
      time_varying = time_varying,
      threshold_parameter = thresh$threshold_parameter,
      threshold_function = thresh$threshold_func,
      time = time,
      time_lu = time_df,
      # Keep the public field faithful to the user's input. Internal code uses
      # `family_spec`, so a multi-family fit must not expose only family one.
      family = family_spec$family_input,
      family_spec = family_spec,
      distribution_column = family_spec$distribution_column,
      likelihood_weights = likelihood_weights,
      smoothers = sm,
      response = response$response,
      tmb_data = tmb_data,
      tmb_params = tmb_params,
      tmb_map = tmb_map,
      tmb_random = tmb_random,
      reml_random = reml_random,
      backend = backend,
      spatial_varying = spatial_varying,
      nonlocal_formula = nonlocal_formula,
      nonlocal_formula_parsed = nonlocal_formula_parsed,
      nonlocal_parsed = nonlocal_parsed,
      nonlocal_grid_supplied = nonlocal_grid_supplied,
      spatial = spatial_user,
      spatiotemporal = spatiotemporal,
      spatial_varying_formula = spatial_varying_formula,
      svc_omega_is_intercept = svc_omega_is_intercept,
      reml = reml,
      priors = priors,
      nlminb_control = .control,
      control = control,
      contrasts = lapply(X_ij, attr, which = "contrasts"),
      # encodings of auxiliary designs for prediction (see `.apply_design()`)
      design_specs = list(
        spatial_varying = svc_design$spec,
        time_varying = tv_design$spec,
        dispformula = disp_design$spec
      ),
      terms = lapply(mf, attr, which = "terms"),
      extra_time = extra_time,
      fitted_time = sort(unique(data[[time]])),
      xlevels = lapply(mf, function(x) stats::.getXlevels(attr(x, "terms"), x)),
      call = match.call(expand.dots = TRUE),
      version = utils::packageVersion("sdmTMB"),
      experimental = experimental
    ),
    class = "sdmTMB"
  )

  tmb_obj <- make_sdmTMB_adfun(
    data = tmb_data, parameters = tmb_params, map = tmb_map,
    profile = control$profile,
    random = tmb_random, backend = backend, silent = silent
  )
  lim <- set_limits(tmb_obj,
    lower = control$lower, upper = control$upper,
    mesh = if (is_areal) NULL else spde$mesh,
    spatial_model = tmb_data$spatial_model,
    silent = FALSE
  )

  out_structure$tmb_obj <- tmb_obj
  out_structure$lower <- lim$lower
  out_structure$upper <- lim$upper

  if (!do_fit) {
    return(out_structure)
  }

  if (isTRUE(domain$normalize_in_r == 1L) && control$normalize) {
    tmb_obj <- TMB::normalize(tmb_obj, flag = "flag", value = 0)
  }

  if (length(tmb_obj$par)) {
    tmb_opt <- maybe_suppress_warnings(control$suppress_nlminb_warnings)(stats::nlminb(
      start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
      lower = lim$lower, upper = lim$upper, control = .control
    ))
  } else {
    tmb_opt <- list(par = tmb_obj$par, objective = tmb_obj$fn(tmb_obj$par))
  }

  tmb_opt <- run_nlminb_loops(
    nlminb_loops = control$nlminb_loops - 1, opt = tmb_opt, obj = tmb_obj,
    lower = lim$lower, upper = lim$upper, control = .control,
    silent = silent, suppress_warnings = isTRUE(control$suppress_nlminb_warnings)
  )

  if ("censored_betabinomial" %in% family_spec$components$family_name) {
    refit <- refit_censored_betabinomial(tmb_obj, tmb_opt, tmb_data, tmb_map,
      tmb_random, lim, .control, profile = control$profile, silent = silent,
      suppress_warnings = isTRUE(control$suppress_nlminb_warnings))
    tmb_obj <- out_structure$tmb_obj <- refit$obj
    tmb_data <- out_structure$tmb_data <- refit$data
    tmb_opt <- refit$opt
  }

  check_bounds(tmb_opt$par, lim$lower, lim$upper)

  # Check fields before Newton steps to avoid Hessian issues after collapse
  if ((control$collapse_spatial_variance || control$collapse_spatiotemporal_ar1) &&
      length(tmb_obj$par) > 0) {
    report_vals <- tmb_obj$report()
    collapse_result <- check_and_collapse_spatial_fields(
      report_vals = report_vals,
      spatial = spatial,
      spatiotemporal = spatiotemporal,
      n_m = n_m,
      n_t = n_t,
      omit_spatial_intercept = omit_spatial_intercept,
      collapse_spatial_variance = control$collapse_spatial_variance,
      collapse_spatial_variance_threshold = control$collapse_spatial_variance_threshold,
      collapse_spatiotemporal_ar1 = control$collapse_spatiotemporal_ar1,
      collapse_ar1_threshold = control$collapse_ar1_threshold,
      has_two_components = has_two_components,
      silent = silent
    )
    if (collapse_result$do_refit) {
      if (!silent) {
        cli_inform(c(
          "i" = "Refitting model with collapsed spatial field(s)..."
        ))
      }

      # Re-run the original call with only the field arguments changed, so
      # missing-value rows are omitted once, from the original inputs
      call <- out_structure$call
      call$spatial <- collapse_result$spatial_arg
      call$spatiotemporal <- collapse_result$spatiotemporal_arg
      for (nm in names(refit_inputs)) call[nm] <- list(refit_inputs[[nm]])
      return(eval(call, envir = parent.frame()))
    }
  }

  # We only end up here if no collapse was detected
  tmb_opt <- run_newton_loops(
    newton_loops = control$newton_loops, tmb_opt, tmb_obj, silent,
    lower = lim$lower, upper = lim$upper
  )

  if (!silent && control$getsd) cli_inform("running TMB sdreport\n")
  set_fit_results(out_structure, tmb_obj, tmb_opt, getsd = control$getsd,
    get_joint_precision = control$get_joint_precision)
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

check_and_collapse_spatial_fields <- function(
    report_vals,
    spatial,
    spatiotemporal,
    n_m,
    n_t,
    omit_spatial_intercept,
    collapse_spatial_variance,
    collapse_spatial_variance_threshold,
    collapse_spatiotemporal_ar1,
    collapse_ar1_threshold,
    has_two_components = FALSE,
    silent = FALSE) {

  do_refit <- FALSE
  spatial_updated <- spatial
  spatiotemporal_updated <- spatiotemporal

  if (collapse_spatial_variance && any(spatial == "on") &&
      !omit_spatial_intercept) {
    est_sigma_O <- c(report_vals$sigma_O)
    which_sigma <- which(
      spatial == "on" & is.finite(est_sigma_O) &
        est_sigma_O < collapse_spatial_variance_threshold
    )

    if (length(which_sigma) > 0L) {
      if (!silent) {
        cli_inform(c(
          "!" = "Spatial variance below threshold ({collapse_spatial_variance_threshold}) detected",
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

  if (collapse_spatial_variance && !all(spatiotemporal == "off")) {
    est_sigma_E <- report_vals$sigma_E

    if (length(est_sigma_E) > 0) {
      est_sigma_E_by_model <- rep(NA_real_, n_m)

      for (m in seq_len(n_m)) {
        idx_start <- (m - 1) * n_t + 1
        idx_end <- m * n_t
        if (idx_end <= length(est_sigma_E)) {
          est_sigma_E_by_model[m] <- min(est_sigma_E[idx_start:idx_end])
        }
      }

      which_sigma <- which(
        spatiotemporal != "off" & is.finite(est_sigma_E_by_model) &
          est_sigma_E_by_model < collapse_spatial_variance_threshold
      )

      if (length(which_sigma) > 0L) {
        if (!silent) {
          cli_inform(c(
            "!" = "Spatiotemporal variance below threshold ({collapse_spatial_variance_threshold}) detected",
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

  if (collapse_spatiotemporal_ar1 &&
      any(spatiotemporal_updated == "ar1")) {
    reported_rho <- c(report_vals$rho)
    est_rho <- rep(NA_real_, n_m)
    n_rho <- min(length(reported_rho), n_m)
    if (n_rho > 0L) {
      est_rho[seq_len(n_rho)] <- reported_rho[seq_len(n_rho)]
    }

    which_iid <- which(
      spatiotemporal_updated == "ar1" & is.finite(est_rho) &
        abs(est_rho) <= collapse_ar1_threshold
    )
    which_rw <- which(
      spatiotemporal_updated == "ar1" & is.finite(est_rho) &
        est_rho >= 1 - collapse_ar1_threshold
    )

    if (length(which_iid) > 0L) {
      if (!silent) {
        cli_inform(c(
          "!" = "Spatiotemporal AR1 correlation within {collapse_ar1_threshold} of zero detected",
          "i" = "Affected model(s): {paste(which_iid, collapse = ', ')}",
          ">" = "Refitting with spatiotemporal field(s) set to IID"
        ))
      }
      spatiotemporal_updated[which_iid] <- "iid"
      do_refit <- TRUE
    }

    if (length(which_rw) > 0L) {
      if (!silent) {
        cli_inform(c(
          "!" = "Spatiotemporal AR1 correlation within {collapse_ar1_threshold} of one detected",
          "i" = "Affected model(s): {paste(which_rw, collapse = ', ')}",
          ">" = "Refitting with spatiotemporal field(s) set to RW"
        ))
      }
      spatiotemporal_updated[which_rw] <- "rw"
      do_refit <- TRUE
    }
  }

  if (has_two_components) {
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

set_limits <- function(tmb_obj, lower, upper, mesh = NULL, spatial_model = 0L,
                       silent = TRUE) {
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
  if ("logit_rho_sar" %in% names(tmb_obj$par) &&
    !"logit_rho_sar" %in% union(names(lower), names(upper))) {
    if (identical(spatial_model, 2L)) {
      .upper["logit_rho_sar"] <- stats::qlogis(0.999)
    } else {
      .lower["logit_rho_sar"] <- stats::qlogis((-0.999 + 1) / 2)
      .upper["logit_rho_sar"] <- stats::qlogis((0.999 + 1) / 2)
    }
  }
  is_kappaS_nl <- names(.lower) == "log_kappaS_nl"
  if (any(is_kappaS_nl) && !is.null(mesh)) {
    bounds <- .nonlocal_log_kappaS_bounds(mesh)
    if (!"log_kappaS_nl" %in% names(lower)) .lower[is_kappaS_nl] <- bounds[[1L]]
    if (!"log_kappaS_nl" %in% names(upper)) .upper[is_kappaS_nl] <- bounds[[2L]]
  }

  list(lower = .lower, upper = .upper)
}

check_spatiotemporal_arg <- function(x, time, .which = 1) {
  .spatiotemporal <- tolower(as.character(x[[.which]]))
  assert_that(.spatiotemporal %in% c("iid", "ar1", "rw", "off", "true", "false"),
    msg = "`spatiotemporal` argument value not valid"
  )
  if (.spatiotemporal == "true") .spatiotemporal <- "iid"
  if (.spatiotemporal == "false") .spatiotemporal <- "off"
  .spatiotemporal <- match.arg(tolower(.spatiotemporal), choices = c("iid", "ar1", "rw", "off"))
  if (is.null(time) && .spatiotemporal != "off") {
    cli_abort("`spatiotemporal` is set but the `time` argument is missing.")
  }
  .spatiotemporal
}

map_all_params <- function(x) {
  m <- list()
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

# With `time_varying_type = "rw"` the first time slice is unconstrained, so any
# combination of `time_varying` columns that lies in the span of the fixed
# effects (e.g., `~ 0 + factor(x)` with `time_varying = ~ 1`) is not identified.
rw_confounded <- function(X, X_rw) {
  if (!ncol(X) || !ncol(X_rw)) return(FALSE)
  rank <- function(x) qr(x)$rank
  rank(cbind(X, X_rw)) < rank(X) + rank(X_rw)
}

check_irregalar_time <- function(data, time, spatiotemporal, time_varying, extra_time,
                                 nonlocal_temporal = FALSE) {
  has_ar1_rw <- any(spatiotemporal %in% c("ar1", "rw")) || !is.null(time_varying)
  has_nl_temporal <- isTRUE(nonlocal_temporal)
  if (!(has_ar1_rw || has_nl_temporal)) return(invisible(NULL))

  if (!is.numeric(data[[time]])) {
    cli_abort("Time column should be integer or numeric if using AR(1), random walk, or temporal nonlocal processes.")
  }

  ti_data <- sort(unique(data[[time]]))
  ti <- sort(union(ti_data, extra_time))
  irregular_fit <- length(unique(diff(ti))) > 1L
  irregular_nl <- has_nl_temporal && length(unique(diff(ti_data))) > 1L

  if (irregular_fit || irregular_nl) {
    missed <- find_missing_time(data[[time]])
    msg <- c(
      "Detected irregular time spacing with an AR(1), random walk, or time-indexed nonlocal process.",
      if (has_ar1_rw) {
        "Consider filling in the missing time slices with `extra_time`."
      },
      if (irregular_nl) {
        c(
          "For time-indexed `nonlocal_formula` terms, include rows in `data` for missing time slices so nonlocal covariates are defined.",
          "Using `extra_time` alone is not sufficient for time-indexed nonlocal terms."
        )
      },
      if (length(missed)) {
        paste0("`extra_time = c(", paste(missed, collapse = ", "), ")`")
      }
    )
    if (irregular_nl) cli_abort(msg) else cli_inform(msg)
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

# Map off the slices of parameter `x` (last dimension = model component) for
# the components in logical vector `off`
map_off_components <- function(x, off) {
  m <- matrix(seq_along(x), ncol = length(off))
  m[, off] <- NA
  factor(m)
}

get_spde_matrices <- function(x) {
  x <- x$spde[c("c0", "g1", "g2")]
  names(x) <- c("M0", "M1", "M2") # legacy INLA names needed!
  x
}
