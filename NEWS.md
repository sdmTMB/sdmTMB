# sdmTMB 0.8.1

## Minor improvements and fixes

* Refactor `get_cog()`, `get_weighted_average()`, and `get_eao()`. 
  `get_cog()` should be more memory efficient and bias correction is now
  enabled for `get_eao()` and `get_weighted_average()`. These all now work
  with subsets of the full set of fitted time elements as `get_index()` has
  for a while.

* Optimize memory usage and speed of `get_range_edge()`.

* Remove extraDistr import to avoid CRAN warning. 

* Switch from `lme4` to `reformulas` package for random effects parsing functions
  (`findbars()`, `nobars()`, `subbars()`, `mkReTrms()`). The `lme4` maintainers
  moved these formula utilities to the separate `reformulas` package.

* Fix `predict()` with `se_fit = TRUE` for delta models to correctly return
  standard errors on combined predictions when `model = NA`. Previously,
  standard errors were only calculated for the first component (`est1`) instead
  of the combined estimate (`est`). #496

* Fix typo in `predict()` for delta models that incorrectly calculated combined
  predictions as `p1 * p1` instead of `p1 * p2` in the specific case of
  `re_form = NA` and `type = "link"`. #496

* Dramatically improve `print()` performance for models with large embedded
  objects in calls (e.g., cross-validation scenarios with replicated data).
  The fix avoids expensive `deparse()` operations on large data frames and
  mesh objects by extracting only the argument names needed for display. #492

* Fix error when parsing model when including `extra_time` in delta-family models
  with a list formula. #493

# sdmTMB 0.8.0

## New features

* `sdmTMB_cv()` now supports the `weights` argument. User-supplied weights are
  combined with the internal fold-assignment mechanism (held-out data are
  assigned weight 0). Weights must be positive (> 0). #486

* Add experimental `collapse_spatial_variance` option to `sdmTMBcontrol()` to
  automatically detect and turn off spatial and/or spatiotemporal random fields
  when their SD parameters are estimated to be very small (below a threshold).
  When enabled, the model will automatically refit with the appropriate fields
  turned off if SD parameters fall below `collapse_threshold` (default 0.01).
  This can help avoid boundary issues and improve model stability when random
  fields are not needed. Set to `FALSE` by default. #263

* Add new experimental function `get_range_edge()` to calculate range edges as
  density-weighted quantiles along a spatial axis (e.g., latitude, coastal distance).
  Range edges are calculated as positions where cumulative density equals specified
  quantiles. Uses linear interpolation for accurate
  quantile estimation and simulation from the joint precision matrix for uncertainty
  quantification. Implements a similar approach to VAST range edge calculations
  following Fredston et al. (2021) <https://doi.org/10.1111/gcb.15614>.
  See the new range edge vignette at <https://sdmTMB.github.io/sdmTMB/articles/>

* The Student-t degrees of freedom parameter is now
  estimated by default in `student()`. Previously it was fixed at 3. To fix it
  at a specific value, set `df` to a numeric value (e.g., `student(df = 3)`).
  To estimate it (new default), set `df = NULL` or omit the argument. The
  parameter is constrained to be > 1. The function now prints an informative
  message about the df parameter behavior.

* Add deviance residuals for left-out data in the cross validation output.
  This can be used to calculate deviance explained of left-out data.
  See `cv_deviance_resid` column in the `data` element of the output of
  `sdmTMB_cv()`.

* Add `save_models` argument to `sdmTMB_cv()` (defaults to `TRUE`). Setting
  `save_models = FALSE` prevents storing fitted model objects for each fold,
  which can substantially reduce memory usage for large datasets or many folds.
  When `FALSE`, functions requiring model access (`tidy()`, `cv_to_waywiser()`)
  will error with informative messages. Essential CV metrics (predictions, log
  likelihoods, deviance residuals, convergence info) remain available.

* Extend `sdmTMB_simulate()` to support time-varying effects with vector `sigma_V`
  inputs and AR1 correlation (`rho_time`). #447

* Add new function `cv_to_waywiser()` to convert cross-validation results to sf
  format for use with the waywiser package. This enables multi-scale spatial
  assessment of model predictions. #193

* Add vignette demonstrating how to fit zero-one-inflated beta (ZOIB) models
  by fitting three separate model components and combining predictions. #440 #441

* Add argument to fix probability of extreme events for `*_mix()` families.
  Note that the internal parameter name has also changed from `p_mix` to
  `p_extreme` and from `logit_p_mix` to `logit_p_extreme`. #318 #474

* Add beta-binomial family (`betabinomial()`) for modeling overdispersed binomial
  data. Supports logit and cloglog links, and includes proper residuals support.

* Add new function `get_weighted_average()` to calculate biomass-weighted averages
  of user-supplied vectors (e.g., depth, temperature). This function follows the
  same pattern as `get_cog()` but allows users to specify any variable for
  weighted averaging. Supports bias correction and area weighting.

* Add vignette on age (or length) composition standardization. To help with this
  add a new experimental function `make_category_svc()`.

* Add `emmeans` support for delta/hurdle models. Use `model = 1` for the binomial
  component or `model = 2` for the positive component when calling `emmeans()`.
  Example: `emmeans(fit, ~ predictor, model = 2)`. #247 #249

## Minor improvements and fixes

* **Fix barrier model implementation**. The SPDE input matrices for the barrier
  model from INLA and INLAspacetime had changed. sdmTMB now appropriately
  uses these new matrices and unit tests in sdmTMBextra should catch such a
  change in the future. #457

* Fix Student-t deviance residuals, which were incorrectly returning `NaN`s.

* Fix `sign()` utility function to avoid `NaN` when `x = 0`. Now returns
  standard mathematical sign function behavior: `sign(0) = 0`.

* Fix bug in `sdmTMB_cv()` automatic fold generation that could result in
  unbalanced folds with duplicate and missing fold IDs. The bug was most severe
  with large `k_folds` values (e.g., leave-one-out cross-validation with
  `k_folds = nrow(data)`), which could cause errors when folds had no data.
  User-supplied `fold_ids` were not affected.

* Add check if newdata has been filtered after prediction and before passing
  to a `get_*()` function.

* Fix `tidy()` to only include the `model` column for delta models. For non-delta
  models, the `model` column is no longer included in the output for
  `effects = "ran_pars"` and `effects = "ran_vals"`, making the output cleaner
  and more consistent.

* Update package logo.

* Add residuals for truncated negative binomial families. #481
  Thanks to @Joseph-Barss

* Fix an issue with residuals for delta models by consistently using `get_par()`,
  and another issue specifically for delta truncated negative binomial models by
  replacing NaN with NA. #484

* Fix `tidy()` with `effects = "ran_pars"` to report min/max anisotropic ranges
  (e.g., `range_min`, `range_max`) for models fit with `anisotropy = TRUE`,
  matching the values shown in `print_anisotropy()`. Standard errors and
  confidence intervals are set to NA since uncertainty in both the range parameter
  and H matrix cannot be easily propagated.

* Fix issue with ggeffects with multiple smoothers + offsets. #450

* Improve `t2()` printing and appearance in `tidy.sdmTMB()`. #415 #472

* Fix `emmeans` support for models with smoothers (`s()` terms). Previously,
  `emmeans` would fail with "Non-conformable elements in reference grid" when
  smoothers were included in the model formula.

# sdmTMB 0.7.4

## Minor improvements and fixes

* Let `simulate.sdmTMB()` work with binomial GLMs with size specified via
  `weights` and `newdata` supplied. #465

* Fix issue with fold logic in LFO (leave-future-out) cross validation 
  for `lfo_forecast > 1`. #454 Thanks to @Joseph-Barss.

* Add `update.sdmTMB()` so that the mesh argument doesn't have to be 
  specified if model is loaded in a fresh session. #461

* Change default in `get_index()` etc. to `bias_correct = TRUE`. This is the
  recommended setting for final inference and speed improvements within TMB
  have made it more viable to include the bias correction by default. #456

* Only retain Newton update parameters if they improve the objective function.
  #455

* Only run Newton updates if maximum absolute gradient is `>= 1e-9` to save
  time. #455
  
* Suppress `nlminb()` warnings by default, which can usually be ignored by the
  user and may be confusing. This can be controlled via
  `sdmTMB(..., control = sdmTMBcontrol(suppress_nlminb_warnings = FALSE))`.
  This option now mirrors tinyVAST.

* Round time-varying AR(1) rho to 2 decimals in model printing/summary.

# sdmTMB 0.7.2

## New features

* Add deviance residuals (`residuals(fit, type = "deviance")`) and
  `deviance.sdmTMB()` method (`deviance(fit)`). Proportion deviance explained
  can be calculated as `1 - deviance(fit) / deviance(fit_null)` where 
  `fit_null` is a null model, e.g., fit with `formula = ~ 1` and turning 
  off any random fields as desired
  (e.g., `spatial = "off", spatiotemporal = "off"`).

* Add `observation_error` argument to `simulate.sdmTMB()` to allow
  turning off observation error simulation. The intended use-case
  is for simulating from random effects but not adding observation
  error. #431

## Minor improvements and fixes

* Change the default in `dharma_residuals()` to
  `test_uniformity = FALSE`. Based on simulation testing, we
  generally do not recommend using these p-values to reject models.

* Fix a bug introduced in version 0.7.0 where printing of the 2nd
  linear predictor smoother fixed effects (`bs`) was accidentally a copy
  of the 1st linear predictor smoother fixed effects.

* Fix bug in simulation with time-varying AR(1) when using the 
  `project()` function. Thanks to A. Allyn for pointing out the bug.

* Fix reporting of converged models with `sdmTMB_cv()`. A recent change
  resulted in reporting only 1 model converged if all models converged.

* Remove warning about old default residuals type.

* Fix `project()` and `simulate.sdmTMB(..., newdata = ...)` when 
  random intercepts/slopes are present. #431
  
* Remove extra TMB data slots for `project()` and
  `simulate.sdmTMB(..., newdata = ...)` to save memory. #431

# sdmTMB 0.7.0

## New features

* Add option for random slopes, or random intercepts to be passed in in 
  `lme4` style formulas, `density ~ (1 | fyear)` or `density ~ (depth | fyear)`,
  Matches output of `lme4` and `glmmTMB`, and summarizes output with `tidy()`.

* Add `project()` experimental function.

* Add `get_eao()` to calculate effective area occupied.

* Allow predicting on new data with `t2()` smoothers. #413

* Add priors for `breakpt()` and `logistic()` parameters. #403

* Add priors on time-varying SD parameters (`sigma_V`).

* Add `cAIC()` for calculating *conditional* AIC. Theory based on
  <https://arxiv.org/abs/2411.14185>; also see
  <https://doi.org/10.1002/ecy.4327>. J.T. Thorson wrote the function code.
  EDF (effective degrees of freedom) will ultimately be further split
  (e.g., split by smoothers) and added to `summary.sdmTMB()`. #383 #387

* Add EDF (effective degrees of freedom) printing to smoothers with 
  `print.sdmTMB()` and `summary.sdmTMB()`. Set argument `edf = TRUE`.
  E.g. `print(fit, edf = TRUE)`. #383 #387

* At experimental function `get_index_split()`, which takes care of 
  splitting a prediction grid by time, undoing the prediction and 
  area-integration index calculations for each chunk to save memory.

* Add `newdata` argument to `simulate.sdmTMB()`. This enables simulating on
  a new data frame similar to how one would predict on new data.

* Add `mle_mvn_samples` argument to `simulate.sdmTMB()`. Defaults to "single".
  If "multiple", then a sample from the random effects is taken for each
  simulation iteration.

* Allow for specifying only lower or upper limits. #394

* `sdmTMB_cv()` gains a `tidy()` and `print()` method for output. #319

* `simulate.sdmTMB()` method now has an `return_tmb_report` argument.

## New vignettes/articles

* Add forecasting and presence-only article vignettes. See
  <https://sdmTMB.github.io/sdmTMB/articles/>

* Add vignette on multispecies models with sdmTMB (or any case where one wants
  additional spatial and or spatiotemporal fields by some group).
  See <https://sdmTMB.github.io/sdmTMB/articles/>

## Minor improvements and fixes

* Add a useful error if memory error occurs on index calculation.

* Fix bug in a check in `make_mesh()` around if coordinates look
  overly large. #427

* Re-enable bias correction for `get_cog()` (get center of gravity).

* Add check for `Inf`/`-Inf` values before fitting. #408

* Add linear component of smoothers to `tidy()`. #90

* Add time varying AR(1) correlation to `tidy()` and `print()`. #374

* Warn if parameter limits are set with `newton_loops > 0`. #394

* Fix bug in `est` column when predicting on new data with Poisson-link
  delta models with `type = "link"` and `re_form = NA`. #389

* Fix bug in `s95` parameter reporting from the `tidy()` method. `s95` is
  present in the logistic threshold models. The model itself was fine but the
  `s95` parameter was supposed to be reported by `tidy()` as a combination of two
  other parameters. This also affected the output in `print()`/`summary()`.

* Add progress bar to `simulate.sdmTMB()`. #346

* Add AUC and TSS examples to cross validation vignette. #268

* Add `model` (linear predictor number) argument to `coef()` method. Also,
  write documentation for `?coef.sdmTMB`. #351

* Add helpful error message if some coordinates in `make_mesh()` are `NA`. #365

* Add informative message if fitting with an offset but predicting with offset
  argument left at `NULL` on `newdata`. #372

* Fix passing of `offset` argument through in `sdmTMB_cv()`. Before it was being
  omitted in the prediction (i.e., set to 0). #372

* Fig bug in `exponentiate` argument for `tidy()`. Set `conf.int = TRUE` as
  default. #353

* Fix bug in prediction from `delta_truncated_nbinom1()` and 
  `delta_truncated_nbinom2()` families. The positive component
  needs to be transformed to represent the mean of the *un*truncated
  distribution first before multiplying by the probability of a non-zero.
  Thanks to @tom-peatman #350

* Add option for `area` to be passed in as the name of a column in the 
  data frame to be used for area weighting. Used in `get_index()`, 
  `get_cog()`, `get_eao()`, etc.
  
# sdmTMB 0.6.0

* Pass several arguments to `DHARMa::plotQQunif()`.

* Add `silent` option in `simulate.sdmTMB()`. Setting it to `FALSE` allows
  monitoring simulations from larger models.

* Fix bug in `est_non_rf1` and `est_non_rf2` columns when all the following
  conditions were true:
  - predicting on new data
  - using a delta model
  - including IID random intercepts or time-varying coefficients
  See #342. Thanks to @tom-peatman for the issue report.

* Fix delta-gamma binomial link printing for `type = 'poisson-link'` #340

* Add suggestion to use an optimized BLAS library to README.

* Add warning if it's detected that there were problems reloading (e.g., with
  `readRDS()`) a fitted model. Simultaneously revert the approach to 
  how reloaded models are reattached.

* Move `log_ratio_mix` parameter to 2nd phase with starting value of -1 instead
  of 0 to improve convergence.

* Fix bugs for `nbinom1()` and `nbinom2_mix()` simulation.

* Allow `profile` argument in the control list to take a character vector of
  parameters. This move these parameters from the outer optimization problem to
  the inner problem (but omits from the from the Laplace approximation). See
  documentation in TMB. This can considerably speed up fitting models with many
  fixed effects.
  
* Add theoretical quantile residuals for the generalized gamma distribution.
  Thanks to J.C. Dunic. #333 
  
* Add `"poisson-link"` option to delta-mixture lognormal.

* Fix bug in simulation from Poisson-link delta models.

* Simplify the internal treatment of extra time slices (`extra_time`). #329
  This is much less bug prone and also fixes a recently introduced bug. #335
  This can slightly affect model results compared to the previous approach if
  extra time was used along with smoothers since the 'fake' extra data
  previously used was included when mgcv determined knot locations for
  smoothers.

# sdmTMB 0.5.0

* Overhaul residuals vignette ('article') 
  <https://sdmTMB.github.io/sdmTMB/articles/residual-checking.html>
  including brief intros to randomized quantile residuals, simulation-based
  residuals, 'one-sample' residuals, and uniform vs. Gaussian residuals.

* Add check if prediction coordinates appear outside of fitted coordinates. #285

* Fix memory issue with Tweedie family on large datasets. #302

* Add experimental option to return standard normal residuals from 
  `dharma_residuals()`.

* Make `simulate.sdmTMB()` not include `extra_time` elements.

* Improved re-initialization of saved fitted model objects in new sessions.

* Fix important bug in `simulate.sdmTMB()` method for delta families where
  the positive linear predictor was only getting simulated for observations
  present in the fitted data.

* Add new `"mle-mvn"` type to `residuals.sdmTMB()` and make it the default.
  This is a fast option for evaluating goodness of fit that should be better
  than the previous default. See the details section in `?residuals.sdmTMB`
  for details. The previous default is now called `"mvn-eb"` but is not
  recommended.
  
* Bring `dharma_residuals()` back over from sdmTMBextra to sdmTMB. Add a new
  option in the `type` argument (`"mle-mvn"`) that should make the
  simulation residuals consistent with the expected distribution.
  See the same new documentation in `?residuals.sdmTMB`. The examples
  in `?dharma_residuals` illustrate suggested use.

* Fix bug in `sanity()` where gradient checks were missing `abs()` such that
  large negative gradients weren't getting caught. #324

* Return `offset` vector in fitted object as an element. Ensure any extra time
  rows of data in the `data` element of the fitted object do not include the
  extra time slices.

* Add experimental residuals option "mle-mvn" where a single approximate 
  posterior sample of the random effects is drawn and these are combined
  with the MLE fixed effects to produce residuals. This may become the
  default option.

* Add the generalized gamma distribution (thanks to J.T. Thorson with additional
  work by J.C. Dunic.) See `gengamma()`. This distribution is still in a testing
  phase and is not recommended for applied use yet. #286
  
* Detect possible issue with factor(time) in formula if same column name is used
  for `time` and `extra_time` is specified. #320

* Improve `sanity()` check output when there are NA fixed effect standard
  errors.

* Set `intern = FALSE` within index bias correction, which seems to be
  considerably faster when testing with most models.

# sdmTMB 0.4.3

* Fix a bug likely introduced in July 2023 that caused issues when
  `extra_time` was specified. This is an important bug and models fit with
  `extra_time` between that date (if using the GitHub version) and v0.4.2.9004
  (2024-02-24) should be checked against a current version of sdmTMB
  (v0.4.2.9005 or greater). On CRAN, this affected v0.4.0 (2023-10-20) to 
  v0.4.2. Details:
  
  * The essence of the bug was that `extra_time` works by padding the data
    with a fake row of data for every extra time element (using the first row of
    data as the template). This is supposed to then be omitted from the 
    likelihood so it has no impact on model fitting beyond spacing
    time-series processes appropriately and setting up internal structures for
    forecasting. Unfortunately, a bug was introduced that caused these fake data
    (1 per extra time element) to be included in the likelihood.

* Issue error if `time` column has NAs. #298 #299

* Fix bug in `get_cog(..., format = "wide")` where the time column was
  hardcoded to `"year"` by accident.

* Poisson-link delta models now use a `type` argument in `delta_gamma()` and
  `delta_lognormal()`. `delta_poisson_link_gamma()` and
  `delta_poisson_link_lognormal()` are deprecated. #290
  
* Delta families can now pass links that are different from the default 
  `"logit"` and `"log"`. #290

# sdmTMB 0.4.2

* Force rebuild of CRAN binaries to fix issue with breaking Matrix ABI change
  causing `NaN gradient` errors. #288 #287

* Fix crash in if `sdmTMB(..., do_index = TRUE)` and `extra_time` supplied along
  with `predict_args = list(newdata = ...)` that lacked `extra_time` elements.

* Allow `get_index()` to work with missing time elements.

* Add the ability to pass a custom randomized quantile function `qres_func`
  to `residuals.sdmTMB()`.

* Add check for factor random intercept columns in `newdata` to avoid a crash.
  #278 #280

* Improve warnings/errors around use of `do_index = TRUE` and `get_index()`
  if `newdata = NULL`. #276

* Fix prediction with `offset` when `newdata` is `NULL` but `offset` is
  specified. #274

* Fix prediction failure when both `offset` and `nsim` are provided and
  model includes `extra_time`. #273

# sdmTMB 0.4.1

* Fix memory issues detected by CRAN 'Additional issues' clang-UBSAN, valgrind.

* Fix a bug predicting on new data with a specified offset and `extra_time`. 
  #270

* Add warning around non-factor handling of the `spatial_varying` formula. #269

* Add experimental `set_delta_model()` for plotting delta models with
  `ggeffects::ggpredict()` (GitHub version only until next CRAN version).

# sdmTMB 0.4.0

* Move add_barrier_mesh() to sdmTMBextra to avoid final INLA dependency.
  https://github.com/sdmTMB/sdmTMBextra

* Switch to using the new fmesher package for all mesh/SPDE calculations. INLA
  is no longer a dependency.
  
* Switch to `diagonal.penalty = FALSE` in mgcv::smoothCon(). 
  This changes the scale of the linear component of the smoother, but
  should result in the same model.
  https://github.com/glmmTMB/glmmTMB/issues/928#issuecomment-1642862066
  
* Implement cross validation for delta models #239

* Remove ELPD from cross validation output. Use sum_loglik instead. #235

* Turn on Newton optimization by default. #182

* print() now checks sanity() and issues a warning if there may be issues. #176

* Poisson-link delta models and censored likelihood distributions have been made
  considerably more robust. #186
  
* Standard errors are now available on SD parameters etc. in tidy() #240

* Fix bug in print()/tidy() for delta-model positive model component sigma_E.
  A recently introduce bug was causing sigma_E for the 2nd model to be reported
  as the 1st model component sigma_E.
  
* Add new anisotropy plotting function.

* Add anisotropic range printing. #149 by @jdunic

# sdmTMB 0.3.0

* Create the sdmTMBextra package to remove rstan/tmbstan helpers, which
  were causing memory sanitizer errors on CRAN.
  https://github.com/sdmTMB/sdmTMBextra
  
* The following functions are affected:

  - `predict.sdmTMB()` now takes `mcmc_samples`, which is output from
    `sdmTMBextra::extract_mcmc()`.
  - `simulate.sdmTMB()` now takes `mcmc_samples`, which is output from
    `sdmTMBextra::extract_mcmc()`.
  - `residuals.sdmTMB()` now takes `mcmc_samples`, which is output
    `sdmTMBextra::predict_mle_mcmc()`. This only affects 
    `residuals(..., type = "mle-mcmc")`.

* Move `dharma_residuals()` to 
  [sdmTMBextra](https://github.com/sdmTMB/sdmTMBextra) to reduce heavy
  dependencies.
  
* See examples in the Bayesian and residuals vignettes or in the help files for
  those functions within sdmTMBextra.

# sdmTMB 0.2.2

* Various fixes to pass CRAN checks. #158

* Fix memory issue highlighted by Additional issues CRAN checks. #158

* 'offset' argument can now be a character value indicating a column name. This
  is the preferred way of using an offset with parallel cross validation. #165

* Fix parallel cross validation when using an offset vector. #165

* Add leave-future-out cross validation functionality. #156

* Example data `qcs_grid` is no longer replicated by year to save package
  space. #158

* Add message with `tidy(fit, "ran_pars")` about why SEs are NA.

* Add anisotropy to `print()` #157

* Fix `predict(..., type = "response", se_fit = TRUE)`, which involves issuing
  a warning and sticking to link space. #140
  
# sdmTMB 0.2.1

* Fixes for resubmission to CRAN.

# sdmTMB 0.2.0

* Initial submission to CRAN.

# sdmTMB 0.1.4

* Relax range parameter `sanity()` check from 1x to 1.5x the greatest
  distance in the data.

* Add Pearson residuals for several families. `residuals(fit, type = "pearson")`
  Useful for checking for overdispersion with N > 1 binomial or Poisson
  families, among other uses. See the `overdisp_fun()` function at:
  https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor

* Fix bug when using `residuals()` or `simulate()` with binomial families
  specified via `cbind()` or `weights = N`. The binomial sample size wasn't
  being passed through typically resulting in Inf/-Inf.

* Add mixture families: `gamma_mix()`, `lognormal_mix()` and associated
  delta/hurdle families: `delta_gamma_mix()`, `delta_lognormal_mix()`. These
  families feature a mixture of two distributions with different means but
  shared variance parameters.
  
* Add `delta_beta()` family.

# sdmTMB 0.1.3

* Tweak `sanity()` checking of standard error size.

* Export previously experimental `plot_anisotropy()` function. The old function
is now `plot_anisotropy2()`.

* Allow passing offset data through `predict.sdmTMB()` via `offset` argument.

# sdmTMB 0.1.2

* Switch `effects = 'ran_vals'` for random intercept values from `tidy.sdmTMB()`
  to match the broom.mixed package.

* Make `tidy.sdmTMB()` return a tibble if the tibble package is installed. Note 
  this could affect old code since `drop = FALSE` is the default for tibbles
  but `drop = TRUE` is the default for data frames (i.e., tibbles always return
  a data frame when subsetted).

* Fix longstanding issue with predicting on newdata with mgcv's `t2()`. 
  Previously this was disabled because of issues. It now works as expected.
  
* Add `knots` argument in `sdmTMB()`, which is passed to mgcv. A common use
  would be to specify end points in a cyclical spline 
  (e.g., `s(x, bs = 'cc', k = 4), knots = list(x = c(1, 3, 5, 7))`) when the
  data don't extend fully to the boundaries that should match up.

# sdmTMB 0.1.1

* Preparing for release on CRAN.

* Add time-varying AR1 option (originally was always a random walk). See 
  `time_varying_type` argument in `?sdmTMB`.

* Allow prediction on `newdata` with missing time elements. #130

* Add check for `offset()` (which *does not* work in sdmTMB, use the `offset`
  argument instead). #131
  
* Add check for random slopes (sdmTMB currently only does random intercepts,
  although slopes can vary spatially). #131

# sdmTMB 0.1.0

* ADREPORT several parameters in natural space.
  <https://github.com/sdmTMB/sdmTMB/discussions/113>

* Improve robustness of model `print()` to more esoteric mgcv smoothers.

* Let `sims_var` work with multiple spatially varying slopes (`zeta_s`); return
  output in named list by coefficients. #107
  
* Add `threshold_coefs` to `sdmTMB_simulate()`.

* Don't make a fake mesh for non-spatial model (faster).

# sdmTMB 0.0.26.9001

* Add vignettes on visreg, ggeffects, and delta families (thanks J. Indivero!)
  #83 #87 #89 Forecasting and presence-only vignettes to be merged in soon.

* Add support for emmeans package. See `?emmeans.sdmTMB` for examples.

* Add support for effects package. The `ggeffects::ggeffect()` function
  can be used to make fast marginal effects plots. `ggeffects::ggpredict()`
  works with a custom fork of ggeffects. A pull request will be made shortly.
  #101 

* Add `vcov()`, `fixef()`, `df.residual`(), `formula()`, `terms()`, and
  `model.frame()` methods.

* Add support for `"cloglog"` link. Code adapted from glmmTMB for robust
  likelihood implementation.

* For delta models, by default share the anisotropy parameters as in VAST. 
  Separate anisotropy (old behavior) can be estimated with
  `control = sdmTMBcontrol(map = list(ln_H_input = factor(c(1, 2, 3, 4))))`
  
* Add experimental `do_index`, `predict_args`, and `index_args` in `sdmTMB()`.
  These can be used to perform prediction and index calculation at the same
  time as fitting. For very large datasets or meshes this can save time
  compared to fitting, predicting, and index calculation in 3 separate steps
  since the TMB AD object doesn't have to be rebuilt. This will somewhat slow
  down the initial fitting.
  
* Remove `max_gradient` and `bad_eig` from `get_index()` output.

* Use unique locations on prediction for huge speedups on large `newdata`
  gridded data.

* Fix bug where in rare cases `get_index()` would return gibberish small values.

* Add `bayesian` argument, which when `TRUE` adds Jacobian adjustments for
  non-linear transformed parameters. This should be `TRUE` if the model 
  will be passed to tmbstan, but `FALSE` otherwise. #95
  
* Add experimental and not-yet-exported `sdmTMB:::plot_anisotropy2()`.

* Add many anisotropy, delta model, and index calculation unit tests.

# sdmTMB 0.0.24.9001

* Enable random walk random field TMB simulation in `sdmTMB_simulate()`.

* Add check for irregular time with AR1 or random walk processes.

* Fix bugs introduced by delta model code (offsets with `extra_time` and
  threshold model prediction).
  
* Fix bug in `sanity()` message with small random field SDs.

# sdmTMB 0.0.24.9000

* Add support for 'delta' (or 'hurdle') models. See examples and documentation
  in `?sdmTMB`. This has resulted in a substantial restructuring of the
  internal model code. By default both model components (e.g., binomial & Gamma)
  share the same formula, spatial, and spatiotemporal structure, but these
  can be separated by supplying argument values in lists where the first
  element corresponds to the first model and the second element corresponds to
  the second model (with some limitations as described in `?sdmTMB`
  documentation 'Details').

* Add support for multiple spatially varying coefficients (used to be limited to
  a single variable).

* Add compatibility with the 'visreg' package for visualizing conditional
  effects of parameters. See `?visreg_delta` for examples.

* Add MCMC residual type to `residuals.sdmTMB()`. These are a 'better' residuals
  but slower to calculate. See documentation 'Details' in `?residuals.sdmTMB`.

* Make `offset` an argument in `sdmTMB()`. Using the reserved word `offset` in
  the formula is now deprecated.

* Add `sanity()` function to perform some basic sanity checks on model fits.

* Make an `sdmTMB()` model object compatible with `update()` method.

* Remove several deprecated arguments.

* Overhaul examples in `?sdmTMB`.

* Use faster "low-rank sparse hessian bias-correction" TMB bias correction.

* Add parallel processing support. See `parallel` argument in `sdmTMBcontrol`.
  By default, grabs value of `sdmTMB.cores` option. E.g.
  `options(sdmTMB.cores = 4)`. Only currently enabled on Mac/Linux.
  Using too many cores can be much slower than 1 core.
  
* Use 'cli' package `cli_abort()`/`cli_warn()`/`cli_inform()` over
  `stop()`/`warning()`/`message()`.

* Add many unit tests.

# sdmTMB 0.0.23.9000

* A package version number that was used for internal testing in the 'delta'
  branch by several people.

# sdmTMB 0.0.22.9001

* Switch to TMBad library for ~3-fold speedup(!)

# sdmTMB 0.0.22.9000

* Fix bug in predictions with `poly(..., raw = FALSE)` on newdata. #77

# sdmTMB 0.0.21.9009

* Add experimental `sdmTMB_stacking()` for ensemble model stacking weights.

* Add fake mesh if random fields are all off. #59

* Make `predict(..., newdata = NULL)` also use `last.par.best` instead of
  `last.par` to match `newdata = df`.

* Fix bug in MVN fixed-effect prior indexing

* `sims` and `n_sims` arguments have been deprecated and standardized
  to `nsim` to match the `simulate()` S3 method.

* Bias correction on `get_index()` and `get_cog()` is now selective and 
  is just applied to the necessary derived parameters.

* INLA projection matrix 'A' is now shared across spatial and spatiotemporal
  fields.

* Add `add_utm_columns()` to ease adding UTM columns.

# sdmTMB 0.0.20.9001

* Add `dharma_residuals()`.

* Fix bug in `simulate.sdmTMB()` and `residuals.sdmTMB()` for binomial family.

# sdmTMB 0.0.20.9000

*  Smoothers now appear in `print()` output. The format should roughly match brms.
   The main-effect component (e.g., `sdepth` for `s(depth)`) represents the
   linear component and the random effect (e.g., `sds(depth)`) component in
   the output corresponds to the standard deviation of the penalized weights.

*  Add `censored_poisson(link = 'log')` family; implemented by @joenomiddlename

* `fields` in `sdmTMB()` is now deprecated and replaced by `spatiotemporal`.

* `include_spatial` in `sdmTMB()` is now deprecated and replaced by `spatial`.

* `spatial_only` in `sdmTMB()` is now deprecated and replaced by `spatiotemporal`. 
   E.g. `spatial_only = TRUE` is now `spatiotemporal = 'off'` or leaving 
   `time = NULL`.
   
* `spde` in `sdmTMB()` is now deprecated and replaced by `mesh`.

* `sdmTMB_simulate()` is new and will likely replace `sdmTMB_sim()` eventually. 
  `sdmTMB_simulate()` is set up to take a formula and a data frame and is easier
  to use if you want different spatial observations (and covariates) for each
  time slice. It can also take a fitted model and modify parts of it to simulate.
  Finally, this function uses TMB for simulation and so is much faster and
  more flexible in what it can simulate (e.g., anisotropy) than the previous version.
  
* `spatial_trend` is now `spatial_varying` and accepts a one-sided formula
  *with a single predictor* of any coefficient that should varying in space as a
  random field. Note that you may want to include a fixed effect for the same
  variable to improve interpretability. If the (scaled) time column is used, it will
  represent a local-time-trend model as before.
  
* The Tweedie power (p) parameter is now in `print()` and `tidy()` output.

* `thetaf` is now `tweedie_p` in `sdmTMB_sim()`.

# sdmTMB 0.0.19.9003

* Fix bug affecting prediction with `se_fit = TRUE` for breakpoint models.

# sdmTMB 0.0.19.9002

* Simulation from the parameter covariance matrix works if random effects
  are turned off. #57

# sdmTMB 0.0.19.9000

* Smoothers `s()` are now *penalized* smoothers: they determine the 
  degree of wiggliness (as in mgcv) and it is no longer necessary to
  choose an appropriate `k` value a priori. Models fit with previous
  versions of sdmTMB with  `s(x, k = ...)` will not match models
  specified the same way in version >= 0.0.19 since the basis functions
  are now penalized. All the various `mgcv::s()` options should be supported
  but `t2()` (and `ti()` and `te()`) are not supported.

# sdmTMB 0.0.18.9001

* Add ELPD (expected log predictive density) to `sdmTMB_cv()`
  <https://arxiv.org/abs/1507.04544>
  
* Fix bug evaluating `...` when `sdmTMB_cv()` was called within a function. #54

# sdmTMB 0.0.18.9000

* Fix minor error in PC Matern prior

# sdmTMB 0.0.17.9000

* Add random walk option: `fields = "RW"`.

* Depreciate `ar1_fields` argument. See new `fields` argument in `sdmTMB().

* Many packages moved from 'Imports' to 'Suggests'

# sdmTMB 0.0.16.9000

* Lower default `nlminb()` `eval.max` and `iter.max` to 1000 and 2000.

* Added `profile` option in `sdmTMBcontrol()`. This can dramatically
  improve model fitting speed with many fixed effects. Note the
  result is likely to be slightly different with `TRUE` vs. `FALSE`.

* Added simulation from the MVN precision matrix to `predict.sdmTMB()`. 
  See the `sims` argument.
  
* Added `gather_sims()` and `spread_sims()` to extract parameter
  simulations from the joint precision matrix in a format that
  matches the tidybayes package.

* Added `get_index_sims()` for a population index calculated from
  the MVN simulation draws.
  
* Added `extract_mcmc()` to extract MCMC samples if the model is
  passed to tmbstan.
  
* Added the ability to predict from a model fitted with tmbstan.
  See the `tmbstan_model` argument in `predict.sdmTMB()`.

* Allowed for separate random field Matern range parameters for 
  spatial and spatiotemporal fields. E.g. `sdmTMB(shared_range = FALSE)`

* Bounded the AR1 rho parameter between -0.999 and 0.999 to improve 
  convergence; was -1 to 1. Please post an issue if this creates
  problems for your model.

* Added `map`, `start`, `lower`, and `upper` options to control
  model fitting. See `sdmTMBcontrol()`.

* Added priors for all parameters. See `?sdmTMB::priors` and the
  `priors` argument in `sdmTMB()`. PC priors are available for
  the random fields. See `?pc_matern` and the details there.
  
* Moved many less-common arguments from `sdmTMB()` to `sdmTMBcontrol()`.

* Fix bug in `sdmTMB_cv()` where fitting and testing data splits
  were reversed. I.e., the small chunk was fit; the big chunk was tested.

# sdmTMB 0.0.15.9000

* Added experimental penalized complexity (PC) prior as used in INLA.
  See arguments `matern_prior_O` and `matern_prior_E`.

* Added back `normalize` argument to `sdmTMB()` and default to `FALSE`.
  Setting to `TRUE` can dramatically speed up some model fits
  (~4 times for some test models).

# sdmTMB 0.0.14.9003

* Added vignette on making pretty maps of the output

# sdmTMB 0.0.14.9001

* Added some protections for possible user errors:
  * AR1 with a spatial-only model
  * Missing factor levels in time
  * Coordinate systems that are too big

# sdmTMB 0.0.14.9000

* Add `re_form_iid` to `predict.sdmTMB()`.

* Add `map_rf` option to `sdmTMB()`. This lets you map (fix at
  their starting values of zero) all random fields to produce a
  classic GLM/GLMM.

# sdmTMB 0.0.13.9000

* Add IID random intercepts interface. E.g. `... + (1 | g)` #34

# sdmTMB 0.0.12.9000

* Add `epsilon_predictor` argument in `sdmTMB()` to allow a model of the
  spatiotemporal variance through time.

# sdmTMB 0.0.11.9000

* Add `penalties` argument to allow for regularization.

# sdmTMB 0.0.10.9001

* Fix Student-t degrees of freedom in the randomized quantile residuals

# sdmTMB 0.0.10.9000

* Fixed parameter initialization for inverse links #35

* Switched Gamma 'phi' parameter to representing shape instead of CV to
  match glm(), glmmTMB(), etc.

# sdmTMB 0.0.9.9000

* Switched the density/abundance index calculation to use the link function as
  opposed to a hardcoded log() so that the `get_generic()` function can be used
  to grab things like standardized average values of the response across a grid.
  What used to be `log_total` in the raw TMB output is now `link_total` but most
  users you shouldn't notice any difference.

# sdmTMB 0.0.8.9000

* Overhauled the simulation function. The function is now called `sdmTMB_sim()`
  and uses INLA functions instead of RandomFields functions for simulating
  the random fields.

* The simulation function can now accommodate all families and links and takes
  an INLA mesh as input.

# sdmTMB 0.0.7.9001

* Allow specifying degrees of freedom in the Student-t family #29

# sdmTMB 0.0.7.9000

* Added a `tidy()` method (from broom and broom.mixed) to return a data frame
  of parameter estimates. The function can extract the fixed effects or the
  random effect parameters (variances, AR1 correlation, spatial range).

* Added an argument `extra_time` to `sdmTMB()`. This introduces additional time
  slices that you can then predict on if you want to interpolate or forecast.
  Internally, it uses Eric Ward's 'weights hack'. This is also useful if you
  have data unevenly spaced in time and you want the gaps evenly spaced for a
  random walk or AR1 process (add any missing years to `extra_time`).

* `make_spde()` is now replaced with `make_mesh()` and `make_spde()` has been
  soft deprecated. `make_mesh()` carries through the x and y column names to
  the predict function and is more in line with the tidyverse style of taking a
  data frame first.

* `make_mesh()` can accept `cutoff` as an argument (as in INLA), which is
  likely a better default way to specify the mesh since it scales across
  regions better and is line with the literature on INLA.

* `make_mesh()` can use a binary search algorithm to find a cutoff that best
  matches a desired number of knots (thanks to Kelli Johnson for the idea).

* Barrier meshes are now possible. See `add_barrier_mesh()` for an example.

* There is a pkgdown website now that gets auto generated with GitHub actions.

* There is the start of a model description vignette.
  It is very much a work in progress.

# sdmTMB 0.0.6.9009

* Fixed bug in dlnorm

# sdmTMB 0.0.6.9005

* Fixed bug in predictions with standard errors where one(?)
  parameter (a breakpoint parameter) would be passed in at its initial
  instead of MLE value.

# sdmTMB 0.0.6.9004

* Fixed bug with predictions on new data in models with break points

* Overhauled cross validation function. The function now:
    * uses Eric's weights hack so it can also be used for forecasting
    * initializes subsequent folds at the MLE of the first fold for
      considerable speed increases
    * works in parallel if a future plan initialized; see examples

* Added threshold parameters to the print method

* Added forecasting example with the weights hack

* Fixed bug in linear break point models

#  sdmTMB 0.0.6.9002

* Fixed GAM predictions with all 0s in new data.

* Add linear and logistic threshold models. #17

# sdmTMB 0.0.5.9000

* Added parsing of mgcv formulas for splines. #16

* Added ability to predict with standard errors at the population level. This
  helps with making marginal-effect plots. #15

* Added optimization options to aid convergence. Also added
  `run_extra_optimization()` to run these on already fit models. Default is for
  no extra optimization.

* Added binomial likelihood to cross validation. Git hash `ee3f3ba`.

* Started keeping track of news in `NEWS.md`.
