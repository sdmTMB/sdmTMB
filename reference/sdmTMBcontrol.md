# Optimization control options

[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) and
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) control
options.

## Usage

``` r
sdmTMBcontrol(
  eval.max = 2000L,
  iter.max = 1000L,
  normalize = FALSE,
  nlminb_loops = 1L,
  newton_loops = 1L,
  getsd = TRUE,
  quadratic_roots = FALSE,
  start = NULL,
  map = NULL,
  lower = NULL,
  upper = NULL,
  censored_upper = NULL,
  multiphase = TRUE,
  profile = FALSE,
  get_joint_precision = TRUE,
  parallel = getOption("sdmTMB.cores", 1L),
  suppress_nlminb_warnings = TRUE,
  collapse_spatial_variance = FALSE,
  collapse_threshold = 0.01,
  ...
)
```

## Arguments

- eval.max:

  Maximum number of evaluations of the objective function allowed.

- iter.max:

  Maximum number of iterations allowed.

- normalize:

  Logical: use
  [`TMB::normalize()`](https://rdrr.io/pkg/TMB/man/normalize.html) to
  normalize the process likelihood using the Laplace approximation? Can
  result in a substantial speed boost in some cases. This used to
  default to `FALSE` prior to May 2021. Currently not working for models
  fit with REML or random intercepts.

- nlminb_loops:

  How many times to run
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) optimization.
  Sometimes restarting the optimizer at the previous best values aids
  convergence. If the maximum gradient is still too large, try
  increasing this to `2`.

- newton_loops:

  How many Newton optimization steps to try after running
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). This
  sometimes aids convergence by further reducing the log-likelihood
  gradient with respect to the fixed effects. This calculates the
  Hessian at the current MLE with
  [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) using a
  finite-difference approach and uses this to update the fixed effect
  estimates.

- getsd:

  Logical indicating whether to call
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- quadratic_roots:

  Experimental feature for internal use right now; may be moved to a
  branch. Logical: should quadratic roots be calculated? Note: on the
  sdmTMB side, the first two coefficients are used to generate the
  quadratic parameters. This means that if you want to generate a
  quadratic profile for depth, and depth and depth^2 are part of your
  formula, you need to make sure these are listed first and that an
  intercept isn't included. For example,
  `formula = cpue ~ 0 + depth + depth2 + as.factor(year)`.

- start:

  A named list specifying the starting values for parameters. You can
  see the necessary structure by fitting the model once and inspecting
  `your_model$tmb_obj$env$parList()`. Elements of `start` that are
  specified will replace the default starting values.

- map:

  A named list with factor `NA`s specifying parameter values that should
  be fixed at a constant value. See the documentation in
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html). This
  should usually be used with `start` to specify the fixed value.

- lower:

  An optional named list of lower bounds within the optimization.
  Parameter vectors with the same name (e.g., `b_j` or `ln_kappa` in
  some cases) can be specified as a numeric vector. E.g.
  `lower = list(b_j = c(-5, -5))`. Note that
  [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) does not
  implement lower and upper bounds, so you must set `newton_loops = 0`
  if setting limits.

- upper:

  An optional named list of upper bounds within the optimization.

- censored_upper:

  An optional vector of upper bounds for `sdmTMBcontrol()`. Values of
  `NA` indicate an unbounded right-censored to distribution, values
  greater that the observation indicate and upper bound, and values
  equal to the observation indicate no censoring.

- multiphase:

  Logical: estimate the fixed and random effects in phases? Phases are
  usually faster and more stable.

- profile:

  Logical: should population-level/fixed effects be profiled out of the
  likelihood? These are then appended to the random effects vector
  without the Laplace approximation. See
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).
  *This can dramatically speed up model fit if there are many fixed
  effects but is experimental at this stage.*

- get_joint_precision:

  Logical. Passed to `getJointPrecision` in
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html). Must
  be `TRUE` to use simulation-based methods in
  [`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
  or `[get_index_sims()]`. If not needed, setting this `FALSE` will
  reduce object size.

- parallel:

  Argument currently ignored. For parallel processing with 3 cores, as
  an example, use `TMB::openmp(n = 3, DLL = "sdmTMB")`. But be careful,
  because it's not always faster with more cores and there is definitely
  an upper limit.

- suppress_nlminb_warnings:

  Suppress uninformative warnings from
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) arising when
  a function evaluation is `NA`, which are then replaced with `Inf` and
  avoided during estimation?

- collapse_spatial_variance:

  Logical: should spatial and/or spatiotemporal random fields be
  automatically dropped if their estimated standard deviation is
  effectively zero (i.e., below `collapse_threshold`)? This helps
  prevent overfitting and numerical instability when the data provide
  little evidence for spatial or spatiotemporal variation. I.e., when
  the variance parameter is estimated on or near the boundary of zero.
  When enabled, the model will be automatically refitted via
  [`update.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/update.sdmTMB.md)
  with the corresponding field(s) disabled. This adds a computational
  cost (a single model refit if collapsing occurs) but can yield a
  simpler, more stable model and more reliable inference. Default is
  `FALSE` for backwards compatibility.

- collapse_threshold:

  Numeric: the standard deviation threshold below which random fields
  are considered to be collapsing to zero. Only used when
  `collapse_spatial_variance = TRUE`. Values are on the standard
  deviation scale (i.e., square root of variance). Default is 0.01.

- ...:

  Anything else. See the 'Control parameters' section of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list of control arguments

## Details

Usually used within
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). For
example:

    sdmTMB(..., control = sdmTMBcontrol(newton_loops = 2))

## Examples

``` r
sdmTMBcontrol()
#> $eval.max
#> [1] 2000
#> 
#> $iter.max
#> [1] 1000
#> 
#> $normalize
#> [1] FALSE
#> 
#> $nlminb_loops
#> [1] 1
#> 
#> $newton_loops
#> [1] 1
#> 
#> $getsd
#> [1] TRUE
#> 
#> $profile
#> [1] FALSE
#> 
#> $quadratic_roots
#> [1] FALSE
#> 
#> $start
#> NULL
#> 
#> $map
#> NULL
#> 
#> $lower
#> NULL
#> 
#> $upper
#> NULL
#> 
#> $censored_upper
#> NULL
#> 
#> $multiphase
#> [1] TRUE
#> 
#> $parallel
#> [1] 1
#> 
#> $get_joint_precision
#> [1] TRUE
#> 
#> $collapse_spatial_variance
#> [1] FALSE
#> 
#> $collapse_threshold
#> [1] 0.01
#> 
```
