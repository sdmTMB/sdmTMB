# Cross validation with sdmTMB models

Performs k-fold or leave-future-out cross validation with sdmTMB models.
Returns the sum of log likelihoods of held-out data (log predictive
density), which can be used to compare models—higher values indicate
better out-of-sample prediction. By default, creates folds randomly and
stratified by time (set a seed for reproducibility), but folds can be
manually assigned via `fold_ids`. See Ward and Anderson (2025) in the
References and the [cross-validation
vignette](https://sdmTMB.github.io/sdmTMB/articles/cross-validation.html).

## Usage

``` r
sdmTMB_cv(
  formula,
  data,
  mesh_args,
  mesh = NULL,
  time = NULL,
  k_folds = 8,
  fold_ids = NULL,
  lfo = FALSE,
  lfo_forecast = 1,
  lfo_validations = 5,
  parallel = TRUE,
  use_initial_fit = FALSE,
  save_models = TRUE,
  future_globals = NULL,
  ...
)
```

## Arguments

- formula:

  Model formula.

- data:

  A data frame.

- mesh_args:

  Arguments for
  [`make_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md).
  If supplied, the mesh will be reconstructed for each fold.

- mesh:

  Output from
  [`make_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md).
  If supplied, the same mesh will be used for all folds. This is faster
  and usually what you want.

- time:

  The name of the time column. Leave as `NULL` if this is only spatial
  data.

- k_folds:

  Number of folds.

- fold_ids:

  Optional vector containing user fold IDs. Can also be a single string,
  e.g. `"fold_id"` representing the name of the variable in `data`.
  Ignored if `lfo` is TRUE

- lfo:

  Logical. Use leave-future-out (LFO) cross validation? If `TRUE`, data
  from earlier time steps are used to predict future time steps. The
  `time` argument must be specified. See Details section below.

- lfo_forecast:

  If `lfo = TRUE`, number of time steps ahead to forecast. For example,
  `lfo_forecast = 1` means fitting to time steps 1 to T and validating
  on T + 1. See Details section below.

- lfo_validations:

  If `lfo = TRUE`, number of times to step through the LFO process
  (i.e., number of validation folds). Defaults to 5. See Details section
  below.

- parallel:

  If `TRUE` and a
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  is supplied, will be run in parallel.

- use_initial_fit:

  Fit the first fold and use those parameter values as starting values
  for subsequent folds? Can be faster with many folds.

- save_models:

  Logical. If `TRUE` (default), the fitted model object for each fold is
  stored in the output. If `FALSE`, models are not saved, which can
  substantially reduce memory usage for large datasets or many folds.
  When `FALSE`, functions that require access to the fitted models
  (e.g., [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`cv_to_waywiser()`](https://sdmTMB.github.io/sdmTMB/reference/cv_to_waywiser.md))
  will not work.

- future_globals:

  A character vector of global variables used within arguments if an
  error is returned that future.apply can't find an object. This vector
  is appended to `TRUE` and passed to the argument `future.globals` in
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).
  Useful if global objects are used to specify arguments like priors,
  families, etc.

- ...:

  All other arguments required to run the
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
  model. The `weights` argument is supported and will be combined with
  the internal fold-assignment mechanism (held-out data are assigned
  weight 0).

## Value

A list:

- `data`: Original data plus columns for fold ID (`cv_fold`), CV
  predicted value (`cv_predicted`), CV log likelihood (`cv_loglik`), and
  CV deviance residuals (`cv_deviance_resid`).

- `models`: A list of fitted models, one per fold. `NULL` if
  `save_models = FALSE`.

- `fold_loglik`: Sum of log likelihoods of held-out data per fold (log
  predictive density per fold). More positive values indicate better
  out-of-sample prediction.

- `sum_loglik`: Sum of `fold_loglik` across all folds (total log
  predictive density). Use this to compare models—more positive values
  are better.

- `pdHess`: Logical vector: was the Hessian positive definite for each
  fold?

- `converged`: Logical: did all folds converge (all `pdHess` `TRUE`)?

- `max_gradients`: Maximum absolute gradient for each fold.

## Details

**Parallel processing**

Parallel processing can be used by setting a
[`future::plan()`](https://future.futureverse.org/reference/plan.html).

For example:

    library(future)
    plan(multisession)
    # now use sdmTMB_cv() ...

**Leave-future-out cross validation (LFOCV)**

An example of LFOCV with 9 time steps, `lfo_forecast = 1`, and
`lfo_validations = 2`:

- Fit data to time steps 1 to 7, predict and validate step 8.

- Fit data to time steps 1 to 8, predict and validate step 9.

An example of LFOCV with 9 time steps, `lfo_forecast = 2`, and
`lfo_validations = 3`:

- Fit data to time steps 1 to 5, predict and validate step 7.

- Fit data to time steps 1 to 6, predict and validate step 8.

- Fit data to time steps 1 to 7, predict and validate step 9.

Note these are time steps as they are presented in order in the data.
For example, in the `pcod` data example below steps between data points
are not always one year but an `lfo_forecast = 2` forecasts 2 time steps
as presented not two years.

See example below.

## References

Ward, E.J., and S.C. Anderson. 2025. Approximating spatial processes
with too many knots degrades the quality of probabilistic predictions.
bioRxiv 2025.11.14.688354.
[doi:10.1101/2025.11.14.688354](https://doi.org/10.1101/2025.11.14.688354)
.

## Examples

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)

# Set parallel processing first if desired with the future package.
# See the Details section above.

m_cv <- sdmTMB_cv(
  density ~ 0 + depth_scaled + depth_scaled2,
  data = pcod, mesh = mesh, spatial = "off",
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

m_cv$fold_loglik
#> [1] -4292.305 -4321.288
m_cv$sum_loglik
#> [1] -8613.593

head(m_cv$data)
#> # A tibble: 6 × 16
#>    year     X     Y depth density present   lat   lon depth_mean depth_sd
#>   <int> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl> <dbl>      <dbl>    <dbl>
#> 1  2003  446. 5793.   201   113.        1  52.3 -130.       5.16    0.445
#> 2  2003  446. 5800.   212    41.7       1  52.3 -130.       5.16    0.445
#> 3  2003  449. 5802.   220     0         0  52.4 -130.       5.16    0.445
#> 4  2003  437. 5802.   197    15.7       1  52.4 -130.       5.16    0.445
#> 5  2003  421. 5771.   256     0         0  52.1 -130.       5.16    0.445
#> 6  2003  418. 5772.   293     0         0  52.1 -130.       5.16    0.445
#> # ℹ 6 more variables: depth_scaled <dbl>, depth_scaled2 <dbl>, cv_fold <int>,
#> #   cv_predicted <dbl>, cv_deviance_resid <dbl>, cv_loglik <dbl>
m_cv$models[[1]]
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ 0 + depth_scaled + depth_scaled2
#> Mesh: sdmTMBmesh (isotropic covariance)
#> Data: tbl_df
#> Family: tweedie(link = 'log')
#>  
#> Conditional model:
#>               coef.est coef.se
#> depth_scaled     -1.83    0.15
#> depth_scaled2     2.00    0.21
#> 
#> Dispersion parameter: 76.49
#> Tweedie p: 1.90
#> ML criterion at convergence: 4280.128
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.
m_cv$max_gradients
#> [1] 2.272052e-11 2.044996e-11


# \donttest{
# Create mesh each fold:
m_cv2 <- sdmTMB_cv(
  density ~ 0 + depth_scaled + depth_scaled2,
  data = pcod, mesh_args = list(xy_cols = c("X", "Y"), cutoff = 20),
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# Use fold_ids:
m_cv3 <- sdmTMB_cv(
  density ~ 0 + depth_scaled + depth_scaled2,
  data = pcod, mesh = mesh,
  family = tweedie(link = "log"),
  fold_ids = rep(seq(1, 3), nrow(pcod))[seq(1, nrow(pcod))]
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# LFOCV:
m_lfocv <- sdmTMB_cv(
  present ~ s(year, k = 4),
  data = pcod,
  lfo = TRUE,
  lfo_forecast = 2,
  lfo_validations = 3,
  family = binomial(),
  mesh = mesh,
  spatial = "off", # fast example
  spatiotemporal = "off", # fast example
  time = "year" # must be specified
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# See how the LFOCV folds were assigned:
fold_table <- table(m_lfocv$data$cv_fold, m_lfocv$data$year)
fold_table
#>    
#>     2013 2015 2017
#>   3  240    0    0
#>   4    0  238    0
#>   5    0    0  240
# }
```
