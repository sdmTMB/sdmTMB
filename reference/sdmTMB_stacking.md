# Perform stacking with log scores on `sdmTMB_cv()` output

**\[experimental\]**

This approach is described in Yao et al. (2018)
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091) . The general
method minimizes (or maximizes) some quantity across models. For simple
models with normal error, this may be the root mean squared error
(RMSE), but other approaches include the log score. We adopt the latter
here, where log scores are used to generate the stacking of predictive
distributions

## Usage

``` r
sdmTMB_stacking(model_list, include_folds = NULL)
```

## Arguments

- model_list:

  A list of models fit with
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
  to generate estimates of predictive densities. You will want to set
  the seed to the same value before fitting each model or manually
  construct the fold IDs so that they are the same across models.

- include_folds:

  An optional numeric vector specifying which folds to include in the
  calculations. For example, if 5 folds are used for k-fold cross
  validation, and the first 4 are needed to generate these weights,
  `include_folds = 1:4`.

## Value

A vector of model weights.

## References

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. 2018. Using Stacking
to Average Bayesian Predictive Distributions (with Discussion). Bayesian
Analysis 13(3): 917–1007. International Society for Bayesian Analysis.
[doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

## Examples

``` r
# \donttest{
# Set parallel processing if desired. See 'Details' in ?sdmTMB_cv

# Depth as quadratic:
set.seed(1)
m_cv_1 <- sdmTMB_cv(
  density ~ 0 + depth_scaled + depth_scaled2,
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
# Depth as linear:
set.seed(1)
m_cv_2 <- sdmTMB_cv(
  density ~ 0 + depth_scaled,
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# Only an intercept:
set.seed(1)
m_cv_3 <- sdmTMB_cv(
  density ~ 1,
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie(link = "log"), k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

models <- list(m_cv_1, m_cv_2, m_cv_3)
weights <- sdmTMB_stacking(models)
weights
#> [1] 0.87857001 0.01750509 0.10392489
# }
```
