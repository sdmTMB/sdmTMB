# Turn sdmTMB model output into a tidy data frame

Turn sdmTMB model output into a tidy data frame

## Usage

``` r
# S3 method for class 'sdmTMB'
tidy(
  x,
  effects = c("fixed", "ran_pars", "ran_vals", "ran_vcov"),
  model = 1,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  silent = FALSE,
  ...
)

# S3 method for class 'sdmTMB_cv'
tidy(x, ...)
```

## Arguments

- x:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- effects:

  A character value. One of `"fixed"` ('fixed' or main-effect
  parameters), `"ran_pars"` (standard deviations, spatial range, and
  other random effect and dispersion-related terms), `"ran_vals"`
  (individual random intercepts or slopes, if included; behaves like
  `ranef()`), or `"ran_vcov"` (list of variance covariance matrices for
  the random effects, by model and group).

- model:

  Which model to tidy if a delta model (1 or 2). The `model` will be
  ignored when effects is `"ran_vals"` (all returned in a single
  dataframe)

- conf.int:

  Include a confidence interval?

- conf.level:

  Confidence level for CI.

- exponentiate:

  Whether to exponentiate the fixed-effect coefficient estimates and
  confidence intervals.

- silent:

  Omit any messages?

- ...:

  Extra arguments (not used).

## Value

A data frame

## Details

Follows the conventions of the broom and broom.mixed packages.

Currently, `effects = "ran_pars"` also includes dispersion-related terms
(e.g., `phi`), which are not actually associated with random effects.

Standard errors for spatial variance terms fit in log space (e.g.,
variance terms, range, or parameters associated with the observation
error) are omitted to avoid confusion. Confidence intervals are still
available.

## Examples

``` r
fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE),
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie()
)
tidy(fit)
#> # A tibble: 3 × 5
#>   term                               estimate std.error conf.low conf.high
#>   <chr>                                 <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)                            3.65     0.281     3.10     4.20 
#> 2 poly(depth_scaled, 2, raw = TRUE)1    -1.54     0.186    -1.90    -1.17 
#> 3 poly(depth_scaled, 2, raw = TRUE)2    -1.11     0.101    -1.31    -0.913
tidy(fit, conf.int = TRUE)
#> # A tibble: 3 × 5
#>   term                               estimate std.error conf.low conf.high
#>   <chr>                                 <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)                            3.65     0.281     3.10     4.20 
#> 2 poly(depth_scaled, 2, raw = TRUE)1    -1.54     0.186    -1.90    -1.17 
#> 3 poly(depth_scaled, 2, raw = TRUE)2    -1.11     0.101    -1.31    -0.913
tidy(fit, "ran_pars", conf.int = TRUE)
#> # A tibble: 4 × 5
#>   term      estimate std.error conf.low conf.high
#>   <chr>        <dbl>     <dbl>    <dbl>     <dbl>
#> 1 range        19.1    14.0       4.58      80.0 
#> 2 phi          14.0     0.677    12.8       15.4 
#> 3 sigma_O       2.14    0.941     0.906      5.07
#> 4 tweedie_p     1.58    0.0153    1.55       1.61

pcod_2011$fyear <- as.factor(pcod_2011$year)
fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE) + (1 | fyear),
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie()
)
tidy(fit, "ran_vals")
#> # A tibble: 4 × 7
#>   group_name term        level_ids estimate std.error conf.low conf.high
#>   <chr>      <chr>       <chr>        <dbl>     <dbl>    <dbl>     <dbl>
#> 1 fyear      (Intercept) 2011        0.0163     0.187   -0.351    0.384 
#> 2 fyear      (Intercept) 2013        0.177      0.188   -0.192    0.545 
#> 3 fyear      (Intercept) 2015        0.232      0.189   -0.139    0.603 
#> 4 fyear      (Intercept) 2017       -0.431      0.205   -0.833   -0.0287
```
