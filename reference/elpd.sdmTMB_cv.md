# Expected log predictive density from cross validation for loo

Converts the pointwise out-of-sample log predictive densities
(`cv_loglik`) from
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
into a loo `elpd_generic` object. This lets
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
compare models, including a standard error of the difference in expected
log predictive density.

## Usage

``` r
# S3 method for class 'sdmTMB_cv'
elpd(x, ...)
```

## Arguments

- x:

  Output from
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md).

- ...:

  Not used.

## Value

An object of class `elpd_generic` and `loo`. See
[`loo::elpd()`](https://mc-stan.org/loo/reference/elpd.html).

## Details

Models compared with
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
must be fit to the same data with the same folds (e.g., via `fold_ids`)
and should use the same `predictive` argument. loo only checks that the
number of observations matches. The standard errors reflect variation
across observations, not the Monte Carlo error from `nsim`.

The printed object refers to a "1 by N log-likelihood matrix" because
the pointwise values are already integrated over draws within
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md).

## See also

[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md),
[`compare_deviance()`](https://sdmTMB.github.io/sdmTMB/reference/compare_deviance.md)
for in-sample deviance explained.

## Examples

``` r
# \donttest{
if (requireNamespace("loo", quietly = TRUE)) {
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  folds <- rep(1:4, length.out = nrow(pcod_2011))
  m1 <- sdmTMB_cv(density ~ 1, data = pcod_2011, mesh = mesh,
    family = tweedie(), fold_ids = folds)
  m2 <- sdmTMB_cv(density ~ depth_scaled + depth_scaled2,
    data = pcod_2011, mesh = mesh, family = tweedie(), fold_ids = folds)
  loo::loo_compare(list(null = loo::elpd(m1), depth = loo::elpd(m2)))
}
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
#>  model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  depth       0.0     0.0      NA                    
#>   null     -73.0    15.1    1.00                    
# }
```
