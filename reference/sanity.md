# Sanity check of an sdmTMB model

Sanity check of an sdmTMB model

## Usage

``` r
sanity(object, big_sd_log10 = 2, gradient_thresh = 0.001, silent = FALSE)
```

## Arguments

- object:

  Fitted model from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- big_sd_log10:

  Value to check size of standard errors against. A value of 2 would
  indicate that standard errors greater than `10^2` (i.e., 100) should
  be flagged.

- gradient_thresh:

  Gradient threshold to issue warning.

- silent:

  Logical: suppress messages? Useful to set to `TRUE` if running large
  numbers of models and just interested in returning sanity list
  objects.

## Value

An invisible named list of checks.

## Details

If `object` is `NA`, `NULL`, or of class `"try-error"`, `sanity()` will
return `FALSE`. This is to facilitate using `sanity()` on models with
[`try()`](https://rdrr.io/r/base/try.html) or
[`tryCatch()`](https://rdrr.io/r/base/conditions.html). See the examples
section.

## Examples

``` r
fit <- sdmTMB(
  present ~ s(depth),
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = binomial()
)
sanity(fit)
#> ✔ Non-linear minimizer suggests successful convergence
#> ✔ Hessian matrix is positive definite
#> ✔ No extreme or very small eigenvalues detected
#> ✔ No gradients with respect to fixed effects are >= 0.001
#> ✔ No fixed-effect standard errors are NA
#> ✔ No standard errors look unreasonably large
#> ✔ No sigma parameters are < 0.01
#> ✔ No sigma parameters are > 100
#> ✔ Range parameter doesn't look unreasonably large

s <- sanity(fit)
#> ✔ Non-linear minimizer suggests successful convergence
#> ✔ Hessian matrix is positive definite
#> ✔ No extreme or very small eigenvalues detected
#> ✔ No gradients with respect to fixed effects are >= 0.001
#> ✔ No fixed-effect standard errors are NA
#> ✔ No standard errors look unreasonably large
#> ✔ No sigma parameters are < 0.01
#> ✔ No sigma parameters are > 100
#> ✔ Range parameter doesn't look unreasonably large
s
#> $hessian_ok
#> [1] TRUE
#> 
#> $eigen_values_ok
#> [1] TRUE
#> 
#> $nlminb_ok
#> [1] TRUE
#> 
#> $range_ok
#> [1] TRUE
#> 
#> $gradients_ok
#> [1] TRUE
#> 
#> $se_magnitude_ok
#> [1] TRUE
#> 
#> $se_na_ok
#> [1] TRUE
#> 
#> $sigmas_ok
#> [1] TRUE
#> 
#> $all_ok
#> [1] TRUE
#> 

# If fitting many models in a loop, you may want to wrap
# sdmTMB() in try() to handle errors. sanity() will take an object
# of class "try-error" and return FALSE.
# Here, we will use stop() to simulate a failed sdmTMB() fit:
failed_fit <- try(stop())
#> Error in try(stop()) : 
s2 <- sanity(failed_fit)
all(unlist(s))
#> [1] TRUE
all(unlist(s2))
#> [1] FALSE
```
