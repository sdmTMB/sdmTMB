# Calculate conditional AIC

Calculates the conditional Akaike Information criterion (cAIC).

## Usage

``` r
cAIC(object, what = c("cAIC", "EDF"), ...)
```

## Arguments

- object:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- what:

  Whether to return the cAIC or the effective degrees of freedom (EDF)
  for each group of random effects.

- ...:

  Other arguments for specific methods. Not used.

## Value

Either the cAIC or the effective degrees of freedom (EDF) by group of
random effects depending on the argument `what`.

## Details

cAIC is designed to optimize the expected out-of-sample predictive
performance for new data that share the same random effects as the
in-sample (fitted) data, e.g., spatial interpolation. In this sense, it
should be a fast approximation to optimizing the model structure based
on k-fold cross-validation.

By contrast, [`AIC()`](https://rdrr.io/r/stats/AIC.html) calculates the
marginal Akaike Information Criterion, which is designed to optimize
expected predictive performance for new data that have new random
effects, e.g., extrapolation, or inference about generative parameters.

cAIC also calculates the effective degrees of freedom (EDF) as a
byproduct. This is the number of fixed effects that would have an
equivalent impact on model flexibility as a given random effect.

Both cAIC and EDF are calculated using Eq. 6 of Zheng, Cadigan, and
Thorson (2024).

For models that include profiled fixed effects, these profiles are
turned off.

## References

**Deriving the general approximation to cAIC used here:**

Zheng, N., Cadigan, N., & Thorson, J. T. (2024). A note on numerical
evaluation of conditional Akaike information for nonlinear mixed-effects
models (arXiv:2411.14185). arXiv.
[doi:10.48550/arXiv.2411.14185](https://doi.org/10.48550/arXiv.2411.14185)

**The utility of EDF to diagnose hierarchical model behaviour:**

Thorson, J. T. (2024). Measuring complexity for hierarchical models
using effective degrees of freedom. Ecology, 105(7), e4327
[doi:10.1002/ecy.4327](https://doi.org/10.1002/ecy.4327)

## Examples

``` r
mesh <- make_mesh(dogfish, c("X", "Y"), cutoff = 15)
fit <- sdmTMB(catch_weight ~ s(log(depth)),
  time_varying = ~1,
  time_varying_type = "ar1",
  time = "year",
  spatiotemporal = "off",
  mesh = mesh,
  family = tweedie(),
  data = dogfish,
  offset = log(dogfish$area_swept)
)
#> Detected irregular time spacing with an AR(1) or random walk process.
#> Consider filling in the missing time slices with `extra_time`.
#> `extra_time = c(2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2020)`
cAIC(fit)
#> [1] 12071.43
cAIC(fit, what = "EDF")
#> __s(log(depth))          b_rw_t         omega_s 
#>        6.613376        9.043457       38.730229 
AIC(fit)
#> [1] 12192.96
```
