# Calculate conditional AIC

Calculates the conditional Akaike Information criterion (cAIC).

## Usage

``` r
cAIC(object, ..., what = c("cAIC", "EDF"))
```

## Arguments

- object:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) or
  [`profile_breakpt()`](https://sdmTMB.github.io/sdmTMB/reference/profile_breakpt.md).

- ...:

  Optionally, more fitted models from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) or
  [`profile_breakpt()`](https://sdmTMB.github.io/sdmTMB/reference/profile_breakpt.md)
  to compare, as with [`AIC()`](https://rdrr.io/r/stats/AIC.html).

- what:

  Whether to return the cAIC or the effective degrees of freedom (EDF)
  for each group of random effects. Must be named. `"EDF"` is only
  available for a single model.

## Value

For a single model, either the cAIC or the effective degrees of freedom
(EDF) by group of random effects depending on the argument `what`. For
several models, a data frame with a row per model and columns `df` (the
effective degrees of freedom: fixed effects plus total EDF) and `cAIC`,
as with [`AIC()`](https://rdrr.io/r/stats/AIC.html).

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

When comparing several models, models without random effects contribute
their marginal AIC, which is what cAIC reduces to with no random
effects.

For output `p` from
[`profile_breakpt()`](https://sdmTMB.github.io/sdmTMB/reference/profile_breakpt.md),
cAIC uses the final fit and counts the profiled cutpoint as one more
fixed effect, as `AIC(p)` does. Use `cAIC(p)`, not `cAIC(p$fit)`:
`p$fit` holds the cutpoint fixed, so it isn't counted. The EDF
(`what = "EDF"`) are those of `p$fit`.

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
#> Detected irregular time spacing with an AR(1), random walk, or time-indexed
#> nonlocal process.
#> Consider filling in the missing time slices with `extra_time`.
#> `extra_time = c(2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2020)`
cAIC(fit)
#> [1] 12071.43
cAIC(fit, what = "EDF")
#> __s(log(depth))          b_rw_t         omega_s 
#>        6.615134        9.043454       38.730196 
AIC(fit)
#> [1] 12192.96

fit2 <- update(fit, catch_weight ~ 1)
#> Detected irregular time spacing with an AR(1), random walk, or time-indexed
#> nonlocal process.
#> Consider filling in the missing time slices with `extra_time`.
#> `extra_time = c(2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2020)`
cAIC(fit, fit2)
#>            df     cAIC
#> fit  63.38878 12071.43
#> fit2 58.61878 12184.28
```
