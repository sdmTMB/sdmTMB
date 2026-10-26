# Compare deviance among models

Calculates the proportion of deviance explained by one or more models
relative to a reduced model (a pseudo-\\R^2\\), and the gain from each
added term or component, such as a spatial random field.

## Usage

``` r
compare_deviance(object, ...)
```

## Arguments

- object:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). The
  most complex model. Its shape parameters (see Details) are used for
  all models.

- ...:

  Other
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
  models fit to the same response with the same family, ordered from
  most to least complex. The last is the reduced model, usually
  intercept-only with no random effects or random fields. Name the
  arguments to label the output rows. Names are arbitrary and only to
  help reading the output. These models can be set up with
  `do_fit = FALSE` and are then fit here. Alternatively, they can be fit
  first, in which case they are refit (if needed, see below) starting
  from the previously estimated parameter values for speed.

## Value

A data frame with one row per model with columns `model`, `deviance`,
`deviance_explained` (relative to the last model), and `gain` (over the
next simpler model). Any fixed shape parameter values (on their internal
scale) are in attribute `"fixed"`.

## Details

Deviance explained, \\1 - D / D_0\\, is a common pseudo-\\R^2\\ for
generalized linear models (Cameron and Windmeijer 1997) and is what mgcv
reports. It reduces to \\R^2\\ for a Gaussian linear model.

**Sequences of models.** `gain` is the increase in deviance explained
over the next simpler model. Like sequential (type I) sums of squares,
it depends on order: a spatial field added after covariates usually
gains less than one added to an intercept-only model, since both can
explain the same spatial pattern.

**Shape parameters.** For some families, the deviance depends on
parameters other than the mean: `phi` for negative binomial,
[`compois()`](https://sdmTMB.github.io/sdmTMB/reference/families.md),
and
[`lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md);
`p` for
[`tweedie()`](https://sdmTMB.github.io/sdmTMB/reference/families.md);
and `phi` and `Q` for
[`gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md).
Deviances computed with different values of these are not comparable, so
the other models are refit with them fixed at their values in `object`,
as mgcv does for its null deviance. Previously fit models are refit
starting from their estimates for speed. Reordering the models changes
the answer because it changes the comparisons.

**Random effects and random fields.** The deviance is conditional on the
predicted random effects, so deviance explained measures in-sample fit
of the fixed and random effects together. Because it does not penalize
complexity, more flexible fields almost always increase it. Do not use
it to choose among random effect structures or meshes; use
[`cAIC()`](https://sdmTMB.github.io/sdmTMB/reference/cAIC.md),
[`AIC()`](https://rdrr.io/r/stats/AIC.html), or
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
instead. For families with shape parameters, fixing `phi` can shift
variation between random effects and observation error in refit models,
so gains from random effects are approximate.

For delta models, the deviance is summed over both components.

## References

Cameron, A. C., and Windmeijer, F. A. G. (1997). An R-squared measure of
goodness of fit for some common nonlinear regression models. Journal of
Econometrics, 77(2), 329–342.
[doi:10.1016/S0304-4076(96)01818-0](https://doi.org/10.1016/S0304-4076%2896%2901818-0)

## See also

[`stats::deviance()`](https://rdrr.io/r/stats/deviance.html),
[`cAIC()`](https://sdmTMB.github.io/sdmTMB/reference/cAIC.md)

## Examples

``` r
fit <- sdmTMB(density ~ depth_scaled + depth_scaled2,
  data = pcod_2011, spatial = "off", family = tweedie())
fit0 <- sdmTMB(density ~ 1,
  data = pcod_2011, spatial = "off", family = tweedie())
compare_deviance(fit, fit0)
#> Refitting the other models with Tweedie p fixed at the value from `object`,
#> since the deviance depends on it.
#>   model deviance deviance_explained      gain
#> 1   fit 14590.31          0.1770922 0.1770922
#> 2  fit0 17730.19          0.0000000        NA

# Not comparable, since each model estimates its own Tweedie p:
1 - deviance(fit) / deviance(fit0)
#> [1] 0.114655

# \donttest{
# Gain from adding a spatial random field after depth:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit_sp <- update(fit, spatial = "on", mesh = mesh)
compare_deviance(fit_sp, depth = fit, intercept = fit0)
#> Refitting the other models with Tweedie p fixed at the value from `object`,
#> since the deviance depends on it.
#>       model deviance deviance_explained      gain
#> 1    fit_sp 11880.74          0.3664711 0.1863758
#> 2     depth 15375.90          0.1800953 0.1800953
#> 3 intercept 18753.27          0.0000000        NA

# Models after the first can be left unfitted
# with `do_fit = FALSE`, since they are fit here with
# Tweedie p fixed anyway:
fit_sp_only <- update(fit0, spatial = "on", mesh = mesh, do_fit = FALSE)
compare_deviance(fit_sp, spatial = fit_sp_only, intercept = fit0)
#> Refitting the other models with Tweedie p fixed at the value from `object`,
#> since the deviance depends on it.
#>       model deviance deviance_explained      gain
#> 1    fit_sp 11880.74          0.3664711 0.1036665
#> 2   spatial 13824.83          0.2628045 0.2628045
#> 3 intercept 18753.27          0.0000000        NA
# }
```
