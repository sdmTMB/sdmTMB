# Compare deviance between two models

Calculates the proportion of deviance explained by a model relative to a
simpler model, usually an intercept-only model. This is a
pseudo-\\R^2\\. Any parameters that the deviance depends on, other than
the mean, are held fixed at their values from the first model.

## Usage

``` r
compare_deviance(object, reduced)
```

## Arguments

- object:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). The
  model of interest. Its shape parameters (see Details) are used for
  both models.

- reduced:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) fit
  to the same response with the same family. For a pseudo-\\R^2\\, an
  intercept-only model with no random effects or random fields (see
  Details).

## Value

A data frame with the deviance of each model, their difference
(`reduced` minus `object`), and the proportion of the deviance of
`reduced` explained by `object`. The refitted `reduced` model is
included as attribute `"reduced"`, and the fixed parameter values (on
their internal scale) as attribute `"fixed"`.

## Details

The proportion of deviance explained, \\1 - D / D_0\\, is a common
pseudo-\\R^2\\ for generalized linear models (Cameron and Windmeijer
1997). It is the proportion of the Kullback-Leibler divergence between
the data and the null model that the fitted model removes. It reduces to
the usual \\R^2\\ for a Gaussian linear model and is the "deviance
explained" reported by mgcv.

**Shape parameters.** The deviance is twice the log-likelihood
difference between a saturated model and the fitted model. For some
families the saturated model, and hence the deviance, depends on
parameters other than the mean:

- negative binomial families: the dispersion parameter `phi`,

- [`tweedie()`](https://sdmTMB.github.io/sdmTMB/reference/families.md):
  the power parameter `p` (but not `phi`),

- [`lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md):
  `phi` (through the bias correction to the mean),

- [`gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md):
  `phi` and `Q`.

If these parameters are estimated separately for each model, the two
deviances are measured against different saturated models and their
ratio is not meaningful. A simpler model typically estimates more
overdispersion (e.g., a smaller NB2 `phi`), which lowers its deviance
and understates the deviance explained by the more complex model. The
naive ratio can even be negative.

This function therefore refits `reduced` with these parameters fixed at
the values estimated in `object`. This matches how mgcv computes the
null deviance for its negative binomial and Tweedie families. Families
whose deviance depends only on the mean (e.g.,
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
[`Gamma()`](https://rdrr.io/r/stats/family.html),
[`gaussian()`](https://rdrr.io/r/stats/family.html)) need no refit.
Because `object` supplies these parameters, swapping the arguments gives
a different answer.

**Random effects and random fields.** For models with random effects,
including spatial and spatiotemporal random fields, the deviance is
conditional on the predicted random effects. The deviance explained is
therefore a conditional pseudo-\\R^2\\: the in-sample variation captured
by the fixed effects and the predicted random effects together, similar
in spirit to the conditional \\R^2\\ for mixed models. It is a
reasonable summary of in-sample fit, but note that:

- Like \\R^2\\, it does not penalize complexity. More flexible random
  fields (e.g., adding spatiotemporal fields or a finer mesh) almost
  always increase it, and a sufficiently flexible field can approach 1.
  Do not use it to choose among random effect structures, random field
  structures, or meshes. Use
  [`cAIC()`](https://sdmTMB.github.io/sdmTMB/reference/cAIC.md) (the
  conditional likelihood with a penalty for effective degrees of
  freedom), [`AIC()`](https://rdrr.io/r/stats/AIC.html), or
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
  instead.

- It describes fit to the observed data, not predictive performance for
  new locations or times.

- `reduced` should not include random effects or random fields. Shape
  parameters such as `phi` also determine how variation is partitioned
  between random effects and observation error, so fixing them in a
  `reduced` model that has random effects changes its random effect
  estimates and makes the result hard to interpret.

- The difference in deviance does not follow a chi-squared distribution
  and is not a likelihood-ratio test.

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
#> Refitting `reduced` with Tweedie p fixed at the value from `object`, since the
#> deviance depends on it.
#>   deviance deviance_reduced difference deviance_explained
#> 1 14590.31         17730.19   3139.878          0.1770922

# Not comparable, since each model estimates its own Tweedie p:
1 - deviance(fit) / deviance(fit0)
#> [1] 0.114655

# Conditional pseudo-R^2 including a spatial random field:
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit_sp <- sdmTMB(density ~ depth_scaled + depth_scaled2,
  data = pcod_2011, mesh = mesh, family = tweedie())
compare_deviance(fit_sp, fit0)
#> Refitting `reduced` with Tweedie p fixed at the value from `object`, since the
#> deviance depends on it.
#>   deviance deviance_reduced difference deviance_explained
#> 1 11880.74         18753.27   6872.532          0.3664711
```
