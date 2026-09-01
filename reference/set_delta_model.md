# Set delta model for [`ggeffects::ggpredict()`](https://strengejacke.github.io/ggeffects/reference/ggpredict.html)

Set a delta model component to predict from with
[`ggeffects::ggpredict()`](https://strengejacke.github.io/ggeffects/reference/ggpredict.html).

## Usage

``` r
set_delta_model(x, model = c(NA, 1, 2))
```

## Arguments

- x:

  An [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
  model fit with a delta family such as
  [`delta_gamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md).

- model:

  Which delta/hurdle linear predictor to predict/plot with. `NA` does
  the combined prediction, `1` does the binomial part, and `2` does the
  positive part.

## Value

The fitted model with a new attribute named `delta_model_predict`. We
suggest you use `set_delta_model()` in a pipe (as in the examples) so
that this attribute does not persist. Otherwise,
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
will choose this model component by default. You can also remove the
attribute yourself after:

    attr(fit, "delta_model_predict") <- NULL

## Details

A complete version of the examples below would be:

    fit <- sdmTMB(density ~ poly(depth_scaled, 2), data = pcod_2011,
      spatial = "off", family = delta_gamma())

    # binomial part:
    set_delta_model(fit, model = 1) |>
      ggeffects::ggpredict("depth_scaled [all]")

    # gamma part:
    set_delta_model(fit, model = 2) |>
      ggeffects::ggpredict("depth_scaled [all]")

    # combined:
    set_delta_model(fit, model = NA) |>
      ggeffects::ggpredict("depth_scaled [all]")

But cannot be run on CRAN until the next version of ggeffects is
available on CRAN. For now, you can install the GitHub version of
ggeffects. <https://github.com/strengejacke/ggeffects>.

## Examples

``` r
fit <- sdmTMB(density ~ poly(depth_scaled, 2), data = pcod_2011,
  spatial = "off", family = delta_gamma())

# binomial part:
set_delta_model(fit, model = 1)
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ poly(depth_scaled, 2)
#> Mesh:  (isotropic covariance)
#> Data: pcod_2011
#> Family: delta_gamma(link1 = 'logit', link2 = 'log')
#> 
#> Delta/hurdle model 1: -----------------------------------
#> Family: binomial(link = 'logit') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)               -0.48    0.09
#> poly(depth_scaled, 2)1   -23.06    3.15
#> poly(depth_scaled, 2)2   -48.79    4.45
#> 
#> 
#> Delta/hurdle model 2: -----------------------------------
#> Family: Gamma(link = 'log') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)                4.24    0.08
#> poly(depth_scaled, 2)1    -5.49    3.50
#> poly(depth_scaled, 2)2   -13.26    3.23
#> 
#> Dispersion parameter: 0.64
#> 
#> ML criterion at convergence: 2936.579
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.

# gamma part:
set_delta_model(fit, model = 2)
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ poly(depth_scaled, 2)
#> Mesh:  (isotropic covariance)
#> Data: pcod_2011
#> Family: delta_gamma(link1 = 'logit', link2 = 'log')
#> 
#> Delta/hurdle model 1: -----------------------------------
#> Family: binomial(link = 'logit') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)               -0.48    0.09
#> poly(depth_scaled, 2)1   -23.06    3.15
#> poly(depth_scaled, 2)2   -48.79    4.45
#> 
#> 
#> Delta/hurdle model 2: -----------------------------------
#> Family: Gamma(link = 'log') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)                4.24    0.08
#> poly(depth_scaled, 2)1    -5.49    3.50
#> poly(depth_scaled, 2)2   -13.26    3.23
#> 
#> Dispersion parameter: 0.64
#> 
#> ML criterion at convergence: 2936.579
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.

# combined:
set_delta_model(fit, model = NA)
#> Model fit by ML ['sdmTMB']
#> Formula: density ~ poly(depth_scaled, 2)
#> Mesh:  (isotropic covariance)
#> Data: pcod_2011
#> Family: delta_gamma(link1 = 'logit', link2 = 'log')
#> 
#> Delta/hurdle model 1: -----------------------------------
#> Family: binomial(link = 'logit') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)               -0.48    0.09
#> poly(depth_scaled, 2)1   -23.06    3.15
#> poly(depth_scaled, 2)2   -48.79    4.45
#> 
#> 
#> Delta/hurdle model 2: -----------------------------------
#> Family: Gamma(link = 'log') 
#> Conditional model:
#>                        coef.est coef.se
#> (Intercept)                4.24    0.08
#> poly(depth_scaled, 2)1    -5.49    3.50
#> poly(depth_scaled, 2)2   -13.26    3.23
#> 
#> Dispersion parameter: 0.64
#> 
#> ML criterion at convergence: 2936.579
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.
```
