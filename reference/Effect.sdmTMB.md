# Calculate effects

Used by effects package

## Usage

``` r
Effect.sdmTMB(focal.predictors, mod, ...)
```

## Arguments

- focal.predictors:

  a character vector of one or more predictors in the model in any
  order.

- mod:

  a regression model object. If no specific method exists for the class
  of `mod`, `Effect.default` will be called.

- ...:

  arguments to be passed down.

## Value

Output from
[`effects::effect()`](https://rdrr.io/pkg/effects/man/effect.html). Can
then be plotted with with associated
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

## Examples

``` r
fit <- sdmTMB(present ~ depth_scaled, data = pcod_2011, family = binomial(),
  spatial = "off")
effects::effect("depth_scaled", fit)
#> 
#>  depth_scaled effect
#> depth_scaled
#>        -3        -2      -0.1         1         3 
#> 0.7511661 0.6634557 0.4673280 0.3544365 0.1897228 
plot(effects::effect("depth_scaled", fit))
```
