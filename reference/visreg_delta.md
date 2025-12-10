# Plot sdmTMB models with the visreg package

sdmTMB models fit with regular (non-delta) families can be passed to
[`visreg::visreg()`](https://pbreheny.github.io/visreg/reference/visreg.html)
or
[`visreg::visreg2d()`](https://pbreheny.github.io/visreg/reference/visreg2d.html)
directly. Examples are shown below. Delta models can use the helper
functions `visreg_delta()` or `visreg2d_delta()` described here.

## Usage

``` r
visreg_delta(object, ..., model = c(1, 2))

visreg2d_delta(object, ..., model = c(1, 2))
```

## Arguments

- object:

  Fit from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)

- ...:

  Any arguments passed to
  [`visreg::visreg()`](https://pbreheny.github.io/visreg/reference/visreg.html)
  or
  [`visreg::visreg2d()`](https://pbreheny.github.io/visreg/reference/visreg2d.html)

- model:

  1st or 2nd delta model

## Value

A plot from the visreg package. Optionally, the data plotted invisibly
if `plot = FALSE`. This is useful if you want to make your own plot
after.

## Details

Note the residuals are currently randomized quantile residuals, *not*
deviance residuals as is usual for GLMs with visreg.

## Examples

``` r
if (require("ggplot2", quietly = TRUE) &&
  require("visreg", quietly = TRUE)) {

# \donttest{
  fit <- sdmTMB(
    density ~ s(depth_scaled),
    data = pcod_2011,
    spatial = "off",
    family = tweedie()
  )
  visreg::visreg(fit, xvar = "depth_scaled")

  visreg::visreg(fit, xvar = "depth_scaled", scale = "response")
  v <- visreg::visreg(fit, xvar = "depth_scaled")
  head(v$fit)
  # now use ggplot2 etc. if desired

  # Delta model example:
  fit_dg <- sdmTMB(
    density ~ s(depth_scaled, year, k = 8),
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatial = "off",
    family = delta_gamma()
  )
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 1, gg = TRUE)
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 2, gg = TRUE)
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 1,
    scale = "response", gg = TRUE
  )
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response"
  )
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response", gg = TRUE, rug = FALSE
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response"
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 1, scale = "response", plot.type = "persp"
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response", plot.type = "gg"
  )
  # }
}



#> These are residuals for delta model component 1. Use the `model` argument to
#> select the other component.
#> These are residuals for delta model component 2. Use the `model` argument to
#> select the other component.
#> These are residuals for delta model component 1. Use the `model` argument to
#> select the other component.
#> These are residuals for delta model component 2. Use the `model` argument to
#> select the other component.

#> These are residuals for delta model component 2. Use the `model` argument to
#> select the other component.


```
