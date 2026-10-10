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
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- ...:

  Arguments passed to
  [`visreg::visreg()`](https://pbreheny.github.io/visreg/reference/visreg.html)
  or
  [`visreg::visreg2d()`](https://pbreheny.github.io/visreg/reference/visreg2d.html).

- model:

  First or second delta-model component.

## Value

A plot from the visreg package. If `plot = FALSE`, the plotted data are
returned invisibly. This is useful if you want to make your own plot.

## Details

Note that the residuals are currently randomized quantile residuals,
*not* the deviance residuals usually used for GLMs with visreg.

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
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 1)
  visreg_delta(fit_dg, xvar = "depth_scaled", model = 2)
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 1,
    scale = "response"
  )
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response"
  )
  visreg_delta(fit_dg,
    xvar = "depth_scaled", model = 2,
    scale = "response", rug = FALSE
  )
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response"
  )
  v2d <- visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 1, scale = "response", plot = FALSE
  )
  if (!is.null(utils::getS3method("persp", "visreg2d", optional = TRUE))) {
    persp(v2d)
  }
  visreg2d_delta(fit_dg,
    xvar = "depth_scaled", yvar = "year",
    model = 2, scale = "response"
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
