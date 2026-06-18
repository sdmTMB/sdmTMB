# Plot PC Matérn priors

Plot PC Matérn priors

## Usage

``` r
plot_pc_matern(
  range_gt,
  sigma_lt,
  range_prob = 0.05,
  sigma_prob = 0.05,
  range_lims = c(range_gt * 0.1, range_gt * 10),
  sigma_lims = c(0, sigma_lt * 2),
  plot = TRUE
)
```

## Arguments

- range_gt:

  A value one expects the spatial or spatiotemporal range is **g**reater
  **t**han with `1 - range_prob` probability.

- sigma_lt:

  A value one expects the spatial or spatiotemporal marginal standard
  deviation (`sigma_O` or `sigma_E` internally) is **l**ess **t**han
  with `1 - sigma_prob` probability.

- range_prob:

  Probability. See description for `range_gt`.

- sigma_prob:

  Probability. See description for `sigma_lt`.

- range_lims:

  Plot range variable limits.

- sigma_lims:

  Plot sigma variable limits.

- plot:

  Logical controlling whether plot is drawn (defaults to `TRUE`).

## Value

A plot from [`image()`](https://rdrr.io/r/graphics/image.html).
Invisibly returns the underlying matrix data. The rows are the sigmas.
The columns are the ranges. Column and row names are provided.

## See also

[`pc_matern()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)

## Examples

``` r
plot_pc_matern(range_gt = 5, sigma_lt = 1)

plot_pc_matern(range_gt = 5, sigma_lt = 10)

plot_pc_matern(range_gt = 5, sigma_lt = 1, sigma_prob = 0.2)

plot_pc_matern(range_gt = 5, sigma_lt = 1, range_prob = 0.2)
```
