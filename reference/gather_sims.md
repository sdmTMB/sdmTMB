# Extract parameter simulations from the joint precision matrix

`spread_sims()` returns a wide-format data frame. `gather_sims()`
returns a long-format data frame. The format matches the format in the
tidybayes `spread_draws()` and `gather_draws()` functions.

## Usage

``` r
spread_sims(object, nsim = 200)

gather_sims(object, nsim = 200)
```

## Arguments

- object:

  Output from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- nsim:

  The number of simulation draws.

## Value

A data frame. `gather_sims()` returns a long-format data frame:

- `.iteration`: the sample ID

- `.variable`: the parameter name

- `.value`: the parameter sample value

`spread_sims()` returns a wide-format data frame:

- `.iteration`: the sample ID

- columns for each parameter with a sample per row

## Examples

``` r
m <- sdmTMB(density ~ depth_scaled,
  data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie())
head(spread_sims(m, nsim = 10))
#>   .iteration X.Intercept. depth_scaled    range      phi tweedie_p  sigma_O
#> 1          1     2.109968   -0.6737406 21.78563 15.20790  1.600079 2.462644
#> 2          2     2.414470   -0.4543877 22.34274 14.84304  1.562145 3.497004
#> 3          3     2.971140   -0.7952331 23.69878 15.41066  1.614860 2.188129
#> 4          4     2.655215   -0.6236876 36.94339 15.45250  1.586901 1.954655
#> 5          5     2.857308   -0.4080048 62.60743 14.88711  1.581576 1.959511
#> 6          6     3.111335   -0.6012389 42.53989 16.30831  1.573077 1.540739
head(gather_sims(m, nsim = 10))
#>   .iteration    .variable   .value
#> 1          1 X.Intercept. 3.094201
#> 2          2 X.Intercept. 3.024605
#> 3          3 X.Intercept. 2.703684
#> 4          4 X.Intercept. 3.114728
#> 5          5 X.Intercept. 3.254302
#> 6          6 X.Intercept. 3.406452
samps <- gather_sims(m, nsim = 1000)

if (require("ggplot2", quietly = TRUE)) {
  ggplot(samps, aes(.value)) + geom_histogram() +
    facet_wrap(~.variable, scales = "free_x")
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
