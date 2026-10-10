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
#> 1          1     2.903231   -0.5205483 28.06440 16.02173  1.569018 2.095025
#> 2          2     2.594193   -0.5693215 40.17931 16.17853  1.556120 1.826335
#> 3          3     3.197205   -0.8617446 47.15718 14.71530  1.597992 1.925214
#> 4          4     3.510531   -0.6283991 17.08680 15.70606  1.583781 2.350573
#> 5          5     2.713482   -0.7645729 24.53404 14.34188  1.572464 2.425859
#> 6          6     2.789035   -0.6301272 33.83201 14.27102  1.602406 1.769855
head(gather_sims(m, nsim = 10))
#>   .iteration    .variable   .value
#> 1          1 X.Intercept. 3.423201
#> 2          2 X.Intercept. 3.167973
#> 3          3 X.Intercept. 2.762833
#> 4          4 X.Intercept. 3.459236
#> 5          5 X.Intercept. 2.705132
#> 6          6 X.Intercept. 2.739150
samps <- gather_sims(m, nsim = 1000)

if (require("ggplot2", quietly = TRUE)) {
  ggplot(samps, aes(.value)) + geom_histogram() +
    facet_wrap(~.variable, scales = "free_x")
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
