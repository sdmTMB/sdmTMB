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
#> 1          1     2.771541   -0.7828212 41.11743 14.69055  1.590188 1.432009
#> 2          2     2.540419   -0.6472292 42.70610 14.03787  1.607405 1.975027
#> 3          3     3.110238   -0.4588569 35.97312 15.42229  1.560272 2.403535
#> 4          4     3.177685   -0.3777292 25.74392 15.20962  1.608361 2.257550
#> 5          5     2.750845   -0.8808476 26.61882 14.00816  1.600663 2.057045
#> 6          6     3.015164   -0.8293287 16.39538 15.38875  1.572145 2.849170
head(gather_sims(m, nsim = 10))
#>   .iteration    .variable   .value
#> 1          1 X.Intercept. 3.159876
#> 2          2 X.Intercept. 3.452272
#> 3          3 X.Intercept. 2.649328
#> 4          4 X.Intercept. 3.272258
#> 5          5 X.Intercept. 2.843023
#> 6          6 X.Intercept. 2.729322
samps <- gather_sims(m, nsim = 1000)

if (require("ggplot2", quietly = TRUE)) {
  ggplot(samps, aes(.value)) + geom_histogram() +
    facet_wrap(~.variable, scales = "free_x")
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
