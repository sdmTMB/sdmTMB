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
#> 1          1     2.951100   -0.6300220 24.16212 14.94730  1.579240 2.567088
#> 2          2     2.862951   -0.7076326 43.50890 15.65875  1.583221 1.943565
#> 3          3     2.568716   -0.5655881 75.22831 13.87196  1.582173 1.227205
#> 4          4     3.306620   -0.5419931 26.18201 16.21668  1.588701 1.671767
#> 5          5     2.871752   -0.7513382 34.94243 15.06824  1.595574 2.008109
#> 6          6     2.903716   -0.7612297 26.66178 14.38083  1.614983 2.700071
head(gather_sims(m, nsim = 10))
#>   .iteration    .variable   .value
#> 1          1 X.Intercept. 3.121410
#> 2          2 X.Intercept. 3.072362
#> 3          3 X.Intercept. 2.692724
#> 4          4 X.Intercept. 3.011796
#> 5          5 X.Intercept. 3.123790
#> 6          6 X.Intercept. 2.508618
samps <- gather_sims(m, nsim = 1000)

if (require("ggplot2", quietly = TRUE)) {
  ggplot(samps, aes(.value)) + geom_histogram() +
    facet_wrap(~.variable, scales = "free_x")
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
