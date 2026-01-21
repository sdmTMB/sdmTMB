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
#> 1          1     2.190191   -0.5833938 30.48107 16.16591  1.566084 2.456067
#> 2          2     2.603178   -0.7394658 17.31691 16.46800  1.598118 3.104577
#> 3          3     2.688482   -1.0642911 36.94842 15.06568  1.607437 1.908867
#> 4          4     2.808461   -0.6082838 27.88899 15.28096  1.600632 2.433739
#> 5          5     2.945574   -0.6886035 39.42264 15.79456  1.595986 2.384444
#> 6          6     3.218579   -0.7272557 54.07811 15.34240  1.613165 1.701300
head(gather_sims(m, nsim = 10))
#>   .iteration    .variable   .value
#> 1          1 X.Intercept. 3.421495
#> 2          2 X.Intercept. 2.847998
#> 3          3 X.Intercept. 2.929768
#> 4          4 X.Intercept. 2.815933
#> 5          5 X.Intercept. 2.103157
#> 6          6 X.Intercept. 2.806020
samps <- gather_sims(m, nsim = 1000)

if (require("ggplot2", quietly = TRUE)) {
  ggplot(samps, aes(.value)) + geom_histogram() +
    facet_wrap(~.variable, scales = "free_x")
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
