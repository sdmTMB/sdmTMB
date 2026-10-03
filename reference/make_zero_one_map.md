# Fix factor-level intercepts for groups with all zeros or ones

A grouping factor level (e.g. factor(year), but also possibly a region,
stratum, region-by-year interaction, or any other factor level) with
all-zero observations can drive the corresponding factor-coded
coefficient to `-Inf`, causing convergence problems. In delta/hurdle and
binomial or beta-binomial models, an all-positive or all-success level
can similarly drive the encounter coefficient to `+Inf`. Tweedie will
have the same issue with all zeros.This helper returns `map` and `start`
lists suitable for passing through
[`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
to fix those problem coefficients at large positive/negative values.
This mirrors an approach taken in VAST.

## Usage

``` r
make_zero_one_map(formula, data, group, family, weights = NULL, value = 20)
```

## Arguments

- formula:

  The fixed-effect formula to be used in
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). For
  delta models, either a single formula applied to both linear
  predictors, or a [`list()`](https://rdrr.io/r/base/list.html) of two
  formulas. The formula should contain `factor(group)` (or
  `as.factor(group)`, or a pre-converted factor column) so each level
  has its own column, typically `response ~ 0 + factor(year)`.

- data:

  The data frame passed to
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- group:

  Character. Name of the grouping column in `data` whose levels should
  be inspected for all-zero/all-positive observations. Most commonly the
  time column.

- family:

  The `sdmTMB` family object (e.g.,
  [`tweedie()`](https://sdmTMB.github.io/sdmTMB/reference/families.md),
  [`delta_gamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)).

- weights:

  Optional numeric vector of trial sizes for non-delta binomial or
  beta-binomial models fit with proportion data. This should match the
  `weights` argument passed to
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).
  Ignored for other families.

- value:

  Positive numeric. Magnitude at which to fix the intercept;
  `plogis(20)`, `plogis(-20)`, and `exp(-20)` are all within ~2e-9 of 0
  or 1, which is plenty.

## Value

A named list with elements:

- `map`:

  A list with `b_j` (and `b_j2` for delta families) factors suitable for
  `sdmTMBcontrol(map = ...)`.

- `start`:

  A list with `b_j` (and `b_j2`) starting-value vectors suitable for
  `sdmTMBcontrol(start = ...)`.

- `all_zero_levels`:

  Group levels with all-zero observations.

- `all_one_levels`:

  Group levels with all-positive observations for delta families, or
  all-success observations for non-delta binomial/beta-binomial
  families. Empty for other families.

## Details

Inspect the returned `all_zero_levels` and `all_one_levels` elements
before fitting to make sure the detected problem levels match your
expectation.

Although the motivating use case is years in fisheries index
standardization (e.g., `density ~ 0 + factor(year)`), the same mechanism
applies to any factor whose levels each get a column in the design
matrix. For an interaction such as region by year, create a combined
column (e.g.,
`data$region_year <- interaction(data$region, data$year, drop = TRUE)`)
and use `~ 0 + region_year`.

For Poisson-link delta families (e.g.,
`delta_gamma(type = "poisson-link")`), the first linear predictor
controls both the encounter probability and the positive rate, so fixing
it shifts the meaning of the second-linear-predictor coefficients in
those levels. The helper still works for all-zero levels (the large
negative value makes encounter probability ~0, and the level contributes
~0 to the index), but emits a message for all-positive levels. Although
this affects interpretation of the coefficients, the overall combined
prediction remains OK.

## See also

[`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md),
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)

## Examples

``` r
d <- pcod
d$density[d$year == 2003] <- 0 # force an all-zero year

m <- make_zero_one_map(
  formula = density ~ 0 + factor(year),
  data = d,
  group = "year",
  family = tweedie()
)
m$all_zero_levels
#> [1] 2003
m$all_one_levels
#> integer(0)
m$start
#> $b_j
#> [1] -20   0   0   0   0   0   0   0   0
#> 
m$map
#> $b_j
#> [1] <NA> 2    3    4    5    6    7    8    9   
#> Levels: 2 3 4 5 6 7 8 9
#> 

# \donttest{
mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
fit <- sdmTMB(
  density ~ 0 + factor(year),
  data = d,
  mesh = mesh,
  time = "year",
  family = tweedie(),
  spatial = "on",
  spatiotemporal = "off", # off for example speed
  control = sdmTMBcontrol(map = m$map, start = m$start) #<
)
#> ℹ Initiating `b_j` at specified starting value(s) of:
#> -20, 0, 0, 0, 0, 0, 0, 0, 0
#> ℹ Fixing or mirroring `b_j`
fit
#> Spatial model fit by ML ['sdmTMB']
#> Formula: density ~ 0 + factor(year)
#> Mesh: mesh (isotropic covariance)
#> Time column: character
#> Data: d
#> Family: tweedie(link = 'log')
#>  
#> Conditional model:
#>                  coef.est coef.se
#> factor(year)2003   -20.00      NA
#> factor(year)2004     3.48    0.40
#> factor(year)2005     3.39    0.40
#> factor(year)2007     2.04    0.41
#> factor(year)2009     2.47    0.40
#> factor(year)2011     3.12    0.40
#> factor(year)2013     2.99    0.40
#> factor(year)2015     3.17    0.40
#> factor(year)2017     2.61    0.40
#> 
#> Dispersion parameter: 14.55
#> Tweedie p: 1.61
#> Matérn range: 29.25
#> Spatial SD: 2.12
#> ML criterion at convergence: 5926.241
#> 
#> See ?tidy.sdmTMB to extract these values as a data frame.
coef(fit)
#> factor(year)2003 factor(year)2004 factor(year)2005 factor(year)2007 
#>       -20.000000         3.483091         3.390535         2.040455 
#> factor(year)2009 factor(year)2011 factor(year)2013 factor(year)2015 
#>         2.469601         3.115998         2.992516         3.172495 
#> factor(year)2017 
#>         2.610296 
# }
```
