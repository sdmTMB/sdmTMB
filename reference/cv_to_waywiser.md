# Convert [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md) objects to sf format for spatial assessment with waywiser

**\[experimental\]** Converts cross-validation results to an
[`sf::sf()`](https://r-spatial.github.io/sf/reference/sf.html) object
for use with spatial model assessment tools such as those in the
waywiser package. This enables multi-scale spatial assessment of model
predictions.

## Usage

``` r
cv_to_waywiser(
  object,
  ll_names = c("longitude", "latitude"),
  ll_crs = 4326,
  utm_crs = get_crs(object$data, ll_names)
)
```

## Arguments

- object:

  An object of class `sdmTMB_cv` from
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md).

- ll_names:

  Column names for longitude and latitude in the original data. **Note
  the order: longitude first, then latitude.**

- ll_crs:

  The coordinate reference system (CRS) for the `ll_names` columns.
  Defaults to `4326` (WGS84 lon/lat).

- utm_crs:

  The projected coordinate reference system (CRS) for the output sf
  object. By default (if you're feeling lucky!) automatically detected
  using
  [`get_crs()`](https://sdmTMB.github.io/sdmTMB/reference/add_utm_columns.md)
  based on `ll_names`. Can be manually specified as an EPSG code (e.g.,
  `32609`) or any format accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).

## Value

An [`sf::sf()`](https://r-spatial.github.io/sf/reference/sf.html) object
with POINT geometry containing:

- truth:

  The observed response values

- estimate:

  The cross-validated predictions

- geometry:

  Spatial point locations

## Details

This function is particularly useful for assessing spatial models at
multiple scales using the waywiser package. After converting to sf
format, you can use functions like
[`waywiser::ww_multi_scale()`](https://docs.ropensci.org/waywiser/reference/ww_multi_scale.html)
to evaluate how model performance changes when predictions are
aggregated to different spatial scales.

For delta/hurdle models, the combined predictions are returned (e.g.,
the product of the encounter probability and positive catch rate).

## See also

[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md),
[`get_crs()`](https://sdmTMB.github.io/sdmTMB/reference/add_utm_columns.md),
<https://sdmTMB.github.io/sdmTMB/articles/cross-validation.html>

## Examples

``` r
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 12)

# Run cross-validation
set.seed(123)
m_cv <- sdmTMB_cv(
  density ~ depth_scaled,
  data = pcod_2011,
  mesh = mesh,
  family = tweedie(),
  k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# Convert with default auto-detected CRS based on lon/lat columns:
cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"))
#> Detected UTM zone 9N; CRS = 32609.
#> Visit https://epsg.io/32609 to verify.

# Or manually specify the desired UTM CRS:
cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"), utm_crs = 32609)

# Use with waywiser for multi-scale assessment
waywiser::ww_multi_scale(
  cv_sf,
  truth,    # column name (unquoted)
  estimate, # column name (unquoted)
  n = list(c(5, 5), c(2, 2)) # 5x5 and 2x2 grid blocks
)
#> # A tibble: 4 × 6
#>   .metric .estimator .estimate .grid_args       .grid         .notes          
#>   <chr>   <chr>          <dbl> <list>           <list>        <list>          
#> 1 rmse    standard       10.8  <tibble [1 × 1]> <sf [25 × 5]> <tibble [0 × 2]>
#> 2 mae     standard        7.87 <tibble [1 × 1]> <sf [25 × 5]> <tibble [0 × 2]>
#> 3 rmse    standard        5.70 <tibble [1 × 1]> <sf [4 × 5]>  <tibble [0 × 2]>
#> 4 mae     standard        4.93 <tibble [1 × 1]> <sf [4 × 5]>  <tibble [0 × 2]>
```
