# Add UTM coordinates to a data frame

Add UTM (Universal Transverse Mercator) coordinates to a data frame.
This is useful since geostatistical modeling should generally be
performed in an equal-distance projection. You can do this yourself
separately with the
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html),
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html),
and
[`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
functions in the sf package.

## Usage

``` r
add_utm_columns(
  dat,
  ll_names = c("longitude", "latitude"),
  ll_crs = 4326,
  utm_names = c("X", "Y"),
  utm_crs = get_crs(dat, ll_names),
  units = c("km", "m")
)

get_crs(dat, ll_names = c("longitude", "latitude"))
```

## Arguments

- dat:

  Data frame that contains longitude and latitude columns.

- ll_names:

  Longitude and latitude column names. **Note the order.**

- ll_crs:

  Input CRS value for `ll_names`.

- utm_names:

  Output column names for the UTM columns.

- utm_crs:

  Output CRS value for the UTM zone; tries to detect with `get_crs()`
  but can be specified manually.

- units:

  UTM units.

## Value

A copy of the input data frame with new columns for UTM coordinates.

## Details

**Note that longitudes west of the prime meridian should be encoded as
running from -180 to 0 degrees.**

You may wish to work in km's rather than the standard UTM meters so that
the range parameter estimate is not too small, which can cause
computational issues. This depends on the the scale of your data.

## Examples

``` r
d <- data.frame(lat = c(52.1, 53.4), lon = c(-130.0, -131.4))
get_crs(d, c("lon", "lat"))
#> Detected UTM zone 9N; CRS = 32609.
#> Visit https://epsg.io/32609 to verify.
#> [1] 32609
add_utm_columns(d, c("lon", "lat"))
#> Detected UTM zone 9N; CRS = 32609.
#> Visit https://epsg.io/32609 to verify.
#>    lat    lon        X        Y
#> 1 52.1 -130.0 431.5034 5772.632
#> 2 53.4 -131.4 340.4411 5919.452
```
