# Ohio lung cancer mortality data, reshaped for spatiotemporal examples

A derived version of `geodaData::ohio_lung`, reshaped into county-year
and county geometry formats for package examples and vignettes.

## Usage

``` r
ohio_df

ohio_sf
```

## Format

`ohio_df`: A county-year data frame with columns:

- county:

  County name.

- year:

  Observation year.

- cases:

  Observed lung cancer deaths.

- pop:

  Population at risk.

- pct_male:

  Proportion of the population that is male.

`ohio_sf`: An `sf` object with columns:

- county:

  County name.

- geometry:

  County polygon geometry.

## Source

Derived from `geodaData::ohio_lung`, originally distributed in the CRAN
package `geodaData` under CC0. The original dataset is based on the
GeoDa Center Ohio lung cancer mortality dataset.

## Examples

``` r
data(ohio_df)
head(ohio_df)
#>     county year cases    pop  pct_male
#> 1    Lucas 1968   147 231377 0.4764055
#> 2   Fulton 1968     4  15571 0.5154639
#> 3   Geauga 1968     4  30887 0.4876033
#> 4 Williams 1968     5  16276 0.5454545
#> 5 Cuyahoga 1968   522 834231 0.4706767
#> 6   Ottawa 1968    10  17802 0.4739884
data(ohio_sf)
head(ohio_sf)
#> Simple feature collection with 6 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 183138 ymin: 4569120 xmax: 499922 ymax: 4623140
#> Projected CRS: WGS 84 / UTM zone 17N
#>     county                       geometry
#> 1    Lucas POLYGON ((270074 4622050, 2...
#> 2   Fulton POLYGON ((220512 4622820, 2...
#> 3   Geauga POLYGON ((499922 4617790, 4...
#> 4 Williams POLYGON ((183138 4593200, 1...
#> 5 Cuyahoga POLYGON ((468455 4577270, 4...
#> 6   Ottawa POLYGON ((298927 4596820, 2...
```
