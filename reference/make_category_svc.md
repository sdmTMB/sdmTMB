# Set up spatially varying coefficients for category composition models

This function helps set up the data structure, formula, and mapping
needed for fitting spatially varying coefficient models with categories
(e.g., ages, length bins, species) that have both spatial and
spatiotemporal random fields. It's particularly useful for age or length
composition standardization models.

## Usage

``` r
make_category_svc(
  data,
  category_column,
  time_column,
  share_spatial_sd = TRUE,
  share_spatiotemporal_sd = TRUE
)
```

## Arguments

- data:

  Data frame containing the composition data.

- category_column:

  Character. Name of the category column (e.g., "Age", "length_bin",
  "species").

- time_column:

  Character. Name of the time column (e.g., "Year").

- share_spatial_sd:

  Logical. If `TRUE`, all categories share the same spatial SD. If
  `FALSE`, each category gets its own spatial SD.

- share_spatiotemporal_sd:

  Logical. If `TRUE`, all category-time combinations share the same
  spatiotemporal SD. If `FALSE`, each gets its own.

## Value

A list containing:

- `data_expanded`: Data frame with added model matrix columns for use in
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- `svc_formula`: Formula for the `spatial_varying` argument in
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- `svc_map`: Map list for the `map` argument in
  [`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md).

- `info`: List with summary information about the model structure.

## Details

This function creates spatially varying coefficient structures for
composition models by setting up:

1.  **Spatial fields**: One field per category (e.g., age-specific
    spatial fields) 2. **Spatiotemporal fields**: One field per
    category-time combination (e.g., age-year fields)

The sharing of variance parameters is controlled by `share_spatial_sd`
and `share_spatiotemporal_sd`. When `TRUE`, all fields of that type
share the same variance parameter, which is more parsimonious but
assumes similar variance magnitudes across categories.

The resulting model structure allows each category to have its own
spatial pattern and temporal variation while controlling parameter
sharing for identifiability and computational efficiency.

## Examples

``` r
set.seed(123)
data <- data.frame(
  age = factor(rep(1:3, each = 20)),
  year = rep(2020:2022, 20),
  abundance = rnorm(60),
  x = runif(60), y = runif(60)
)

# Set up model components
setup <- make_category_svc(
  data = data,
  category_column = "age",
  time_column = "year",
  share_spatial_sd = TRUE,
  share_spatiotemporal_sd = TRUE
)

# Check the setup
setup$info
#> $n_categories
#> [1] 3
#> 
#> $n_times
#> [1] 3
#> 
#> $n_spatial_terms
#> [1] 3
#> 
#> $n_spatiotemporal_terms
#> [1] 9
#> 
#> $n_variance_parameters
#> [1] 2
#> 
#> $spatial_terms
#> [1] "age1" "age2" "age3"
#> 
#> $spatiotemporal_terms
#> [1] "factor(year)2020:age1" "factor(year)2021:age1" "factor(year)2022:age1"
#> [4] "factor(year)2020:age2" "factor(year)2021:age2" "factor(year)2022:age2"
#> [7] "factor(year)2020:age3" "factor(year)2021:age3" "factor(year)2022:age3"
#> 
#> $share_spatial_sd
#> [1] TRUE
#> 
#> $share_spatiotemporal_sd
#> [1] TRUE
#> 

# See the age composition standardization vignette for more details
```
