# Create an areal spatial domain for SAR/CAR models

**\[experimental\]**

Build an areal domain object with adjacency weights and unit labels. The
returned object can be supplied to `mesh` in
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) in
areal SAR/CAR workflows.

## Usage

``` r
make_areal_domain(
  spatial_domain,
  space_column = NULL,
  id_column = NULL,
  adjacency = c("rook", "queen")
)
```

## Arguments

- spatial_domain:

  A named `igraph` object or an `sf`/`sfc` polygon object.

- space_column:

  Column name in model data that identifies areal unit membership. For
  `sf` polygon input, defaults to `id_column` when `id_column` is
  supplied; otherwise defaults to `"area"`.

- id_column:

  Optional column name in an `sf` polygon `spatial_domain` containing
  areal unit IDs. If omitted for `sf` polygons, stable IDs are
  generated.

- adjacency:

  Polygon adjacency type for `sf` polygon input: `"rook"` for shared
  edges or `"queen"` for any touching boundary.

## Value

A list with class `c("sdmTMBareal", "sdmTMBdomain")`.

## See also

[`make_areal_grid()`](https://sdmTMB.github.io/sdmTMB/reference/make_areal_grid.md)

## Examples

``` r
data(ohio_df)
data(ohio_sf)

domain <- make_areal_domain(ohio_sf, id_column = "county")
domain$n_s
#> [1] 88
domain$unit_names
#>  [1] "Lucas"      "Fulton"     "Geauga"     "Williams"   "Cuyahoga"  
#>  [6] "Ottawa"     "Wood"       "Lorain"     "Sandusky"   "Trumbull"  
#> [11] "Henry"      "Erie"       "Defiance"   "Summit"     "Portage"   
#> [16] "Huron"      "Medina"     "Seneca"     "Paulding"   "Hancock"   
#> [21] "Putnam"     "Mahoning"   "Ashland"    "Crawford"   "Richland"  
#> [26] "Wayne"      "Van Wert"   "Wyandot"    "Stark"      "Columbiana"
#> [31] "Allen"      "Hardin"     "Mercer"     "Carroll"    "Morrow"    
#> [36] "Marion"     "Auglaize"   "Holmes"     "Tuscarawas" "Jefferson" 
#> [41] "Knox"       "Logan"      "Union"      "Shelby"     "Coshocton" 
#> [46] "Delaware"   "Harrison"   "Darke"      "Licking"    "Champaign" 
#> [51] "Guernsey"   "Miami"      "Belmont"    "Muskingum"  "Franklin"  
#> [56] "Madison"    "Clark"      "Noble"      "Fairfield"  "Perry"     
#> [61] "Montgomery" "Preble"     "Monroe"     "Greene"     "Pickaway"  
#> [66] "Morgan"     "Fayette"    "Hocking"    "Washington" "Warren"    
#> [71] "Butler"     "Clinton"    "Athens"     "Ross"       "Vinton"    
#> [76] "Highland"   "Hamilton"   "Clermont"   "Brown"      "Jackson"   
#> [81] "Meigs"      "Pike"       "Adams"      "Gallia"     "Scioto"    
#> [86] "Lawrence"   "Ashtabula"  "Lake"      

# \donttest{
fit <- sdmTMB(
  cases ~ pct_male,
  data = ohio_df,
  mesh = domain,
  spatial_model = "car",
  family = poisson(link = "log"),
  offset = log(ohio_df$pop)
)
# }
```
