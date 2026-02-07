# Example fish survey data

Various fish survey datasets.

## Usage

``` r
pcod

pcod_2011

pcod_mesh_2011

qcs_grid

dogfish

yelloweye

hbll_s_grid

wcvi_grid
```

## Format

`pcod`: Trawl survey data for Pacific Cod in Queen Charlotte Sound. A
data frame.

`pcod_2011`: A version of `pcod` for years 2011 and after (smaller for
speed). A data frame.

`pcod_mesh_2011`: A mesh pre-built for `pcod_2011` for examples. A list
of class `sdmTMBmesh`.

`qcs_grid` A 2x2km prediction grid for Queen Charlotte Sound. A data
frame.

`dogfish`: Trawl survey data for Pacific Spiny Dogfish on West Coast
Vancouver Island. A data frame.

`yelloweye`: Survey data for Yelloweye Rockfish from the Hard Bottom
Longline Survey (South) off West Coast Vancouver Island.

`hbll_s_grid`: A survey domain grid to go with `yelloweye`. A data
frame.

`wcvi_grid`: A survey domain grid to go with `dogfish`. A data frame.
