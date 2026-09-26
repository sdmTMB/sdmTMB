# Transform a mesh object into a mesh with correlation barriers

Moved to the [sdmTMBextra](https://github.com/sdmTMB/sdmTMBextra)
package. Make sure to load sdmTMBextra *after* sdmTMB.

## Usage

``` r
add_barrier_mesh(
  spde_obj = deprecated(),
  barrier_sf = deprecated(),
  range_fraction = 0.2,
  proj_scaling = 1,
  plot = FALSE
)
```

## Arguments

- spde_obj:

  Output from
  [`make_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md).

- barrier_sf:

  An sf object with polygons defining the barriers. For example, a
  coastline dataset for ocean data. **Note that this object must have
  the same projection as the data used to generate the x and y columns
  in `spde_obj`.**

- range_fraction:

  The fraction of the spatial range that barrier triangles have.

- proj_scaling:

  If `spde_obj` was created with scaling of the coordinates after the
  projection (e.g., dividing UTMs by 1000 so the spatial range is on a
  reasonable scale) the x and y values in `spde_obj` are multiplied by
  this scaling factor before applying the projection from `barrier_sf`.

- plot:

  Logical.

## Value

Deprecated. See the [sdmTMBextra](https://github.com/sdmTMB/sdmTMBextra)
package.

## Examples

``` r
if (FALSE) { # \dontrun{
add_barrier_mesh()
} # }
```
