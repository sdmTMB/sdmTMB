# Construct an SPDE mesh for sdmTMB

Construct an SPDE mesh for use with sdmTMB.

## Usage

``` r
make_mesh(
  data,
  xy_cols,
  type = c("kmeans", "cutoff", "cutoff_search"),
  cutoff,
  n_knots,
  seed = 42,
  mesh = NULL,
  fmesher_func = fmesher::fm_rcdt_2d_inla,
  convex = NULL,
  concave = convex,
  ...
)

# S3 method for class 'sdmTMBmesh'
plot(x, ...)
```

## Arguments

- data:

  A data frame.

- xy_cols:

  A character vector of x and y column names contained in `data`. These
  should likely be in an equal distance projection. For a helper
  function to convert to UTMs, see
  [`add_utm_columns()`](https://sdmTMB.github.io/sdmTMB/reference/add_utm_columns.md).

- type:

  Method to create the mesh. Also see `mesh` argument to supply your own
  mesh.

- cutoff:

  An optional cutoff if type is `"cutoff"`. The minimum allowed triangle
  edge length.

- n_knots:

  The number of desired knots if `type` is not `"cutoff"`.

- seed:

  Random seed. Affects
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) determination
  of knot locations if `type = "kmeans"`.

- mesh:

  An optional mesh created via fmesher instead of using the above
  convenience options.

- fmesher_func:

  Which fmesher function to use. Options include
  [`fmesher::fm_rcdt_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_rcdt_2d.html)
  and
  [`fmesher::fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
  along with version without the `_inla` on the end.

- convex:

  If specified, passed to
  [`fmesher::fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html).
  Distance to extend non-convex hull from data.

- concave:

  If specified, passed to
  [`fmesher::fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html).
  "Minimum allowed reentrant curvature". Defaults to `convex`.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

- x:

  Output from `make_mesh()`.

## Value

`make_mesh()`: A list of class `sdmTMBmesh`. The element `mesh` is the
output from `fmesher_func` (default is
[`fmesher::fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)).
See `mesh$mesh$n` for the number of vertices.

`plot.sdmTMBmesh()`: A plot of the mesh and data points. To make your
own ggplot2 version, pass `your_mesh$mesh` to `inlabru::gg()`.

## Examples

``` r
# Extremely simple cutoff:
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 5, type = "cutoff")
plot(mesh)


# Using a k-means algorithm to assign vertices:
mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
plot(mesh)


# \donttest{
# But, it's better to develop more tailored meshes:

# Pass arguments via '...' to fmesher::fm_mesh_2d_inla():
mesh <- make_mesh(
  pcod, c("X", "Y"),
  fmesher_func = fmesher::fm_mesh_2d_inla,
  cutoff = 8, # minimum triangle edge length
  max.edge = c(20, 40), # inner and outer max triangle lengths
  offset = c(5, 40) # inner and outer border widths
)
plot(mesh)


# Or define a mesh directly with fmesher (formerly in INLA):
inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(pcod$X, pcod$Y), # coordinates
  max.edge = c(25, 50), # max triangle edge length; inner and outer meshes
  offset = c(5, 25),  # inner and outer border widths
  cutoff = 5 # minimum triangle edge length
)
mesh <- make_mesh(pcod, c("X", "Y"), mesh = inla_mesh)
plot(mesh)

# }
```
