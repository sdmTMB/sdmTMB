# Shrink a fitted model before saving it

Removes objects from a fitted model that sdmTMB rebuilds as needed, so
the model takes less space when saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). The stripped model
works the same as one read back from disk: its objective function is
taped again the first time a method needs it.

## Usage

``` r
strip_sdmTMB(object, strip_joint_precision = FALSE)
```

## Arguments

- object:

  A model fit with
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- strip_joint_precision:

  Also remove the joint precision matrix of the fixed and random
  effects? This usually gives the largest saving, but makes some methods
  slower. See Details.

## Value

The fitted model with the removed objects.

## Details

This removes the handle to the in-memory RTMB tape, which doesn't
survive a save and reload and can be the largest part of a saved model,
and the cached sparse Cholesky factor of the random-effects Hessian,
which sdmTMB rebuilds cheaply when needed.

The joint precision matrix is often the largest part of what remains,
but removing it (`strip_joint_precision = TRUE`) has a cost later:

- Simulation-based methods (e.g., `predict(..., nsim = )`,
  [`get_index_sims()`](https://sdmTMB.github.io/sdmTMB/reference/get_index_sims.md),
  [`spread_sims()`](https://sdmTMB.github.io/sdmTMB/reference/gather_sims.md),
  [`gather_sims()`](https://sdmTMB.github.io/sdmTMB/reference/gather_sims.md),
  and
  [`project()`](https://sdmTMB.github.io/sdmTMB/reference/project.md))
  recompute it with `RTMB::sdreport()` each time they are called. It
  isn't stored back in the model, so every call costs about as much as
  the `sdreport()` step of the original fit, which can take minutes for
  large spatiotemporal models.

- [`get_index()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md),
  [`get_cog()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md),
  and
  [`get_eao()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  can't use a shortcut based on the matrix and instead build a new
  objective function and run `RTMB::sdreport()` on it.

Other methods (e.g., [`predict()`](https://rdrr.io/r/stats/predict.html)
without `nsim`,
[`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html)) are unaffected,
as are models fit without the matrix (see `get_joint_precision` in
[`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)).

The original `object` is not changed.

## Examples

``` r
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
  family = tweedie())
fit_small <- strip_sdmTMB(fit)
object.size(serialize(fit, NULL))
#> 4271816 bytes
object.size(serialize(fit_small, NULL))
#> 1099664 bytes
p <- predict(fit_small)
```
