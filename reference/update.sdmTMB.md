# Update an sdmTMB model

This method updates an sdmTMB model with new arguments while
automatically handling the mesh object to avoid environment issues when
loading models from saved files.

## Usage

``` r
# S3 method for class 'sdmTMB'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  An sdmTMB model object.

- formula.:

  Optional updated formula.

- ...:

  Other arguments to update in the model call.

- evaluate:

  If `TRUE` (default), the updated call is evaluated; if `FALSE`, the
  call is returned unevaluated.

## Value

An updated sdmTMB model object (if `evaluate = TRUE`) or an unevaluated
call (if `evaluate = FALSE`).

## Details

Unless `control` is supplied, the updated model reuses the fitted
model's
[`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
settings, including its backend, regardless of the current
`sdmTMB.backend` option. A supplied `control` replaces these settings,
including the backend (e.g.,
`update(fit, control = sdmTMBcontrol(backend = "rtmb"))`).

## Examples

``` r
mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
  family = tweedie(link = "log"))
fit2 <- update(fit, family = delta_gamma())
```
