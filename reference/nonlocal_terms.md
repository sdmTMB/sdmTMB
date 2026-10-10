# Non-local covariate terms

Wrappers that mark covariates in the `nonlocal_formula` argument of
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) and
[`sdmTMB_simulate()`](https://sdmTMB.github.io/sdmTMB/reference/simulate_new.md).
They are not called directly. `diffusion(x)` spreads covariate `x`
across space via an SPDE diffusion operator with estimated scale
`kappaS_nl`. `time_lag(x)` carries `x` forward in time via \\z_t =
(x_t + \kappa_T z\_{t-1}) / (1 + \kappa_T)\\ with estimated `kappaT_nl`
and first-order autocorrelation \\\rho_T = \kappa_T / (1 + \kappa_T)\\.
When both wrap the same covariate (`~ diffusion(x) + time_lag(x)`), they
select one joint space-time operator with a single transformed predictor
and coefficient.

## Usage

``` r
diffusion(x)

time_lag(x, start = c("stationary", "zero"))
```

## Arguments

- x:

  A bare covariate name found in `data` or `nonlocal_data`.

- start:

  The transformed state before the first time slice. `"stationary"`
  (default) assumes the covariate was held at its first time slice
  beforehand, so the state starts at the operator's fixed point: the
  first slice itself for a time lag, or its spatial diffusion for the
  joint operator. Results are then invariant to shifting the covariate
  by a constant. `"zero"` starts from zero as in Thorson et al. (2026),
  which shrinks early time slices toward zero (e.g., \\z_1 = (1 -
  \rho_T) x_1\\ for a time lag alone) and makes results depend on how
  the covariate is centered: centering at a constant `c` is equivalent
  to assuming the covariate equaled `c` before the first time slice.

## Value

These functions error if called outside `nonlocal_formula`.

## See also

[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md),
[`plot_nonlocal_covariate()`](https://sdmTMB.github.io/sdmTMB/reference/nonlocal_formula_plots.md),
and the [non-local covariates
vignette](https://sdmTMB.github.io/sdmTMB/articles/nonlocal-covariates.html).

## Examples

``` r
# Pass as `nonlocal_formula` in sdmTMB(); see ?plot_nonlocal_covariate
nonlocal_formula <- ~ diffusion(x) + time_lag(x, start = "zero")
nonlocal_formula
#> ~diffusion(x) + time_lag(x, start = "zero")
#> <environment: 0x559a2f33e978>
```
