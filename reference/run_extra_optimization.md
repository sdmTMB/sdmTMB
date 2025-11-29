# Run extra optimization on an already fitted object

**\[experimental\]**

## Usage

``` r
run_extra_optimization(object, nlminb_loops = 0, newton_loops = 1)
```

## Arguments

- object:

  An object from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- nlminb_loops:

  How many extra times to run
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) optimization.
  Sometimes restarting the optimizer at the previous best values aids
  convergence.

- newton_loops:

  How many extra Newton optimization loops to try with
  [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html). Sometimes
  aids convergence.

## Value

An updated model fit of class `sdmTMB`.

## Examples

``` r
# Run extra optimization steps to help convergence:
# (Not typically needed)
fit <- sdmTMB(density ~ 0 + poly(depth, 2) + as.factor(year),
  data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie())
fit_1 <- run_extra_optimization(fit, newton_loops = 1)
#> attempting to improve convergence with Newton update(s)
#> retaining parameters from before Newton update
#> and skipping further Newton updates
max(fit$gradients)
#> [1] 5.791606e-09
max(fit_1$gradients)
#> [1] 5.791606e-09
```
