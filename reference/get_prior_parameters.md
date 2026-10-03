# Inspect custom priors

`get_prior_parameters()` lists the parameters that custom priors can
refer to: the elements of `par` and `theta` passed to the `custom` and
`custom_log_jacobian` functions in
[`sdmTMBpriors()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md),
with their values. Use it with a fitted model or one set up with
`do_fit = FALSE`.

## Usage

``` r
get_prior_parameters(object, random = FALSE)

get_prior_densities(object)
```

## Arguments

- object:

  An [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
  model fit with the RTMB backend. For `get_prior_densities()`, it must
  have custom priors.

- random:

  Include the random effects (e.g., random field values)? There are
  often many of them.

## Value

`get_prior_parameters()`: a data frame with one row per element:

- `expression`: how to refer to the element in a custom prior function;

- `list`, `name`: `"par"` or `"theta"`, and the element's name;

- `label`: the coefficient name, where available;

- `value`: the estimate (or starting value if the model isn't fit);

- `status`: for `par`, `"estimated"`, `"random"`, or `"shared"`
  (estimated jointly with other elements via `map`); `"derived"` for
  `theta`.

`get_prior_densities()`: a data frame with one row per term added to the
objective: `type` (`"density"`, or `"log_jacobian"` if
`bayesian = TRUE`), `term` (the returned name, or its position), and
`log_density`.

## Details

`get_prior_densities()` evaluates the custom prior functions at the
estimated parameters, including random effects at their conditional
modes, and returns the log density of each term. Like the objective, it
evaluates `custom_log_jacobian` only if `bayesian = TRUE`.

`par` holds the raw parameters, with the names and shapes used
internally. Matrix columns are model components (e.g., the two parts of
a delta model); the second component's main effects are in `b_j2`.
`theta` holds natural-scale transformations of `par`, such as
`phi = exp(ln_phi)` and `range = sqrt(8) / kappa`.

Elements of `par` that are fixed, by `map` or because the model doesn't
use them, are left out, since a prior on them would have no effect. All
elements of `theta` are listed, including transformations of parameters
this model doesn't use (e.g., `tweedie_p` in a model without a Tweedie
family). Check that the `par` elements underlying a `theta` element are
listed before putting a prior on it. Standard deviations of random
fields the model doesn't include (e.g., `sigma_O` without a spatial
field) are zero.

These names and shapes are internal and could change in future versions
of sdmTMB. Check them with `get_prior_parameters()` when writing custom
priors.

## See also

[`sdmTMBpriors()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)

## Examples

``` r
fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011,
  family = tweedie(), mesh = pcod_mesh_2011,
  control = sdmTMBcontrol(backend = "rtmb"),
  priors = sdmTMBpriors(custom = function(par, theta) {
    c(depth = RTMB::dnorm(par$b_j[2], 0, 1, log = TRUE))
  }))
get_prior_parameters(fit)
#>                 expression  list        name        label      value    status
#> 1               par$b_j[1]   par         b_j  (Intercept)  2.8939851 estimated
#> 2               par$b_j[2]   par         b_j depth_scaled -0.6321196 estimated
#> 3          par$ln_tau_O[1]   par    ln_tau_O         <NA>  0.2689947 estimated
#> 4       par$ln_kappa[1, 1]   par    ln_kappa         <NA> -2.2958680    shared
#> 5       par$ln_kappa[2, 1]   par    ln_kappa         <NA> -2.2958680    shared
#> 6            par$thetaf[1]   par      thetaf         <NA>  0.3632222 estimated
#> 7            par$ln_phi[1]   par      ln_phi         <NA>  2.7275598 estimated
#> 8        theta$kappa[1, 1] theta       kappa         <NA>  0.1006740   derived
#> 9        theta$kappa[2, 1] theta       kappa         <NA>  0.1006740   derived
#> 10       theta$range[1, 1] theta       range         <NA> 28.0949210   derived
#> 11       theta$range[2, 1] theta       range         <NA> 28.0949210   derived
#> 12     theta$sigma_O[1, 1] theta     sigma_O         <NA>  2.1411888   derived
#> 13 theta$log_sigma_O[1, 1] theta log_sigma_O         <NA>  0.7613612   derived
#> 14     theta$sigma_E[1, 1] theta     sigma_E         <NA>  0.0000000   derived
#> 15 theta$log_sigma_E[1, 1] theta log_sigma_E         <NA>  0.0000000   derived
#> 16            theta$rho[1] theta         rho         <NA>  0.0000000   derived
#> 17     theta$sigma_V[1, 1] theta     sigma_V         <NA>  0.0000000   derived
#> 18    theta$rho_time[1, 1] theta    rho_time         <NA>  0.0000000   derived
#> 19            theta$phi[1] theta         phi         <NA> 15.2955166   derived
#> 20      theta$tweedie_p[1] theta   tweedie_p         <NA>  1.5898202   derived
#> 21      theta$p_extreme[1] theta   p_extreme         <NA>  0.5000000   derived
#> 22      theta$mix_ratio[1] theta   mix_ratio         <NA>  1.3678794   derived
get_prior_densities(fit)
#>      type  term log_density
#> 1 density depth   -1.118726
```
