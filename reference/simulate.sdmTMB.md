# Simulate from a fitted sdmTMB model

`simulate.sdmTMB` is an S3 method for producing a matrix of simulations
from a fitted model. This is similar to
[`lme4::simulate.merMod()`](https://rdrr.io/pkg/lme4/man/simulate.merMod.html)
and
[`glmmTMB::simulate.glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/simulate.glmmTMB.html).
It can be used with the DHARMa package among other uses.

## Usage

``` r
# S3 method for class 'sdmTMB'
simulate(
  object,
  nsim = 1L,
  seed = sample.int(1e+06, 1L),
  type = c("mle-eb", "mle-mvn"),
  model = c(NA, 1, 2),
  newdata = NULL,
  re_form = NULL,
  mle_mvn_samples = c("single", "multiple"),
  mcmc_samples = NULL,
  return_tmb_report = FALSE,
  observation_error = TRUE,
  size = NULL,
  silent = FALSE,
  ...
)
```

## Arguments

- object:

  sdmTMB model

- nsim:

  Number of response lists to simulate. Defaults to 1.

- seed:

  Random number seed

- type:

  How parameters should be treated. `"mle-eb"`: fixed effects are at
  their maximum likelihood (MLE) estimates and random effects are at
  their empirical Bayes (EB) estimates. `"mle-mvn"`: fixed effects are
  at their MLEs but random effects are taken from a single approximate
  sample. This latter option is a suggested approach if these
  simulations will be used for goodness of fit testing (e.g., with the
  DHARMa package).

- model:

  If a delta/hurdle model, which model to simulate from? `NA` =
  combined, `1` = first model, `2` = second mdoel.

- newdata:

  Optional new data frame from which to simulate.

- re_form:

  `NULL` to specify a simulation conditional on fitted random effects
  (this only simulates observation error). `~0` or `NA` to simulate new
  random affects (smoothers, which internally are random effects, will
  not be simulated as new).

- mle_mvn_samples:

  Applies if `type = "mle-mvn"`. If `"single"`, take a single MVN draw
  from the random effects. If `"multiple"`, take an MVN draw from the
  random effects for each of the `nsim`.

- mcmc_samples:

  An optional matrix of MCMC samples. See
  [`extract_mcmc()`](https://sdmTMB.github.io/sdmTMB/reference/extract_mcmc.md)
  in the [sdmTMBextra](https://github.com/sdmTMB/sdmTMBextra) package.

- return_tmb_report:

  Return the TMB report from
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)? This lets you
  parse out whatever elements you want from the simulation. Not usually
  needed.

- observation_error:

  Logical. Simulate observation error?

- size:

  A vector of size (trials) in the case of a binomial family with
  `newdata` specified. If left `NULL` with `newdata`, will be assumed to
  be a vector of 1s.

- silent:

  Logical. Silent?

- ...:

  Extra arguments passed to
  [`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md).
  E.g., one may wish to pass an `offset` argument if `newdata` are
  supplied in a model with an offset.

## Value

Returns a matrix; number of columns is `nsim`.

## See also

[`sdmTMB_simulate()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_simulate.md)

## Examples

``` r
# start with some data simulated from scratch:
set.seed(1)
predictor_dat <- data.frame(X = runif(300), Y = runif(300), a1 = rnorm(300))
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
dat <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(),
  range = 0.5,
  sigma_O = 0.2,
  seed = 42,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)
fit <- sdmTMB(observed ~ 1 + a1, data = dat, family = poisson(), mesh = mesh)

# simulate from the model:
s1 <- simulate(fit, nsim = 300)
dim(s1)
#> [1] 300 300

# test whether fitted models are consistent with the observed number of zeros:
sum(s1 == 0)/length(s1)
#> [1] 0.3297667
sum(dat$observed == 0) / length(dat$observed)
#> [1] 0.3466667

# simulate with random effects sampled from their approximate posterior
s2 <- simulate(fit, nsim = 1, params = "mle-mvn")
# these may be useful in conjunction with DHARMa simulation-based residuals

# simulate with new random fields:
s3 <- simulate(fit, nsim = 1, re_form = ~ 0)
```
