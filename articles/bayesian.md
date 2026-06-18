# Bayesian estimation with sdmTMB

**If the code in this vignette has not been evaluated, a rendered
version is available on the [documentation
site](https://sdmTMB.github.io/sdmTMB/index.html) under ‘Articles’.**

``` r

library(ggplot2)
library(dplyr)
library(sdmTMB)
library(rstan) # for plot() method
options(mc.cores = parallel::detectCores()) # use rstan parallel processing
```

Bayesian estimation is possible with sdmTMB by passing fitted models to
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
(Monnahan & Kristensen 2018). All sampling is then done using Stan (Stan
Development Team 2021), and output is returned as a `stanfit` object.

Why might you want to pass an sdmTMB model to Stan?

- to obtain probabilistic inference on parameters
- to avoid the Laplace approximation on the random effects
- to robustly quantify uncertainty on derived quantities not already
  calculated in the model
- in some cases, models that struggle to converge with maximum
  likelihood can be adequately sampled with MCMC given carefully chosen
  priors (e.g., Monnahan *et al.* 2021)

## Simulating data

Here we will demonstrate using a simulated dataset.

``` r

set.seed(123)
predictor_dat <- data.frame(
  X = runif(500), Y = runif(500),
  a1 = rnorm(500)
)
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
# plot(mesh)
# mesh$mesh$n
sim_dat <- sdmTMB_simulate(
  formula = ~a1,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(),
  range = 0.3,
  phi = 0.2,
  sigma_O = 0.2,
  seed = 123,
  B = c(0.8, -0.4) # B0 = intercept, B1 = a1 slope
)
```

Visualize our simulated data:

``` r

ggplot(sim_dat, aes(X, Y, colour = observed)) +
  geom_point() +
  scale_color_viridis_c()
```

## Fitting the model with marginal likelihood

First, fit a spatial random field GLMM with maximum likelihood:

``` r

fit <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh,
  family = gaussian(),
  spatial = "on"
)
fit
```

## Adding priors

In that first model fit we did not use any priors (with maximum
likelihood estimation, these also can be thought of as penalties on the
likelihood). For this first model, the priors are implied as uniform on
the internal parameter space. However, sdmTMB provides the option of
applying additional priors. Here we will show an example of applying a
Normal(0, 5) (mean, SD) prior on the intercept and a Normal(0, 1) prior
on the slope parameter. We could guess at the model matrix structure
based on our formula, but we can verify it by looking at the internal
model matrix from the previous fit (using `do_fit = FALSE` would save
time if you didn’t want to fit it the first time).

``` r

head(fit$tmb_data$X_ij[[1]])
```

Each column corresponds to the order of the `b` priors:

``` r

fit <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh,
  family = gaussian(),
  spatial = "on",
  priors = sdmTMBpriors(
    # location = vector of means; scale = vector of standard deviations:
    b = normal(location = c(0, 0), scale = c(5, 2)),
  )
)
fit
```

## Fixing a spatial correlation parameter to improve convergence

Sometimes some of the spatial correlation parameters can be challenging
to estimate with Stan. One option is to apply penalized complexity (PC)
priors with
[`sdmTMBpriors()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
to the Matérn parameters. Another option, which can also be used in
conjunction with the priors, is to fix one or more parameters at their
maximum likelihood estimate (MLE) values. Frequently, fixing the
parameter `ln_kappa` can help convergence (e.g., Monnahan *et al.*
2021). This estimated parameter is transformed into the range estimate,
so it controls the rate of spatial correlation decay.

Now we will rebuild the fitted object with fixed (‘mapped’) `ln_kappa`
parameters using the [`update()`](https://rdrr.io/r/stats/update.html)
function. We’ll use `do_fit = FALSE` to avoid actually fitting the
updated model since it’s not necessary.

``` r

# grab the internal parameter list at estimated values:
pars <- sdmTMB::get_pars(fit)
# create a 'map' vector for TMB
# factor NA values cause TMB to fix or map the parameter at the starting value:
kappa_map <- factor(rep(NA, length(pars$ln_kappa)))

# rebuild model updating some elements:
fit_mle <- update(
  fit,
  control = sdmTMBcontrol(
    start = list(
      ln_kappa = pars$ln_kappa #<
    ),
    map = list(
      ln_kappa = kappa_map #<
    )
  ),
  do_fit = FALSE #<
)
```

## Jacobian adjustments

Adding priors / penalties and fixing spatial parameters represent
strategies to help successful convergence for maximum likelihood
estimation. If we want to do true Bayesian sampling we need to make one
more adjustment to our function call: accounting for non-linear
transformations of parameters with Jacobian adjustments.

What are Jacobian adjustments and why do we need them? Jacobian
adjustments are necessary when the parameters of a model are transformed
in a way that changes their scale or distribution. A good example of
this is the estimation of variance parameters. Whether we’re interested
in spatial, spatiotemporal, or residual variation, the quantity of
interest is usually the variance or standard deviation $`\sigma`$. These
quantities are constrained to be greater than 0, so a widely used
estimation strategy is to estimate them in log space, which is not
constrained. With ‘ln_sigma’ estimated, ‘sigma = exp(ln_sigma)’ can be
calculated internal to a model and used to calculate the likelihood.
There are a number of helpful references detailing the math behind this
in greater detail including the [Stan
manual](https://mc-stan.org/docs/stan-users-guide/reparameterization.html#changes-of-variables).
Without equations, the Jacobian adjustment can be thought of as properly
stretching the posterior distribution of parameters to account for the
transformation of variables.

In `sdmTMB`, we can turn these Jacobian adjustments on with the flag
`bayesian = TRUE`. Applying this to our `fit_mle` object,

``` r

fit_bayes <- update(fit_mle,
  bayesian = TRUE
)
```

It is important to emphasize that this `bayesian` flag needs to be
enabled to any model passed to Stan; MCMC estimation without it will
lead to biased parameter estimates.

## Passing the model to tmbstan

Now we can pass the `$tmb_obj` element of our model to
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html). We
are only using 1000 iterations and 2 chains so this vignette builds
quickly. In practice, you will likely want to use more (e.g., 2000
iterations, 4 chains).

``` r

fit_stan <- tmbstan::tmbstan(
  fit_bayes$tmb_obj,
  iter = 1000, chains = 2,
  seed = 8217 # ensures repeatability
)
```

Sometimes you may need to adjust the sampler settings such as:

``` r

tmbstan::tmbstan(
  ...,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)
```

See the Details section in
[`?rstan::stan`](https://mc-stan.org/rstan/reference/stan.html).

You can also ‘thin’ samples via the `thin` argument if working with
model predictions becomes cumbersome given a large number of required
samples.

We can look at the model:

``` r

fit_stan
```

The `Rhat` values look reasonable (\< 1.05). The `n_eff` (number of
effective samples) values mostly look reasonable (\> 100) for inference
about the mean for all parameters except the intercept (`b_j[1]`).
Furthermore, we can see correlation in the MCMC samples for `b_j[1]`. We
could try running for more iterations and chains and/or placing priors
on this and other parameters as described below (highly recommended).

Now we can use various functions to visualize the posterior:

``` r

plot(fit_stan)
pars_plot <- c("b_j[1]", "b_j[2]", "ln_tau_O", "omega_s[1]")

bayesplot::mcmc_trace(fit_stan, pars = pars_plot)
bayesplot::mcmc_pairs(fit_stan, pars = pars_plot)
```

## Posterior predictive checks

We can perform posterior predictive checks to assess whether our model
can generate predictive data that are consistent with the observations.
For this, we can make use of
[`simulate.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/simulate.sdmTMB.md)
while passing in our Stan model.
[`simulate.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/simulate.sdmTMB.md)
will take draws from the joint parameter posterior and add observation
error. We need to ensure `nsim` is less than or equal to the total
number of post-warmup samples.

``` r

set.seed(19292)
samps <- sdmTMBextra::extract_mcmc(fit_stan)
s <- simulate(fit_mle, mcmc_samples = samps, nsim = 50)
bayesplot::pp_check(
  sim_dat$observed,
  yrep = t(s),
  fun = bayesplot::ppc_dens_overlay
)
```

See `?bayesplot::pp_check`. The solid line represents the density of the
observed data and the light blue lines represent the density of 50
posterior predictive simulations. In this case, the simulated data seem
consistent with the observed data.

## Plotting predictions

We can make predictions with our Bayesian model by supplying the
posterior samples to the `mcmc_samples` argument in
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md).

``` r

pred <- predict(fit_mle, mcmc_samples = samps)
```

The output is a matrix where each row corresponds to a row of predicted
data and each column corresponds to a sample.

``` r

dim(pred)
```

We can summarize these draws in various ways to visualize them:

``` r

sim_dat$post_mean <- apply(pred, 1, mean)
sim_dat$post_sd <- apply(pred, 1, sd)

ggplot(sim_dat, aes(X, Y, colour = post_mean)) +
  geom_point() +
  scale_color_viridis_c()

ggplot(sim_dat, aes(X, Y, colour = post_sd)) +
  geom_point() +
  scale_color_viridis_c()
```

Or predict on a grid for a given value of `a1`:

``` r

nd <- expand.grid(
  X = seq(0, 1, length.out = 70),
  Y = seq(0, 1, length.out = 70),
  a1 = 0
)
pred <- predict(fit_mle, newdata = nd, mcmc_samples = samps)

nd$post_mean <- apply(pred, 1, mean)
nd$post_sd <- apply(pred, 1, sd)

ggplot(nd, aes(X, Y, fill = post_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed()

ggplot(nd, aes(X, Y, fill = post_sd)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed()
```

## Extracting parameter posterior samples

We can extract posterior samples with
[`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html),

``` r

post <- rstan::extract(fit_stan)
```

The result is a list where each element corresponds to a parameter or
set of parameters:

``` r

names(post)
hist(post$b_j[, 1])
```

As an example of calculating a derived parameter, here we will calculate
the marginal spatial random field standard deviation:

``` r

ln_kappa <- get_pars(fit_mle)$ln_kappa[1] # 2 elements since 2nd would be for spatiotemporal
ln_tau_O <- post$ln_tau_O
sigma_O <- 1 / sqrt(4 * pi * exp(2 * ln_tau_O + 2 * ln_kappa))
hist(sigma_O)
```

## Extracting the posterior of other predicted elements

By default
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
returns the overall prediction in link space when a tmbstan model is
passed in. If instead we want some other element that we might find in
the usual data frame returned by
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
when applied to a regular sdmTMB model, we can specify that through the
`sims_var` argument.

For example, let’s extract the spatial random field values `"omega_s"`.
Other options are documented in `?predict.sdmTMB()`.

``` r

fit_pred <- predict(
  fit_mle,
  newdata = nd,
  mcmc_samples = samps,
  sims_var = "omega_s" #<
)

nd$spatial_rf_mean <- apply(fit_pred, 1, mean)
nd$spatial_rf_sd <- apply(fit_pred, 1, sd)

ggplot(nd, aes(X, Y, fill = spatial_rf_mean)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed()

ggplot(nd, aes(X, Y, fill = spatial_rf_sd)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed()
```

## References

Monnahan, C. & Kristensen, K. (2018). [No-U-turn sampling for fast
bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R
packages](https://doi.org/10.1371/journal.pone.0197954). *PloS one*,
**13**.

Monnahan, C.C., Thorson, J.T., Kotwicki, S., Lauffenburger, N., Ianelli,
J.N. & Punt, A.E. (2021). [Incorporating vertical distribution in index
standardization accounts for spatiotemporal availability to acoustic and
bottom trawl gear for semi-pelagic
species](https://doi.org/10.1093/icesjms/fsab085) (O.R. Godo, Ed.).
*ICES Journal of Marine Science*, **78**, 1826–1839.

Stan Development Team. (2021). RStan: The R interface to Stan. Retrieved
from <https://mc-stan.org/>
