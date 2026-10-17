# sdmTMB model description

## Introduction

This vignette describes the statistical model underlying sdmTMB. sdmTMB
fits spatial and spatiotemporal generalized linear mixed models (GLMMs)
using Template Model Builder (TMB) (Kristensen *et al.* 2016), with the
model written in R through
[RTMB](https://CRAN.R-project.org/package=RTMB) and spatial random
fields approximated using the stochastic partial differential equation
(SPDE) approach (Lindgren *et al.* 2011). Further details are available
in the sdmTMB paper (Anderson *et al.* 2025), which should be referenced
and cited when using sdmTMB in a publication.

**If this vignette is being viewed on CRAN, note that many other
vignettes describing how to use sdmTMB are available on the
[documentation site](https://sdmTMB.github.io/sdmTMB/) under
[Articles](https://sdmTMB.github.io/sdmTMB/articles/).**

## Notation conventions

This vignette uses the following notation conventions, which generally
follow the guidance in Edwards & Auger-Méthé (2019):

- Greek symbols for parameters and random effects, except $`\mathbf{b}`$
  for random intercepts and slopes to match lme4 and glmmTMB,

- the Latin/Roman alphabet for data except for cases where another
  symbol is used by convention,

- bold symbols for vectors or matrices (e.g., $`\boldsymbol{\omega}`$ is
  a vector and $`\omega_{\mathbf{s}}`$ is the value of
  $`\boldsymbol{\omega}`$ at a point in space $`\mathbf{s}`$), with
  Latin/Roman symbols in non-italics (e.g., $`\mathbf{X}`$,
  $`\mathbf{Q}`$),

- $`\phi`$ for all distribution dispersion parameters for consistency
  with the code,

- $`\mathbb{E}[y]`$ to define the expected value (mean) of variable
  $`y`$,

- $`\mathrm{Var}[y]`$ to define the variance of the variable $`y`$,

- $`(\mathbf{s}, t)`$ subscripts as shorthand for observation $`i`$,
  which was sampled at location $`\mathbf{s}_i`$ and time $`t_i`$
  (several observations can share the same location and time),

- a $`^*`$ superscript represents values interpolated to observation or
  prediction locations as opposed to values at mesh vertices (e.g.,
  $`\boldsymbol{\omega}`$ vs. $`\boldsymbol{\omega}^*`$), and

- where possible, notation has been chosen to match VAST (Thorson 2019)
  to maintain consistency (e.g., $`\boldsymbol{\omega}`$ for spatial
  fields and $`\boldsymbol{\epsilon}_t`$ for spatiotemporal fields).

Tables of indices and symbols are summarized in [notation reference
tables](#notation-reference-tables) at the end of this vignette.

## sdmTMB model structure

An sdmTMB model combines relationships with measured predictors,
variation across space and time, and an observation distribution
describing variation around the expected response. For example, a model
of species density might use depth as a predictor, a spatial field for
persistent differences among locations, and a spatiotemporal field for
spatial patterns that change from year to year. With a log link, a
simple version is

``` math
\log(\mu_{\mathbf{s},t}) = \beta_0 + \beta_{\mathrm{depth}}\,\mathrm{depth}_{\mathbf{s},t} + \omega_{\mathbf{s}} + \epsilon_{\mathbf{s},t}.
```

Here, $`\mu_{\mathbf{s},t}`$ is the expected density conditional on the
random fields. The spatial field $`\omega_{\mathbf{s}}`$ describes
spatial variation remaining after accounting for depth that persists
through time, while $`\epsilon_{\mathbf{s},t}`$ describes additional
spatial deviations that change through time. All terms on the right are
on the log scale: a field value of $`\log(2)`$ doubles the expected
density, holding the other terms constant. The observation distribution
describes variation around this conditional mean. The fields describe
spatial and spatiotemporal patterns, not their causes.

The general model structure is

``` math
\begin{aligned}
y_{\mathbf{s},t} &\sim \mathcal{D} \left( \mu_{\mathbf{s},t}, \phi \right),\\
\mu_{\mathbf{s},t} &= f^{-1} \left( \eta_{\mathbf{s},t} \right),\\
\eta_{\mathbf{s},t} &=
\underbrace{\mathbf{X}^{\mathrm{main}}_{\mathbf{s},t} \boldsymbol{\beta}}_{\substack{\text{main}\\\text{effects}}} +
\underbrace{O_{\mathbf{s},t}}_{\text{offset}} +
\underbrace{\sum_{k} \mathbf{z}_{k,\mathbf{s},t}^\top \mathbf{b}_{k,g}}_{\substack{\text{random intercepts}\\\text{and slopes}}} +
\underbrace{\mathbf{X}^{\mathrm{tvc}}_{\mathbf{s},t} \boldsymbol{\gamma}_{t}}_{\substack{\text{time-varying}\\\text{coefficients}}} +
\underbrace{\mathbf{X}^{\mathrm{svc}}_{\mathbf{s},t} \boldsymbol{\zeta}_{\mathbf{s}}}_{\substack{\text{spatially varying}\\\text{coefficients}}} +
\underbrace{\omega_{\mathbf{s}}}_{\substack{\text{spatial}\\\text{field}}} +
\underbrace{\epsilon_{\mathbf{s},t}}_{\substack{\text{spatiotemporal}\\\text{field}}},
\end{aligned}
```

where

- $`y_{\mathbf{s},t}`$ represents the response data at point
  $`\mathbf{s}`$ and time $`t`$;
- $`\mathcal{D}`$ represents the observation distribution
  ([family](#observation-model-families)) with parameter
  $`\mu_{\mathbf{s},t}`$, dispersion parameter $`\phi`$ where
  applicable, and, for some families, additional parameters (e.g., the
  Tweedie power parameter);
- $`\mu = f^{-1}(\eta)`$ represents the family’s main parameter, which
  for most families is the conditional mean of the response (i.e.,
  $`\mathbb{E}[y_{\mathbf{s},t} \mid \text{random effects}] = \mu_{\mathbf{s},t}`$);
- $`f`$ represents a link function (e.g., log or logit) and $`f^{-1}`$
  represents its inverse;
- $`\eta`$ represents the linear predictor (i.e., in link space, before
  applying $`f^{-1}`$);
- $`\mathbf{X}^{\mathrm{main}}`$, $`\mathbf{X}^{\mathrm{tvc}}`$, and
  $`\mathbf{X}^{\mathrm{svc}}`$ represent design matrices (the
  superscript identifiers ‘main’ = main effects, ‘tvc’ = time-varying
  coefficients, and ‘svc’ = spatially varying coefficients);
- $`\boldsymbol{\beta}`$ represents a vector of fixed-effect
  coefficients;
- $`O_{\mathbf{s},t}`$ represents an offset: a covariate with a
  coefficient fixed at one that enters on the linear predictor scale
  (e.g., log effort with a log link);
- $`\mathbf{b}_{k,g}`$ represents the vector of random intercepts and/or
  slopes for level $`g`$ of grouping factor $`k`$ (with $`g`$ being the
  level of factor $`k`$ for this observation) and
  $`\mathbf{z}_{k,\mathbf{s},t}`$ the corresponding covariate values (a
  row of the random-effect design matrix $`\mathbf{Z}`$);
- $`\boldsymbol{\gamma}_{t}`$ represents a vector of time-varying
  coefficients, where each coefficient $`\gamma_{p,t}`$ can follow a
  random walk or AR(1) process;
- $`\boldsymbol{\zeta}_{\mathbf{s}}`$ represents the values of the
  spatially varying coefficients at location $`\mathbf{s}`$, where each
  coefficient field
  $`\boldsymbol{\zeta}_{l} \sim \mathrm{MVNormal}(\boldsymbol{0},\boldsymbol{\Sigma}_{\zeta,l})`$
  is a Gaussian random field;
- $`\omega_{\mathbf{s}}`$ represents a spatial random field,
  $`\boldsymbol{\omega}\sim \mathrm{MVNormal}(\boldsymbol{0},\boldsymbol{\Sigma}_\omega)`$,
  with $`\boldsymbol{\Sigma}_\omega`$ approximating a Matérn covariance
  via the SPDE approach (described [below](#gaussian-random-fields));
  and
- $`\epsilon_{\mathbf{s},t}`$ represents a spatiotemporal random field
  with IID, AR(1), or random walk dynamics through time; its spatial
  covariance is based on the Matérn model, and its marginal variance
  depends on the temporal process.

For binomial and beta-binomial counts, $`\mu`$ is a success probability
and the count mean is $`N\mu`$. For truncated negative binomials,
ordered beta, and mixture families, $`\mu`$ describes an underlying
distribution or component rather than the full response mean; the family
sections explain these distinctions. In the predictor equations, scalar
field terms such as $`\omega_{\mathbf{s}}`$ mean the field evaluated at
the observation location. We use stars explicitly when distinguishing
mesh-vertex values from their interpolated values (e.g.,
$`\boldsymbol{\omega}^* = \mathbf{A}\boldsymbol{\omega}`$).

Penalized smoothers are included in
$`\mathbf{X}^{\mathrm{main}}_{\mathbf{s},t} \boldsymbol{\beta}`$ for
brevity (their penalized coefficients are random effects; see
[Smoothers](#smoothers)), and [nonlocal (distributed lag)
terms](#nonlocal-formula-covariate-diffusion) are omitted from the
overview. Delta (hurdle) families have two linear predictors, each of
which can contain the above components (see [Delta
models](#delta-models)).

A single sdmTMB model will rarely, if ever, contain all of the above
components.

## Linear predictor components

This section describes each component of the linear predictor in more
detail, using ‘$`\ldots`$’ to represent the other optional components.

### Main effects

``` math
\begin{aligned}
\mu_{\mathbf{s},t} &= f^{-1} \left( \mathbf{X}^{\mathrm{main}}_{\mathbf{s},t} \boldsymbol{\beta}+ \ldots \right)
\end{aligned}
```

Within
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md),
$`\mathbf{X}^{\mathrm{main}}_{\mathbf{s},t} \boldsymbol{\beta}`$ is
defined by the `formula` argument and represents the main-effect model
matrix and a corresponding vector of coefficients. This main effect
formula can contain optional penalized smoothers or non-linear functions
as defined below.

#### Smoothers

Smoothers in sdmTMB are implemented with the same formula syntax
familiar to mgcv (Wood 2017) users fitting GAMs (generalized additive
models). Smooths are implemented in the formula using `+ s(x)`, which
implements a smooth from
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html). Within these
smooths, the same syntax commonly used in
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html) can be applied,
e.g. 2-dimensional smooths may be constructed with `+ s(x, y)`; smooths
can be specific to various factor levels, `+ s(x, by = group)`; smooths
can vary according to a continuous variable, `+ s(x, by = x2)`; the
basis function dimensions may be specified, e.g. `+ s(x, k = 4)` (see
[`?mgcv::choose.k`](https://rdrr.io/pkg/mgcv/man/choose.k.html)); and
various types of splines may be constructed such as cyclic splines to
model seasonality, e.g. `+ s(month, bs = "cc", k = 12)`.

sdmTMB supports penalized smooths, which balance fit to the data against
excessive wiggliness. The function
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
represents each supported smooth using unpenalized coefficients and
penalized coefficients treated as random effects, with associated design
matrices. This allows the smooth to be estimated in a mixed-effects
modelling framework. This is the same approach as is implemented in the
R packages gamm4 (Wood & Scheipl 2020) and brms (Bürkner 2017).

#### Linear break-point threshold models

The linear break-point or “hockey stick” model can be used to describe
threshold or asymptotic responses. This function consists of two pieces,
so that for $`x < c`$, $`s(x) = x \cdot \lambda`$, and for $`x \ge c`$,
$`s(x) = c \cdot \lambda`$. Here, $`\lambda`$ represents the slope of
the function up to the threshold (break point) $`c`$, and the product
$`c \cdot \lambda`$ represents the value at the asymptote. No
constraints are placed on $`\lambda`$ or $`c`$.

These models can be fit by including `+ breakpt(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
break-point covariate.

#### Logistic threshold models

Models with logistic threshold relationships between a predictor and the
response can be fit with the form

``` math
s(x)=\tau + \psi\ { \left[ 1+{ e }^{ -\ln \left(19\right) \cdot \left( x-s50 \right)
     / \left(s95 - s50 \right) } \right] }^{-1},
```

where $`s(x)`$ denotes the threshold function, $`\psi`$ is a scaling
parameter (controlling the height and sign of the response curve;
unconstrained), $`\tau`$ is an intercept, $`s50`$ is the value of $`x`$
at which the function has reached 50% of its total change
($`s(x) = \tau + 0.5\psi`$), and $`s95`$ is the value of $`x`$ at which
it has reached 95% of its total change ($`s(x) = \tau + 0.95\psi`$). The
parameter $`s50`$ is unconstrained but $`s95`$ is constrained to be
larger than $`s50`$.

These models can be fit by including `+ logistic(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
logistic covariate.

### Offset terms

Offset terms can be included through the `offset` argument in
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). These
are included in the linear predictor as

``` math
\begin{aligned}
  \mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + O_{\mathbf{s},t} + \ldots \right),
\end{aligned}
```

where $`O_{\mathbf{s},t}`$ is an offset term—a variable on the linear
predictor scale with its coefficient fixed at one. With a log link, the
offset is usually a **log-transformed** measure of effort (e.g.,
`log(area_swept)`), so that $`\mu_{\mathbf{s},t}`$ is proportional to
effort.

When predicting,
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
without `newdata` uses the offset from the fitted model. With `newdata`,
the offset is 0 unless an `offset` vector is supplied to
[`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md).
Therefore, if the offset represents log effort, predictions on `newdata`
are for one unit of effort (`log(1) = 0`) by default.

### Random intercepts and slopes

Multilevel (hierarchical) intercepts and slopes follow the familiar
lme4/glmmTMB formulation. Each grouping factor $`k`$ can carry an
intercept-only term or a vector of random effects (intercept plus
slopes) with a full covariance matrix:

``` math
\begin{aligned}
  \mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \sum_{k} \mathbf{z}_{k,\mathbf{s},t}^\top \mathbf{b}_{k,g} + \ldots \right),\\
  \mathbf{b}_{k,g} &\sim \mathrm{MVNormal} \left( \boldsymbol{0}, \boldsymbol{\Sigma}_{k} \right),
\end{aligned}
```

where $`g`$ is the level of grouping factor $`k`$ for the observation at
$`(\mathbf{s}, t)`$, $`\mathbf{b}_{k,g}`$ is the vector of random
effects for that level, and $`\mathbf{z}_{k,\mathbf{s},t}`$ contains the
corresponding covariate values (1 for an intercept and the covariate
value for each slope). Random effects are independent across levels and
grouping factors. The scalar special case `(1 | g)` is a random
intercept $`b_{g} \sim \mathrm{Normal}(0, \sigma_{b}^2)`$.

Use lme4-style syntax in `formula`, e.g. `(1 | g)` for random
intercepts, `(1 + x | g)` for correlated intercepts and slopes, or
`(1 | g1) + (1 + x | g2)` to give each grouping factor its own
(co)variance.

Internally, standard deviations are estimated on the log scale, and
correlations are represented by unconstrained parameters of a Cholesky
factor of the correlation matrix (TMB’s `UNSTRUCTURED_CORR`), which
guarantees a positive-definite $`\boldsymbol{\Sigma}_{k}`$ for each
grouping factor.

### Time-varying regression parameters

Parameters can be modelled as time-varying according to a random walk or
first-order autoregressive, AR(1), process. The time-series model is
defined by `time_varying_type`. For all types:

``` math
\begin{aligned}
  \mu_{\mathbf{s},t} &= f^{-1} \left( \ldots +  \mathbf{X}^{\mathrm{tvc}}_{\mathbf{s},t} \boldsymbol{\gamma}_{t} + \ldots \right),
\end{aligned}
```
where $`\boldsymbol{\gamma}_{t}`$ is an optional vector of time-varying
regression parameters and $`\mathbf{X}^{\mathrm{tvc}}_{\mathbf{s},t}`$
is the corresponding model matrix with covariate values. This is defined
via the `time_varying` argument, assuming that the `time` argument is
also supplied a column name. `time_varying` takes a *one-sided* formula.
`~ 1` implies a time-varying intercept.

The following equations describe the time-series process for a single
time-varying coefficient $`\gamma_{p,t}`$, where $`p`$ indexes the
time-varying coefficient. When multiple time-varying coefficients are
specified, each follows its time-series process independently with its
own variance parameter $`\sigma_{\gamma,p}^2`$ and (for AR1) correlation
parameter $`\rho_{\gamma,p}`$.

For `time_varying_type = 'rw'`, the initial value is unpenalized
(equivalent to an improper flat prior); only changes between successive
time steps receive a Gaussian penalty:

``` math
\begin{aligned}
  \gamma_{p,t>1} &\sim \mathrm{Normal} \left(\gamma_{p,t-1}, \sigma_{\gamma,p}^2 \right).
\end{aligned}
```

Because the initial value is estimated as part of the time-varying
formula, the same variable should not appear in the fixed-effects
formula. The formula `time_varying = ~ 1` implicitly represents a
time-varying intercept (assuming the `time` argument has been supplied)
and, in this case, the intercept should be omitted from the main-effects
formula (e.g., `formula = y ~ 0 + ...` or `formula = y ~ -1 + ...`).

For `time_varying_type = 'rw0'`, the first time step is estimated from a
mean-zero prior:

``` math
\begin{aligned}
  \gamma_{p,t=1} &\sim \mathrm{Normal} \left(0, \sigma_{\gamma,p}^2 \right),\\
  \gamma_{p,t>1} &\sim \mathrm{Normal} \left(\gamma_{p,t-1}, \sigma_{\gamma,p}^2 \right).
\end{aligned}
```
In this case, the time-varying variable (including the intercept)
*should* be included in the main effects.

For `time_varying_type = 'ar1'`:

``` math
\begin{aligned}
  \gamma_{p,t=1} &\sim \mathrm{Normal} \left(0, \sigma_{\gamma,p}^2 \right),\\
  \gamma_{p,t>1} &\sim \mathrm{Normal} \left(\rho_{\gamma,p}\gamma_{p,t-1}, (1 - \rho_{\gamma,p}^2) \sigma_{\gamma,p}^2 \right),
\end{aligned}
```
where $`\rho_{\gamma,p}`$ is the correlation between subsequent time
steps. Here, $`\sigma_{\gamma,p}`$ is the stationary marginal standard
deviation; the innovation standard deviation is
$`\sigma_{\gamma,p}\sqrt{1 - \rho_{\gamma,p}^2}`$. As with `rw0`,
include the corresponding main effect to allow a non-zero average
coefficient.

### Spatial random fields

Spatial random fields, $`\omega_{\mathbf{s}}`$, are included if
`spatial = 'on'` (or `TRUE`) and omitted if `spatial = 'off'` (or
`FALSE`).

``` math
\begin{aligned}
\mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \omega_{\mathbf{s}} + \ldots \right),\\
\boldsymbol{\omega}&\sim \mathrm{MVNormal} \left( \boldsymbol{0}, \boldsymbol{\Sigma}_\omega \right),
\end{aligned}
```

where $`\boldsymbol{\omega}`$ is the vector of field values at the mesh
vertices. The marginal standard deviation of $`\boldsymbol{\omega}`$,
$`\sigma_\omega`$, is indicated by `Spatial SD` in the printed model
output or as `sigma_O` in the output of `sdmTMB::tidy(fit, "ran_pars")`.
The ‘O’ is for ‘omega’ ($`\omega`$).

Internally, the random field is a Gaussian Markov random field (GMRF)
defined by a sparse precision (inverse covariance) matrix
$`\mathbf{Q}_\omega`$:

``` math
\boldsymbol{\omega}\sim \mathrm{MVNormal}\left(\boldsymbol{0}, \mathbf{Q}^{-1}_\omega\right).
```

$`\mathbf{Q}_\omega`$ depends on the mesh and on the parameters
$`\kappa_\omega`$ and $`\tau_\omega`$, which determine the spatial range
and $`\sigma_\omega`$ (see [Matérn
parameterization](#mat%C3%A9rn-parameterization)).

### Spatiotemporal random fields

Spatiotemporal random fields are included by default if there are
multiple time elements (`time` argument is not `NULL`) and can be set to
IID (independent and identically distributed, `'iid'`; default), AR(1)
(`'ar1'`), random walk (`'rw'`), or off (`'off'`) via the
`spatiotemporal` argument. These text values are case insensitive.

Spatiotemporal random fields are represented by
$`\boldsymbol{\epsilon}_t`$ within sdmTMB. This has been chosen to match
the representation in VAST (Thorson 2019). The standard deviation
parameter $`\sigma_\epsilon`$ is indicated by `Spatiotemporal SD` in the
printed model output or as `sigma_E` in the output of
`sdmTMB::tidy(fit, "ran_pars")`. The ‘E’ is for ‘epsilon’
($`\epsilon`$). For IID and AR(1) fields, this is the marginal standard
deviation; for a random walk, it is the innovation standard deviation,
as described below.

#### IID spatiotemporal random fields

IID spatiotemporal random fields (`spatiotemporal = 'iid'`) can be
represented as

``` math
\begin{aligned}
\mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \epsilon_{\mathbf{s},t} + \ldots \right),\\
\boldsymbol{\epsilon}_{t} &\sim \mathrm{MVNormal} \left( \boldsymbol{0}, \boldsymbol{\Sigma}_{\epsilon} \right).
\end{aligned}
```

where $`\epsilon_{\mathbf{s},t}`$ represents random field deviations at
point $`\mathbf{s}`$ and time $`t`$. The random fields are assumed
independent across time steps.

Similarly to the spatial random fields, these IID spatiotemporal random
fields are parameterized internally with a sparse precision matrix
$`\mathbf{Q}_\epsilon`$ (with parameters $`\kappa_\epsilon`$ and
$`\tau_\epsilon`$)

``` math
\boldsymbol{\epsilon}_{t} \sim \mathrm{MVNormal}\left(\boldsymbol{0}, \mathbf{Q}^{-1}_\epsilon\right).
```

#### AR(1) spatiotemporal random fields

First-order autoregressive, AR(1), spatiotemporal random fields
(`spatiotemporal = 'ar1'`) add a parameter defining the correlation
between random field deviations from one time step to the next. They are
defined as

``` math
\begin{aligned}
\mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \epsilon_{\mathbf{s},t} + \ldots \right),\\
\boldsymbol{\epsilon}_{t=1} &\sim \mathrm{MVNormal} (\boldsymbol{0}, \boldsymbol{\Sigma}_{\epsilon}),\\
\boldsymbol{\epsilon}_{t>1} &= \rho \boldsymbol{\epsilon}_{t-1} + \sqrt{1 - \rho^2} \boldsymbol{\delta}_{t},  \:
\boldsymbol{\delta}_{t} \sim \mathrm{MVNormal} \left(\boldsymbol{0}, \boldsymbol{\Sigma}_{\epsilon} \right),
\end{aligned}
```
where $`\rho`$ is the correlation between subsequent spatiotemporal
random fields and $`\boldsymbol{\delta}_t`$ are IID spatial innovations.
Because the first time step is drawn from the stationary distribution
and the $`\sqrt{1 - \rho^2}`$ factor scales the innovations, every time
step has the same marginal variance, $`\sigma_\epsilon^2`$; the process
is stationary. The correlation $`\rho`$ allows for mean-reverting
spatiotemporal fields, and is constrained to be $`-1 < \rho < 1`$.
Internally, the parameter is estimated as `ar1_phi`, which is
unconstrained. The parameter `ar1_phi` is transformed to $`\rho`$ with
$`\rho = 2 \left( \mathrm{logit}^{-1}(\texttt{ar1\_phi}) \right) - 1`$,
mapping `ar1_phi` to the range $`(-1, 1)`$.

#### Random walk spatiotemporal random fields (RW)

Random walk spatiotemporal random fields (`spatiotemporal = 'rw'`)
represent a model where the difference in spatiotemporal deviations from
one time step to the next are IID. They are defined as

``` math
\begin{aligned}
\mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \epsilon_{\mathbf{s},t} + \ldots \right),\\
\boldsymbol{\epsilon}_{t=1} &\sim \mathrm{MVNormal} (\boldsymbol{0}, \boldsymbol{\Sigma}_{\epsilon}),\\
\boldsymbol{\epsilon}_{t>1} &= \boldsymbol{\epsilon}_{t-1} +  \boldsymbol{\delta}_{t},  \:
\boldsymbol{\delta}_{t} \sim \mathrm{MVNormal} \left(\boldsymbol{0}, \boldsymbol{\Sigma}_{\epsilon} \right),
\end{aligned}
```

where the distribution of the spatiotemporal field in the initial time
step is the same as for the AR(1) model, but the absence of the $`\rho`$
parameter allows the spatiotemporal field to be non-stationary in time.
A random walk has no stationary variance: under this initialization,
$`\mathrm{Cov}(\boldsymbol{\epsilon}_t) = t\boldsymbol{\Sigma}_\epsilon`$,
so the marginal variance grows with each time step. Therefore,
$`\sigma_\epsilon`$ is the standard deviation of the innovations
$`\boldsymbol{\delta}_t`$ (and of the field in the first time step), not
the marginal standard deviation of the field in later time steps.

### Spatially varying coefficients (SVC)

Spatially varying coefficient models are defined as

``` math
\begin{aligned}
  \mu_{\mathbf{s},t} &= f^{-1} \left( \ldots + \mathbf{X}^{\mathrm{svc}}_{\mathbf{s}, t} \boldsymbol{\zeta}_{\mathbf{s}} + \ldots \right),\\
  \boldsymbol{\zeta}_{l} &\sim \mathrm{MVNormal} \left( \boldsymbol{0}, \boldsymbol{\Sigma}_{\zeta,l} \right),
\end{aligned}
```

where $`\zeta_{l,\mathbf{s}}`$ represents the value of the $`l`$-th
spatially varying coefficient field $`\boldsymbol{\zeta}_l`$ at location
$`\mathbf{s}`$. When multiple spatially varying coefficients are
specified, each is estimated as an independent spatial random field with
its own marginal variance $`\sigma_{\zeta,l}^2`$.
$`\mathbf{X}^{\mathrm{svc}}_{\mathbf{s}, t}`$ represents the
corresponding design matrix for the spatially varying coefficient terms,
as defined by a one-sided formula supplied to `spatial_varying`. For
example `spatial_varying = ~ 0 + x`, where `0` omits the intercept.
Because $`\zeta_{l,\mathbf{s}}`$ is a mean-zero deviation, the main
effect of `x` should usually also be included in `formula` so that the
coefficient at location $`\mathbf{s}`$ is
$`\beta_x + \zeta_{l,\mathbf{s}}`$.

The random fields are parameterized internally with a sparse precision
matrix $`\mathbf{Q}_{\zeta,l}`$:

``` math
\boldsymbol{\zeta}_{l} \sim \mathrm{MVNormal}\left(\boldsymbol{0}, \mathbf{Q}^{-1}_{\zeta,l}\right).
```

Each $`\mathbf{Q}_{\zeta,l}`$ has its own $`\tau_{\zeta,l}`$ (and
therefore $`\sigma_{\zeta,l}`$) but shares $`\kappa`$ with the spatial
random field $`\boldsymbol{\omega}`$.

### Nonlocal formula (covariate diffusion)

Nonlocal (distributed lag) terms let the effect of a covariate be spread
over nearby locations and/or previous time steps, rather than acting
only where and when the covariate was measured. This can be useful when
ecological responses are delayed, transported, or accumulated through
space and time. They are added with the `nonlocal_formula` argument to
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) as a
one-sided formula using
[`diffusion()`](https://sdmTMB.github.io/sdmTMB/reference/nonlocal_terms.md)
(space) and/or
[`time_lag()`](https://sdmTMB.github.io/sdmTMB/reference/nonlocal_terms.md)
(time) wrappers, e.g., `~ diffusion(x) + time_lag(x)`.
[`time_lag()`](https://sdmTMB.github.io/sdmTMB/reference/nonlocal_terms.md)
terms require a `time` argument. See (and cite) Lindmark *et al.* (2025)
for the spatial diffusion model and Thorson *et al.* (2026) for the
extension to time; those papers give the full formulation and guidance
on interpretation.

Each covariate $`x_r`$ in `nonlocal_formula` is transformed into a
smoothed covariate $`\tilde{x}_r`$, which enters the linear predictor
with an estimated coefficient:

``` math
\mu_{\mathbf{s},t} = f^{-1} \left(\ldots + \sum_{r} \beta^{\mathrm{nl}}_{r}\,\tilde{x}_{r,\mathbf{s},t} + \ldots \right).
```

If the same covariate appears in both wrappers, as in
`~ diffusion(x) + time_lag(x)`, both apply to one transformed covariate
with one coefficient. If the covariates differ, as in
`~ diffusion(x1) + time_lag(x2)`, the model has two transformed
covariates and two coefficients.

For `diffusion(x)`, the covariate at time $`t`$, $`\mathbf{x}_t`$
(evaluated at the mesh vertices), is smoothed in space by solving

``` math
\left(\mathbf{C} + \kappa_{S}^{-2}\mathbf{G}\right)\tilde{\mathbf{x}}_{t} = \mathbf{C}\mathbf{x}_{t},
```

where $`\mathbf{C}`$ and $`\mathbf{G}`$ are the SPDE mass and stiffness
matrices (see [Matérn parameterization](#mat%C3%A9rn-parameterization))
and $`\kappa_S^{-1}`$ controls the spatial scale of the diffusion. The
output reports this scale as `RMSDK` (the root-mean-squared displacement
of the diffusion kernel, $`2/\kappa_S`$, in the units of the mesh
coordinates).

For `time_lag(x)`, the covariate at each vertex is an exponentially
weighted average of current and past values:

``` math
\tilde{x}_{\mathbf{s},t} = (1 - \rho_T)\,x_{\mathbf{s},t} + \rho_T\,\tilde{x}_{\mathbf{s},t-1}, \qquad \rho_T = \frac{\kappa_T}{1 + \kappa_T} \in [0, 1),
```

By default, `time_lag(x)` uses `start = "stationary"`, assuming that the
covariate held its first observed value before the series began:
$`\tilde{x}_{\mathbf{s},0} = x_{\mathbf{s},1}`$. Use
`time_lag(x, start = "zero")` for $`\tilde{x}_{\mathbf{s},0} = 0`$. The
choice affects the transformed covariate near the start of the series.
Larger $`\rho_T`$ means a longer memory of past covariate values.

When both wrappers are applied to the same covariate, spatial smoothing
and temporal propagation act on the same evolving field (this
corresponds to `kappaST = 0` in Thorson *et al.* (2026)), and the
reported `MSDK` and `RMSDK` account for both. With both wrappers, the
stationary initial value is the spatially smoothed first covariate
slice. A separate $`\kappa_S`$ is estimated for each covariate with
spatial diffusion, and a separate $`\kappa_T`$ for each covariate with a
temporal lag.

## Gaussian random fields

This section describes how the spatial ($`\boldsymbol{\omega}`$),
spatiotemporal ($`\boldsymbol{\epsilon}_t`$), and spatially varying
coefficient ($`\boldsymbol{\zeta}_l`$) random fields are constructed.
How each field enters the model is described under [Linear predictor
components](#linear-predictor-components).

### Interpreting range and standard deviation

A spatial random field is a collection of random effects whose values
are correlated across locations. For the usual isotropic Matérn model,
correlation depends on distance, with nearby locations more strongly
correlated than distant ones. The **range** describes the distance over
which this correlation declines: at one range, correlation is
approximately 0.13, not zero. The **marginal standard deviation**
describes the typical size of field deviations at a location, on the
link scale. A larger range produces broader spatial patterns; a larger
standard deviation allows larger departures from the predictor effects.
For spatiotemporal fields, the reported standard deviation describes the
field itself under IID or AR(1) dynamics, but the changes between time
steps under a random walk.

In a model with both spatial and spatiotemporal fields,
`share_range = TRUE` (the default) estimates a common range with
separate standard deviations for the two fields. Sharing a range can be
useful when the data contain limited information about separate ranges.
Use `share_range = FALSE` to estimate separate ranges. Spatially varying
coefficient fields share the spatial range.

The following construction details are optional for understanding the
model components and interpreting these parameters.

### Matérn parameterization

For continuous-space models without barriers, sdmTMB approximates
Gaussian random fields with Matérn spatial covariance. The Matérn
covariance between the field values at locations $`\mathbf{s}_j`$ and
$`\mathbf{s}_k`$, separated by distance
$`d_{jk} = \lVert \mathbf{s}_j - \mathbf{s}_k \rVert`$, is

``` math
\mathrm{Cov}\left( \omega_{\mathbf{s}_j}, \omega_{\mathbf{s}_k} \right) = \frac{\sigma^2}{2^{\nu - 1}\Gamma(\nu)} (\kappa d_{jk})^\nu K_\nu \left( \kappa d_{jk} \right),
```

where $`\sigma^2`$ is the marginal variance, $`\nu`$ controls the
smoothness of the field, $`\Gamma`$ is the gamma function, $`K_\nu`$ is
the modified Bessel function of the second kind, and $`\kappa`$ is a
scaling parameter (larger $`\kappa`$ means correlation decays faster
with distance).

Computing with a dense Matérn covariance matrix is slow for many
locations. Instead, sdmTMB uses the stochastic partial differential
equation (SPDE) approach of Lindgren *et al.* (2011). A Gaussian field
with Matérn covariance is the solution to the SPDE
$`(\kappa^2 - \Delta)^{\alpha/2}(\tau\,\omega(\mathbf{s})) = \mathcal{W}(\mathbf{s})`$,
where $`\Delta`$ is the Laplacian, $`\mathcal{W}`$ is Gaussian white
noise, and $`\alpha = \nu + d/2`$ for $`d`$ spatial dimensions. sdmTMB
uses $`d = 2`$ and $`\alpha = 2`$, so $`\nu = 1`$. Approximating the
solution on a triangulated mesh with piecewise-linear basis functions
(the finite element method) gives field values at the mesh vertices with
a sparse precision matrix

``` math
\mathbf{Q}= \tau^2 \left( \kappa^4 \mathbf{C} + 2 \kappa^2 \mathbf{G} + \mathbf{G} \mathbf{C}^{-1} \mathbf{G} \right),
```

where $`\mathbf{C}`$ is a diagonal approximation to the finite-element
mass matrix (called mass lumping), and $`\mathbf{G}`$ is a sparse
stiffness matrix. Both matrices depend only on the mesh. The sparsity of
$`\mathbf{Q}`$ (most entries are zero) is what makes the computations
fast.

The parameters $`\kappa`$ and $`\tau`$ determine the reported range and
marginal standard deviation:

``` math
\textrm{range} = \frac{\sqrt{8}}{\kappa}, \qquad
\sigma = \frac{1}{\sqrt{4 \pi}\, \tau \kappa}.
```

These relationships describe the stationary Matérn field; the finite
mesh and its boundaries introduce approximation error. For a fixed
$`\tau`$, increasing $`\kappa`$ decreases $`\sigma`$. Sharing a range
means sharing $`\kappa`$, while separate $`\tau`$ parameters allow
different field standard deviations. Throughout, $`\mathbf{Q}_\omega`$,
$`\mathbf{Q}_\epsilon`$, and $`\mathbf{Q}_{\zeta,l}`$ denote full
precision matrices, including the $`\tau^2`$ scaling (the code applies
this scaling separately when evaluating the field densities).

### Projection $`\mathbf{A}`$ matrix

In the SPDE approach, the random field at any location is defined by the
piecewise-linear basis expansion
$`\omega(\mathbf{s}) = \sum_k \psi_k(\mathbf{s})\, \omega_k`$, where
$`\omega_k`$ is the field value at mesh vertex $`k`$ and $`\psi_k`$ is a
“tent” function that is 1 at vertex $`k`$ and decreases linearly to 0 at
the neighbouring vertices. Evaluating the field at observed or
prediction locations is therefore a multiplication by a sparse
projection matrix $`\mathbf{A}`$ with elements
$`A_{ik} = \psi_k(\mathbf{s}_i)`$(Lindgren & Rue 2015):

``` math
\boldsymbol{\omega}^* = \mathbf{A}\boldsymbol{\omega},
```
where $`\boldsymbol{\omega}^*`$ represents the values of the spatial
random fields at the observed locations or predicted data locations. The
matrix $`\mathbf{A}`$ has a row for each data point or prediction point
and a column for each mesh vertex. Each row has at most three non-zero
elements: the barycentric weights of the vertices of the triangle
containing that location. This interpolation is part of the model
definition rather than a separate approximation step. The same
interpolation happens for any spatiotemporal random fields

``` math
\boldsymbol{\epsilon}_t^* = \mathbf{A}\boldsymbol{\epsilon}_t.
```

### Anisotropy

TMB allows for anisotropy, where spatial covariance may depend on
direction ([full
details](https://kaskr.github.io/adcomp/namespaceR__inla.html)).
Anisotropy can be turned on or off with the logical `anisotropy`
argument to
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). There
are a number of ways to implement anisotropic covariance (Fuglstad *et
al.* 2015), and we adopt geometric anisotropy defined by a 2-parameter
matrix $`\mathbf{H}`$, which replaces Euclidean distance with
$`d_{jk} = \sqrt{(\mathbf{s}_j - \mathbf{s}_k)^\top \mathbf{H}^{-1} (\mathbf{s}_j - \mathbf{s}_k)}`$.
The elements of $`\mathbf{H}`$ are defined by the parameters $`h_1 > 0`$
and $`h_2`$ so that $`H_{1,1} = h_{1}`$, $`H_{1,2} = H_{2,1} = h_{2}`$
and $`H_{2,2} = (1 + h_{2}^2) / h_{1}`$. $`\mathbf{H}`$ is symmetric and
positive definite with determinant 1, so it stretches distances along
one axis and shrinks them along the perpendicular axis without changing
area. The result is a different range in each direction, with the
orientation of the major axis estimated.

Once a model is fitted with
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md), the
anisotropy relationships may be plotted using the
[`plot_anisotropy()`](https://sdmTMB.github.io/sdmTMB/reference/plot_anisotropy.md)
function, which takes the fitted object as an argument. If a barrier
mesh is used, anisotropy is disabled.

### Incorporating physical barriers into the SPDE

In some cases the spatial domain of interest may be complex and bounded
by some barrier such as by land or water (e.g., coastlines, islands,
lakes). SPDE models allow for physical barriers to be incorporated into
the modelling (Bakka *et al.* 2019). With
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
models, the mesh construction occurs in two steps: the user (1)
constructs a mesh with a call to
[`sdmTMB::make_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md),
and (2) passes the mesh to
[`sdmTMBextra::add_barrier_mesh()`](https://rdrr.io/pkg/sdmTMBextra/man/add_barrier_mesh.html).
The barriers must be constructed as `sf` objects (Pebesma 2018) with
polygons defining the barriers. See
[`?sdmTMBextra::add_barrier_mesh`](https://rdrr.io/pkg/sdmTMBextra/man/add_barrier_mesh.html)
for an example.

The barrier implementation requires the user to select a fraction value
(`range_fraction` argument) that defines the fraction of the usual
spatial range when crossing the barrier (Bakka *et al.* 2019). For
example, if the range was estimated at 10 km, `range_fraction = 0.2`
would assign a 2 km range to the triangles marked as barrier. That makes
correlation decay much faster through the barrier than around it, which
effectively discourages spatial smoothing from “jumping” across
landmasses. From experimentation, values around 0.1 or 0.2 seem to work
well but values much lower than 0.1 can result in convergence issues.
The `range_fraction` is fixed, not estimated.

Following Bakka *et al.* (2019), the barrier model lets the range vary
in space through the SPDE

``` math
\omega(\mathbf{s}) - \nabla \cdot \left( \frac{r(\mathbf{s})^2}{8} \nabla \omega(\mathbf{s}) \right) = r(\mathbf{s}) \sqrt{\frac{\pi}{2}}\, \sigma\, \mathcal{W}(\mathbf{s}),
```

where $`\nabla`$ is the gradient, $`r(\mathbf{s}) = r`$ in the normal
part of the domain and $`r(\mathbf{s}) = r \times`$`range_fraction` in
barrier triangles. With a constant range $`r = \sqrt{8}/\kappa`$, this
is the Matérn SPDE given above, with marginal standard deviation
$`\sigma`$. With barriers, the field is no longer a stationary Matérn
field, so its precision matrix is not the $`\mathbf{Q}`$ given above.
sdmTMB builds it with the formulation used in the
[INLAspacetime](https://CRAN.R-project.org/package=INLAspacetime)
package, parameterized directly by the range and marginal standard
deviation. The reported range is still $`\sqrt{8}/\kappa`$, and
$`\sigma`$ is still computed from $`\kappa`$ and $`\tau`$ as above. Both
describe the field away from barriers; near barriers, the variance and
correlation are modified.

[This website](https://haakonbakkagit.github.io/btopic128.html) by
Francesco Serafini and Haakon Bakka provides an illustration with INLA.
The original TMB implementation in sdmTMB was adapted from code written
by Olav Nikolai Breivik and Hans Skaug at the [TMB Case
Studies](https://github.com/skaug/tmb-case-studies) GitHub site.

### Areal autoregressive models (SAR/CAR)

In addition to continuous-space SPDE random fields, sdmTMB can fit areal
spatial autoregressive models with `spatial_model = "sar"` or
`spatial_model = "car"`. These require an areal domain supplied to
`mesh` with the help of
[`make_areal_domain()`](https://sdmTMB.github.io/sdmTMB/reference/make_areal_domain.md).

For both areal options, the latent spatial and spatiotemporal effects
are Gaussian Markov random fields with precision defined by the areal
weight matrix $`\mathbf{W}`$. For the spatial field and IID or AR(1)
spatiotemporal fields,

``` math
\begin{aligned}
\boldsymbol{\omega}&\sim \mathrm{MVNormal}\left(\boldsymbol{0}, \sigma_\omega^2 \mathbf{Q}^{-1}\right),\\
\boldsymbol{\epsilon}_t &\sim \mathrm{MVNormal}\left(\boldsymbol{0}, \sigma_\epsilon^2 \mathbf{Q}^{-1}\right).
\end{aligned}
```

For a random walk, $`\sigma_\epsilon^2 \mathbf{Q}^{-1}`$ instead
describes the covariance of the initial field and subsequent
innovations. Unlike the SPDE case, $`\sigma_\omega`$ and
$`\sigma_\epsilon`$ here are scale parameters, not marginal standard
deviations: the diagonal entries of $`\sigma^2\mathbf{Q}^{-1}`$ depend
on $`\mathbf{W}`$ and the autocorrelation parameter and generally differ
among areas.

For SAR (`spatial_model = "sar"`), the precision matrix is

``` math
\mathbf{Q}_{\mathrm{SAR}} = \left(\mathbf{I} - \rho_{\mathrm{SAR}}\mathbf{W}\right)^\top
\left(\mathbf{I} - \rho_{\mathrm{SAR}}\mathbf{W}\right),
```

where $`\rho_{\mathrm{SAR}} \in (-1, 1)`$ is reported as `rho_sar`. By
default, SAR uses a row-normalized $`\mathbf{W}`$
(`sar_weight_style = "row"`), with an option to use raw weights. With
`sar_weight_style = "raw"`, SAR uses the unnormalized adjacency/weight
matrix directly, so neighbour influence depends on the original edge
weights (and, if unweighted, on the number of neighbours) rather than
being scaled to sum to 1 within each row.

For CAR (`spatial_model = "car"`), the precision matrix is

``` math
\mathbf{Q}_{\mathrm{CAR}} = \mathbf{D} - \alpha_{\mathrm{CAR}}\mathbf{W},
```

where $`\mathbf{D}`$ is diagonal with elements
$`D_{ii} = \sum_j W_{ij}`$ (the number of neighbours for unweighted
adjacency; set to 1 for areas with no neighbours) and
$`\alpha_{\mathrm{CAR}} \in (0, 1)`$ is reported as `alpha_car`. CAR
requires a symmetric adjacency matrix.

## Observation model families

Here we describe the main observation families that are available in
sdmTMB and comment on their parametrization, statistical properties,
utility, and code representation in sdmTMB. Families are grouped by
outcome type to make it easier to locate an appropriate observation
model.

### Bounded or binary outcomes (0–1, proportions)

#### Binomial

``` math
\operatorname{Binomial} \left(N, \mu \right)
```
where $`N`$ is the size or number of trials, and $`\mu`$ is the
probability of success for each trial. If $`N = 1`$, the distribution
becomes the Bernoulli distribution. Internally, the distribution is
parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaecb5a18095a320b42e2d20c4b120f5f5)
in TMB, which is numerically stable when probabilities approach 0 or 1.
Following the structure of
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), lme4, and glmmTMB, a
binomial family can be specified in one of 4 ways:

1.  the response may be a factor (and the model classifies the first
    level versus all others)
2.  the response may be binary (0/1)
3.  the response can be a matrix of form `cbind(success, failure)`, or
4.  the response may be the observed proportions, and the `weights`
    argument is used to specify the Binomial size ($`N`$) parameter
    (`probability ~ ..., weights = N`).

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaee11f805f02bc1febc6d7bf0487671be).

Example: `family = binomial(link = "logit")`

#### Beta-binomial

``` math
\operatorname{BetaBinomial} \left(N, \mu, \phi \right)
```
where $`N`$ is the number of trials, $`\mu`$ is the mean success
probability, and $`\phi`$ is a precision parameter that controls
overdispersion relative to a Binomial distribution. The implied Beta
parameters are $`\alpha = \mu \phi`$ and $`\beta = (1 - \mu)\phi`$, and
the variance is

``` math
\mathrm{Var}[y] = N \mu (1 - \mu)\, \frac{\phi + N}{\phi + 1},
```

For $`N > 1`$ and $`0 < \mu < 1`$, this exceeds the Binomial variance
for finite $`\phi`$ and approaches it as $`\phi \to \infty`$. For
$`N = 1`$, the distribution reduces to Bernoulli. Available links are
logit and cloglog. This family is useful for overdispersed counts of
successes/failures (e.g., aggregated Bernoulli data, proportions with
extra-binomial variation).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L505-L512).

Example: `family = betabinomial(link = "logit")`

#### Beta

``` math
\operatorname{Beta} \left(\mu \phi, (1 - \mu) \phi \right)
```
where $`\mu`$ is the mean and $`\phi`$ is a precision parameter. This
parametrization follows Ferrari & Cribari-Neto (2004) and the betareg R
package (Cribari-Neto & Zeileis 2010). The variance is
$`\mu (1 - \mu) / (\phi + 1)`$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga5324c83759d5211c7c2fbbad37fa8e59).

Example: `family = Beta(link = "logit")`

#### Ordered beta

The ordered beta distribution (Kubinec 2023) is for continuous
proportions on the closed interval $`[0, 1]`$ that include exact zeros
and/or ones. A single linear predictor $`\eta`$ and two estimated
cutpoints $`c_1 < c_2`$ define
$`\Pr(y = 0) = 1 - \mathrm{logit}^{-1}(\eta - c_1)`$ and
$`\Pr(y = 1) = \mathrm{logit}^{-1}(\eta - c_2)`$, and values in
$`(0, 1)`$ follow $`\operatorname{Beta}(\mu \phi, (1 - \mu)\phi)`$ with
$`\mu = \mathrm{logit}^{-1}(\eta)`$. Here, $`\mu`$ is the mean
conditional on $`0 < y < 1`$; the full response mean is
$`\Pr(y = 1) + \Pr(0 < y < 1)\mu`$. It is a parsimonious alternative to
zero-one-inflated beta models because all three components share one
linear predictor.

Example: `family = ordbeta(link = "logit")`

### Count data

#### Poisson

``` math
\operatorname{Poisson} \left( \mu \right)
```
where $`\mu`$ represents the mean and $`\mathrm{Var}[y] = \mu`$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa1ed15503e1441a381102a8c4c9baaf1).

Example: `family = poisson(link = "log")`

#### Censored Poisson

A Poisson distribution in which some observations are censored: the true
count is only known to lie in the interval $`[y, U]`$, where $`U`$ may
be infinite (right censoring). The likelihood of a censored observation
is $`\Pr(y \le Y \le U)`$ under $`\operatorname{Poisson}(\mu)`$. Upper
bounds are supplied with `sdmTMB(censored_upper = ...)`, where `Inf`
indicates no upper bound and a value equal to $`y`$ indicates an
uncensored observation. Watson *et al.* (2023) developed this approach
to account for hook competition in longline surveys, where observed
catch counts are lower bounds on what would have been caught in the
absence of competition for baited hooks. See the [hook competition
article](https://sdmTMB.github.io/sdmTMB/articles/hook-competition.html)
for a worked example.

Example: `family = censored_poisson(link = "log")`

#### Negative Binomial 2 (NB2)

``` math
\operatorname{NB2} \left( \mu, \phi \right)
```

where $`\mu`$ is the mean and $`\phi`$ is the dispersion parameter. The
variance scales quadratically with the mean
$`\mathrm{Var}[y] = \mu + \mu^2 / \phi`$(Hilbe 2011). The NB2
parametrization is more commonly seen in ecology than the NB1.
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga76266c19046e04b651fce93aa0810351).

Example: `family = nbinom2(link = "log")`

#### Negative Binomial 1 (NB1)

``` math
\operatorname{NB1} \left( \mu, \phi \right)
```

where $`\mu`$ is the mean and $`\phi`$ is the dispersion parameter. The
variance scales linearly with the mean
$`\mathrm{Var}[y] = \mu + \mu \phi = \mu (1 + \phi)`$(Hilbe 2011).
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/sdmTMB.cpp#L577-L582)
based on NB2 and borrowed from glmmTMB.

Example: `family = nbinom1(link = "log")`

#### Truncated negative binomial

Zero-truncated versions of the NB2 and NB1 distributions, for positive
counts ($`y \ge 1`$). Here, $`\mu`$ is the mean of the untruncated
distribution; the mean of the truncated distribution is
$`\mu / \Pr(Y > 0)`$. These are mainly used as the positive component of
delta models for counts (e.g.,
[`delta_truncated_nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)).

Example: `family = truncated_nbinom2(link = "log")` or
`family = truncated_nbinom1(link = "log")`

#### Negative binomial 2 mixture

This is a 2 component mixture that extends the NB2 distribution,
following mixture-distribution approaches for aggregation (e.g.,
schooling) in survey data (Thorson *et al.* 2011).

``` math
(1 - p) \cdot \operatorname{NB2} \left( \mu_1, \phi \right) + p \cdot \operatorname{NB2} \left( \mu_2, \phi \right)
```

where $`\mu_1 = f^{-1}(\eta)`$ is the mean of the smaller component,
$`\mu_2 = r\mu_1`$ is the mean of the larger component for an estimated
ratio $`r > 1`$, and $`p`$ is the probability of the larger component.
The full response mean is $`(1-p)\mu_1 + p\mu_2`$.

Example: `family = nbinom2_mix(link = "log")`

### Positive continuous outcomes

#### Gamma

``` math
\operatorname{Gamma} \left( \phi, \frac{\mu}{\phi}  \right)
```
where $`\phi`$ represents the Gamma shape and $`\mu / \phi`$ represents
the scale. The mean is $`\mu`$ and variance is $`\mu^2 / \phi`$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gab0e2205710a698ad6a0ed39e0652c9a3).

Example: `family = Gamma(link = "log")`

#### Lognormal

sdmTMB uses the “bias-corrected” lognormal distribution where $`\phi`$
represents the standard deviation in log-space:

``` math
\operatorname{Lognormal} \left( \log \mu - \frac{\phi^2}{2}, \phi^2 \right).
```
Because of the bias correction, $`\mathbb{E}[y] = \mu`$ and
$`\mathrm{Var}[\log y] = \phi^2`$.

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L47-L54)
based on the TMB [`dnorm()`](https://rdrr.io/r/stats/Normal.html) normal
density.

Example: `family = lognormal(link = "log")`

#### Generalized gamma

``` math
\operatorname{GenGamma} \left( \mu, \phi, Q^{\mathrm{gg}}\right)
```

sdmTMB implements the Prentice (1974) parameterization introduced for
spatiotemporal models and index standardization by Dunic *et al.*
(2025), with parameters mean $`\mu`$, scale $`\phi`$, and a shape
parameter $`Q^{\mathrm{gg}}`$ (reported as `Generalized gamma Q`).

Here, $`\mu`$ is the mean on the data scale, $`\phi`$ is a scale
parameter that equals the log-scale standard deviation in the lognormal
limit, and $`Q^{\mathrm{gg}}`$ controls the shape
($`Q^{\mathrm{gg}}\to 0`$ yields the lognormal;
$`Q^{\mathrm{gg}}= \phi`$ yields the gamma). See Dunic *et al.* (2025)
for the full PDF in the Prentice formulation as implemented. Links
available: identity, log, inverse. This flexibility is useful for
right-skewed positive responses with tails heavier or lighter than
gamma/lognormal, and is often paired in a hurdle/mixture as
[`delta_gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
for zero-inflated biomass or catch data (see Dunic *et al.* (2025)).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L7-L43).

Example: `family = gengamma(link = "log")`; delta/hurdle:
`family = delta_gengamma(link1 = "logit", link2 = "log")`.

#### Gamma mixture

This is a 2 component mixture that extends the Gamma distribution,
motivated by mixture-distribution treatments of aggregation in survey
data (Thorson *et al.* 2011),

``` math
(1 - p) \cdot \operatorname{Gamma} \left( \phi, \frac{\mu_{1}}{\phi}  \right) + p \cdot \operatorname{Gamma} \left( \phi, \frac{\mu_{2}}{\phi}  \right),
```
where $`\phi`$ represents the Gamma shape, $`\mu_{1} / \phi`$ represents
the scale for the first (smaller component) of the distribution,
$`\mu_{2} / \phi`$ represents the scale for the second (larger
component) of the distribution, and $`p`$ controls the contribution of
each component to the mixture (also interpreted as the probability of
larger events).

As with the [NB2 mixture](#negative-binomial-2-mixture),
$`\mu_1 = f^{-1}(\eta)`$ and $`\mu_2 = r\mu_1`$ for an estimated ratio
$`r > 1`$. The full response mean is
$`(1-p) \cdot \mu_{1} + p \cdot \mu_{2}`$. The variance follows the
usual mixture formula:
$`(1-p)\left(\mu_1^2 / \phi + \mu_1^2\right) + p\left(\mu_2^2 / \phi + \mu_2^2\right) - \left[(1-p)\mu_1 + p\mu_2\right]^2`$.

Here, and for the other mixture distributions, the probability of the
larger mean can be obtained from
`plogis(fit$model$par[["logit_p_extreme"]])` and the ratio of the larger
mean to the smaller mean can be obtained from
`1 + exp(fit$model$par[["log_ratio_mix"]])`. The standard errors are
available in the TMB sdreport: `fit$sd_report`.

If you wish to fix the probability of a large (i.e., extreme) mean,
which can be hard to estimate, you can fix this value and pass this to
the family:

``` r

sdmTMB(...,
  family = gamma_mix(link = "log", p_extreme = 0.01)
)
```

See also `family = delta_gamma_mix()` for an extension incorporating
this distribution with delta models.

#### Lognormal mixture

This is a 2 component mixture that extends the lognormal distribution,
again in the spirit of mixture approaches for aggregating/schooling data
(Thorson *et al.* 2011),

``` math
(1 - p) \cdot \operatorname{Lognormal} \left( \log \mu_{1} - \frac{\phi^2}{2}, \phi^2 \right) + p \cdot \operatorname{Lognormal} \left( \log \mu_{2} - \frac{\phi^2}{2}, \phi^2 \right).
```

As with the other mixtures, $`\mu_2 = r\mu_1`$ with
$`\mu_1 = f^{-1}(\eta)`$, and the full response mean is
$`\mathbb{E}[y] = (1-p) \cdot \mu_{1} + p \cdot \mu_{2}`$ (the
$`-\phi^2/2`$ terms make each component’s mean $`\mu_i`$). The log-scale
variance of the mixture is not simply $`\phi^2`$; it can be obtained
with the standard mixture-variance formula using component log-means
$`\log \mu_i - \phi^2/2`$ and log-variance $`\phi^2`$.

As with the Gamma mixture, $`p`$ controls the contribution of each
component to the mixture (also interpreted as the probability of larger
events).

Example: `family = lognormal_mix(link = "log")`. See also
`family = delta_lognormal_mix()` for an extension incorporating this
distribution with delta models. Like with the gamma mixture, fixed
probabilities of extreme events ($`p`$ in notation above) can be passed
in, e.g.

``` r

sdmTMB(...,
  family = delta_lognormal_mix(p_extreme = 0.01)
)
```

### Non-negative continuous outcomes with exact zeros

#### Tweedie

``` math
\operatorname{Tweedie} \left(\mu, p, \phi \right), \: 1 < p < 2
```

where $`\mu`$ is the mean, $`p`$ is the power parameter constrained
between 1 and 2, and $`\phi`$ is the dispersion parameter. The Tweedie
distribution (a compound Poisson-gamma distribution for $`1 < p < 2`$)
can be helpful for modelling data that are non-negative and continuous
with exact zeros. The variance is $`\phi \mu^p`$. [Delta
models](#delta-models) are an alternative that model zeros and positive
values with separate linear predictors.

Internally, $`p = \mathrm{logit}^{-1} (\texttt{thetaf}) + 1`$, which
constrains it between 1 and 2, and `thetaf` is estimated as an
unconstrained parameter.

The [source
code](https://kaskr.github.io/adcomp/tweedie_8cpp_source.html) is
implemented as in the [cplm](https://CRAN.R-project.org/package=cplm)
package (Zhang 2013) and is based on Dunn & Smyth (2005). The TMB
version is defined
[here](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga262f3c2d1cf36f322a62d902a608aae0).

Example: `family = tweedie(link = "log")`

### Continuous real-valued outcomes

#### Gaussian

``` math
\operatorname{Normal} \left( \mu, \phi^2 \right)
```
where $`\mu`$ is the mean and $`\phi`$ is the standard deviation. The
variance is $`\phi^2`$.

Example: `family = gaussian(link = "identity")`

Code defined [within
TMB](https://kaskr.github.io/adcomp/dnorm_8hpp.html).

#### Student-t

``` math
\operatorname{Student-t} \left( \mu, \phi, \nu \right)
```

where $`\mu`$ is the location (mean for $`\nu > 1`$), $`\phi`$ is a
scale parameter (not the standard deviation; the variance is
$`\phi^2 \nu / (\nu - 2)`$ for $`\nu > 2`$), and $`\nu`$, the degrees of
freedom (`df`), is estimated by default and can optionally be fixed by
the user. Lower values of $`\nu`$ result in heavier tails compared to
the Gaussian distribution. Above approximately `df = 20`, the
distribution becomes very similar to the Gaussian. The Student-t
distribution with a low degrees of freedom (e.g., $`\nu \le 7`$) can be
helpful for modelling data that would otherwise be suitable for Gaussian
but needs an approach that is robust to outliers (e.g., Anderson *et
al.* 2017).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L37-L45)
based on the [`dt()`](https://rdrr.io/r/stats/TDist.html) distribution
in TMB.

Example: `family = student(link = "identity", df = 7)`

## Delta models

sdmTMB allows for several different kinds of delta (also known as
hurdle) models. These families are implemented by specifying the family
as a delta distribution. For example:

``` r

sdmTMB(
  ...,
  family = delta_gamma()
)
```

The list of supported families is included in the documentation on
[additional
families](https://sdmTMB.github.io/sdmTMB/reference/families.html),
[delta
models](https://sdmTMB.github.io/sdmTMB/articles/delta-models.html), and
[Poisson-link delta
models](https://sdmTMB.github.io/sdmTMB/articles/poisson-link.html). By
default, the `delta_*` families don’t use the Poisson link, but this
structure can be specified with `delta_gamma(type = "poisson-link")`,
which follows the formulation in Thorson (2018).

In the “standard” delta model implementation, sdmTMB constructs two
internal models (formally, two “linear predictors”), with the first
model representing presence-absence and the second model representing
the positive component (such as catch rates in fisheries applications).
The positive component can be continuous (e.g.,
[`delta_gamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md),
[`delta_lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md),
[`delta_gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)),
a proportion in $`(0, 1)`$
([`delta_beta()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)),
or a positive count from a zero-truncated negative binomial
([`delta_truncated_nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.md),
[`delta_truncated_nbinom1()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)).
The expected value is the product of the probability of a non-zero
observation and the expected value of the positive component. Default
links associated with each family can be inspected with
[`delta_lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
and equivalent functions. For the standard delta models, the default
first linear predictor link is logit. For the Poisson-link type, the
first linear predictor link is log.

The `formula`, `spatial`, `spatiotemporal`, and `share_range` arguments
of [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) can
be specified independently as a 2-element list. For example, the spatial
random field might be estimated for only the first linear predictor
with:

``` r

sdmTMB(
  ...,
  family = delta_gamma(),
  spatial = list("on", "off")
)
```

Or we may want separate main-effects formulas. For example:

``` r

sdmTMB(
  formula = list(
    y ~ depth + I(depth^2),
    y ~ 1),
  family = delta_gamma()
)
```

All other arguments are shared between the linear predictors.

Currently if delta models contain smoothers, both components must have
the same main-effects formula.

## Optimization

### Optimization details

By default, sdmTMB fits models by maximum marginal likelihood: the
likelihood is integrated over the random effects, accounting for their
uncertainty when estimating the remaining parameters. These remaining
parameters include regression coefficients, variance and correlation
parameters, spatial ranges, and observation-distribution parameters. We
call them the *outer* parameters below because they are optimized in an
outer loop, while the random effects are optimized in an inner loop for
each candidate set of outer-parameter values. TMB and glmmTMB call these
‘fixed effects’, which includes the variance parameters, not just the
regression coefficients. With `reml = TRUE`, the regression coefficients
are also integrated out, together with the random effects. The joint
likelihood is written with
[RTMB](https://CRAN.R-project.org/package=RTMB) (the default backend; a
C++ TMB template can be used instead with
`sdmTMBcontrol(backend = "tmb")`), and Template Model Builder (TMB)
(Kristensen *et al.* 2016) uses automatic differentiation and the
Laplace approximation to integrate over the random effects. This yields
an approximation to the marginal log likelihood of the outer parameters
and its gradient, and the negative marginal log likelihood is minimized
with the non-linear optimization routine
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) in R (Gay 1990;
R Core Team 2021). For given outer-parameter values, the random effects
are set to the values that maximize the joint log likelihood (their
conditional modes), around which the Laplace approximation is taken
(Kristensen *et al.* 2016).

Like AD Model Builder (Fournier *et al.* 2012), sdmTMB fits parameters
in phases by default (`multiphase = TRUE` in
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)).
Phased estimation is usually faster and more stable, but not always,
because it requires building the model an extra time. In sdmTMB, the
first phase holds random effects at their initial values and estimates
the regression and observation-distribution parameters to obtain
starting values. The second phase fits the full model, using the
first-phase estimates as starting values and estimating the random
effects and their variance and correlation parameters as well.

In some cases, a single call to
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) may not result
in convergence (e.g., the maximum gradient of the marginal likelihood
with respect to the outer parameters is not small enough yet), and the
algorithm may need to be run multiple times. In the
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
function, we include an argument `nlminb_loops` that will restart the
optimization at the previous best values. The number of `nlminb_loops`
should generally be small (e.g., 2 or 3), and defaults to 1. After
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), sdmTMB takes
Newton steps to further reduce the gradient: it computes the Hessian
$`\mathcal{H}`$ of the negative marginal log likelihood numerically with
[`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) and updates
the outer parameters $`\boldsymbol{\theta}`$ as
$`\boldsymbol{\theta} - \mathcal{H}^{-1} \mathbf{g}`$, where
$`\mathbf{g}`$ is the gradient. The number of Newton steps is set with
`newton_loops` in
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
(default 1). A step is only accepted if it does not increase the
negative marginal log likelihood. If a model is already fit, the
function
[`sdmTMB::run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.md)
can run additional optimization loops with either routine to further
reduce the maximum gradient.

### Assessing convergence

The [`sanity()`](https://sdmTMB.github.io/sdmTMB/reference/sanity.md)
function runs a set of basic convergence checks on a fitted model (e.g.,
the Hessian is positive definite, gradients are small, standard errors
are not `NA` or very large, and random field variances have not
collapsed to zero) and is a good first step. Much of the guidance around
diagnostics and glmmTMB also applies to sdmTMB, e.g. [the glmmTMB
vignette on
troubleshooting](https://CRAN.R-project.org/package=glmmTMB).
Optimization with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) involves
specifying the number of iterations and evaluations (`eval.max` and
`iter.max`) and the tolerances (`abs.tol`, `rel.tol`, `x.tol`,
`xf.tol`)—a greater number of iterations and smaller tolerance
thresholds increase the chance that the optimal solution is found, but
more evaluations translates into longer computation time. Warnings of
non-positive-definite Hessian matrices (accompanied by parameters with
`NA`s for standard errors) often mean models are improperly specified
given the data. Standard errors can be observed in the output of
`print.sdmTMB()` or by checking `fit$sd_report`. The maximum gradient of
the marginal likelihood with respect to the outer parameters can be
checked by inspecting `fit$gradients`. Guidance varies, but gradients
should be very close to zero (e.g., on the order of $`10^{-3}`$ or
smaller after reasonable parameter scaling) before assuming the fitting
routine is consistent with convergence. If maximum gradients are already
relatively small, they can sometimes be reduced further with additional
optimization calls beginning at the previous best parameter vector as
described above with
[`sdmTMB::run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.md).

## Notation reference tables

Tables of all major indices (Table 1) and symbols (Table 2) are provided
here for quick reference.

| Symbol         | Description                                      |
|:---------------|:-------------------------------------------------|
| $`\mathbf{s}`$ | Index for space; a vector of x and y coordinates |
| $`t`$          | Index for time                                   |
| $`g`$          | Level of a random-effect grouping factor         |
| $`k`$          | Random-effect grouping factor                    |
| $`p`$          | Index for time-varying coefficient               |
| $`l`$          | Index for spatially varying coefficient (SVC)    |

Table 1: Subscript notation {.table}

| Symbol | Code | Description |
|:---|:---|:---|
| $`y`$ | `y_i` | Observed response data |
| $`\mu`$ | `mu` | Inverse-linked family parameter (usually the conditional response mean) |
| $`\eta`$ | `eta_i` | Linear predictor before applying the inverse link ($`f^{-1}`$) |
| $`\phi`$ | `phi` | A dispersion parameter for a distribution (estimated as `ln_phi`) |
| $`f`$ | `fit$family$link` | Link function |
| $`f^{-1}`$ | `fit$family$linkinv` | Inverse link function |
| $`\boldsymbol{\beta}`$ | `b_j` | Fixed-effect coefficient vector |
| $`\mathbf{X}`$ | `X_ij` | A predictor model matrix |
| $`\mathbf{Z}`$ | `Zt_list` | Random effect design matrix; $`\mathbf{z}_{k,\mathbf{s},t}`$ is the row for one observation and grouping factor $`k`$ |
| $`O_{\mathbf{s}, t}`$ | `offset_i` | An offset variable at point $`\mathbf{s}`$ and time $`t`$ |
| $`\omega_{\mathbf{s}}`$ | `omega_s` | Spatial random field at point $`\mathbf{s}`$ (vertex) |
| $`\omega_{\mathbf{s}}^*`$ | `omega_s_A` | Spatial random field at point $`\mathbf{s}`$ (interpolated) |
| $`\zeta_{\mathbf{s}}`$ | `zeta_s` | Spatially varying coefficient random field at point $`\mathbf{s}`$ (vertex) |
| $`\zeta_{\mathbf{s}}^*`$ | `zeta_s_A` | Spatially varying coefficient random field at point $`\mathbf{s}`$ (interpolated) |
| $`\epsilon_{\mathbf{s}, t}`$ | `epsilon_st` | Spatiotemporal random field at point $`\mathbf{s}`$ and time $`t`$ (vertex) |
| $`\epsilon_{\mathbf{s}, t}^*`$ | `epsilon_st_A_vec` | Spatiotemporal random field at point $`\mathbf{s}`$ and time $`t`$ (interpolated) |
| $`\delta_{\mathbf{s},t}`$ | \- | AR(1) or random walk spatiotemporal innovations (vertex) |
| $`\gamma_{p,t}`$ | `b_rw_t` | Time-varying coefficients |
| $`\mathbf{b}_{k,g}`$ | `re_b_pars` | Vector of random intercepts/slopes for level $`g`$ of grouping factor $`k`$ |
| $`\boldsymbol{\Sigma}_\omega`$ | \- | Spatial random field covariance matrix |
| $`\boldsymbol{\Sigma}_\zeta`$ | \- | Spatially varying coefficient random field covariance matrix |
| $`\boldsymbol{\Sigma}_\epsilon`$ | \- | Spatiotemporal marginal covariance (IID/AR1) or innovation covariance (RW) |
| $`\mathbf{Q}_\omega`$ | `Q` | Spatial random field precision matrix (code omits $`\tau^2`$) |
| $`\mathbf{Q}_\zeta`$ | `Q` | Spatially varying coefficient precision matrix; shares $`\kappa_\omega`$, with its own $`\tau`$ (code omits $`\tau^2`$) |
| $`\mathbf{Q}_\epsilon`$ | `Q_st` | Spatial precision matrix of the spatiotemporal process; `Q` if the range is shared (code omits $`\tau^2`$) |
| $`\boldsymbol{\Sigma}_{k}`$ | `re_cov_pars` | Random effect covariance matrix for grouping factor $`k`$ (log SDs and correlation parameters) |
| $`\sigma_\omega`$ | `sigma_O` | Spatial marginal SD for SPDE fields; scale parameter for SAR/CAR |
| $`\sigma_\epsilon`$ | `sigma_E` | Spatiotemporal marginal SD (IID/AR1) or innovation SD (RW); scale parameter for SAR/CAR |
| $`\sigma_\zeta`$ | `sigma_Z` | Spatially varying coefficient marginal SD for SPDE fields; scale parameter for SAR/CAR |
| $`\sigma_\gamma`$ | `sigma_V` | Time-varying coefficient innovation SD (RW/RW0) or stationary marginal SD (AR1) |
| $`\tau`$ | `ln_tau_O`, `ln_tau_E`, `ln_tau_Z` | SPDE precision scaling parameters (log scale) |
| $`\kappa_\omega`$ | `kappa[1, ]` | Spatial Matérn scaling parameter (estimated as `ln_kappa`) |
| $`\kappa_\epsilon`$ | `kappa[2, ]` | Spatiotemporal Matérn scaling parameter |
| $`\sqrt{8}/\kappa`$ | `range` | Distance at which correlation drops to approximately 0.13 |
| $`\rho`$ | `rho` | Correlation between spatiotemporal random fields in subsequent time steps (estimated as `ar1_phi`) |
| $`\rho_{\gamma}`$ | `rho_time` | Correlation between time-varying coefficients in subsequent time steps |
| $`\rho_{\mathrm{SAR}}`$ | `rho_sar` | SAR spatial autocorrelation parameter |
| $`\alpha_{\mathrm{CAR}}`$ | `alpha_car` | CAR spatial autocorrelation parameter |
| $`\lambda`$, $`c`$ | `s_slope`, `s_cut` | Breakpoint slope and break point |
| $`s50`$, $`s95`$, $`\psi`$ | `s50`, `s95`, `s_max` | Logistic threshold parameters |
| $`\kappa_S`$, $`\kappa_T`$ | `kappaS_nl`, `kappaT_nl` | Nonlocal (distributed lag) spatial and temporal scale parameters |
| $`p`$ (Tweedie) | `tweedie_p` | Tweedie power parameter (estimated as `thetaf`) |
| $`\nu`$ | `student_df` | Student-t degrees of freedom |
| $`Q^{\mathrm{gg}}`$ | `gengamma_Q` | Generalized gamma shape parameter ($`Q \rightarrow 0`$ gives lognormal; $`Q=\phi`$ gives gamma) |
| $`c_1`$, $`c_2`$ | `psi` | Ordered beta cutpoints |
| $`p`$ (mixture) | `p_extreme` | Probability of the larger mixture component |
| $`\mathbf{A}`$ | `A_station` | Sparse projection matrix to interpolate from mesh vertices to the unique data or prediction locations |
| $`\mathbf{H}`$ | `H` | 2-parameter matrix defining geometric anisotropy (estimated as `ln_H_input`) |

Table 2: Symbol notation, names in the model output or RTMB model code,
and descriptions. {.table}

## References

Anderson, S.C., Branch, T.A., Cooper, A.B. & Dulvy, N.K. (2017).
Black-swan events in animal populations. *Proceedings of the National
Academy of Sciences*, **114**, 3252–3257.

Anderson, S.C., Ward, E.J., English, P.A., Barnett, L.A.K. & Thorson,
J.T. (2025). [sdmTMB: An R package for fast, flexible, and user-friendly
generalized linear mixed effects models with spatial and spatiotemporal
random fields](https://doi.org/10.18637/jss.v115.i02). *Journal of
Statistical Software*, **115**, 1–46.

Bakka, H., Vanhatalo, J., Illian, J., Simpson, D. & Rue, H. (2019).
Non-stationary Gaussian models with physical barriers. *arXiv:1608.03787
\[stat\]*. Retrieved from <https://arxiv.org/abs/1608.03787>

Bürkner, P.-C. (2017). [brms: An R package for Bayesian multilevel
models using Stan](https://doi.org/10.18637/jss.v080.i01). *Journal of
Statistical Software*, **80**, 1–28.

Cribari-Neto, F. & Zeileis, A. (2010). Beta regression in R. *Journal of
Statistical Software*, **34**, 1–24.

Dunic, J.C., Conner, J., Anderson, S.C. & Thorson, J.T. (2025). [The
generalized gamma is a flexible distribution that outperforms
alternatives when modelling catch rate
data](https://doi.org/10.1093/icesjms/fsaf040). *ICES Journal of Marine
Science*, **82**, fsaf040.

Dunn, P.K. & Smyth, G.K. (2005). [Series evaluation of Tweedie
exponential dispersion model
densities](https://doi.org/10.1007/s11222-005-4070-y). *Statistics and
Computing*, **15**, 267–280.

Edwards, A.M. & Auger-Méthé, M. (2019). Some guidance on using
mathematical notation in ecology. *Methods in Ecology and Evolution*,
**10**, 92–99.

Ferrari, S. & Cribari-Neto, F. (2004). [Beta regression for modelling
rates and proportions](https://doi.org/10.1080/0266476042000214501).
*Journal of Applied Statistics*, **31**, 799–815.

Fournier, D.A., Skaug, H.J., Ancheta, J., Ianelli, J., Magnusson, A.,
Maunder, M.N., Nielsen, A. & Sibert, J. (2012). AD model builder: Using
automatic differentiation for statistical inference of highly
parameterized complex nonlinear models. *Optimization Methods and
Software*, **27**, 233–249.

Fuglstad, G.-A., Lindgren, F., Simpson, D. & Rue, H. (2015). Exploring a
new class of non-stationary spatial Gaussian random fields with varying
local anisotropy. *Statistica Sinica*, **25**, 115–133.

Gay, D.M. (1990). Usage summary for selected optimization routines.
*Computing Science Technical Report*, **153**, 1–21.

Hilbe, J.M. (2011). *Negative Binomial Regression*. Cambridge University
Press.

Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H. & Bell, B.M. (2016).
[TMB: Automatic differentiation and Laplace
approximation](https://doi.org/10.18637/jss.v070.i05). *Journal of
Statistical Software*, **70**, 1–21.

Kubinec, R. (2023). [Ordered beta regression: A parsimonious,
well-fitting model for continuous data with lower and upper
bounds](https://doi.org/10.1017/pan.2022.20). *Political Analysis*,
**31**, 519–536.

Lindgren, F. & Rue, H. (2015). [Bayesian spatial modelling with
R-INLA](https://doi.org/10.18637/jss.v063.i19). *Journal of Statistical
Software*, **63**, 1–25.

Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between
Gaussian fields and Gaussian Markov random fields: The stochastic
partial differential equation approach. *J. R. Stat. Soc. Ser. B Stat.
Methodol.*, **73**, 423–498.

Lindmark, M., Anderson, S.C. & Thorson, J.T. (2025). [Estimating
scale-dependent covariate responses using two-dimensional diffusion
derived from the stochastic partial differential equation
method](https://doi.org/10.1111/2041-210X.70177). *Methods in Ecology
and Evolution*, **17**, 207–218.

Pebesma, E. (2018). Simple Features for R: Standardized Support for
Spatial Vector Data. *The R Journal*, **10**, 439–446. Retrieved from
<https://doi.org/10.32614/RJ-2018-009>

Prentice, R.L. (1974). [A log gamma model and its maximum likelihood
estimation](https://doi.org/10.1093/biomet/61.3.539). *Biometrika*,
**61**, 539–544.

R Core Team. (2021). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing, Vienna, Austria.
Retrieved from <https://www.R-project.org/>

Thorson, J.T. (2019). [Guidance for decisions using the Vector
Autoregressive Spatio-Temporal (VAST) package in stock, ecosystem,
habitat and climate
assessments](https://doi.org/10.1016/j.fishres.2018.10.013). *Fisheries
Research*, **210**, 143–161.

Thorson, J.T. (2018). [Three problems with the conventional delta-model
for biomass sampling data, and a computationally efficient
alternative](https://doi.org/10.1139/cjfas-2017-0266). *Canadian Journal
of Fisheries and Aquatic Sciences*, **75**, 1369–1382.

Thorson, J., Anderson, S. & Lindmark, M. (2026). [Temperature carryover
effect revealed for marine fishes using spatio-temporal distributed lag
models](https://doi.org/10.32942/X2W95P). *EcoEvoRxiv*.

Thorson, J.T., Stewart, I.J. & Punt, A.E. (2011). [Accounting for fish
shoals in single- and multi-species survey data using mixture
distribution models](https://doi.org/10.1139/f2011-086). *Canadian
Journal of Fisheries and Aquatic Sciences*, **68**, 1681–1693.

Watson, J., Edwards, A.M. & Auger-Méthé, M. (2023). [A statistical
censoring approach accounts for hook competition in abundance indices
from longline surveys](https://doi.org/10.1139/cjfas-2022-0159).
*Canadian Journal of Fisheries and Aquatic Sciences*, **80**, 468–486.

Wood, S.N. (2017). *Generalized additive models: An introduction with
R*, 2nd ednn. Chapman and Hall/CRC.

Wood, S. & Scheipl, F. (2020). *gamm4: Generalized Additive Mixed Models
using ’mgcv’ and ’lme4’*. Retrieved from
<https://CRAN.R-project.org/package=gamm4>

Zhang, Y. (2013). [Likelihood-based and Bayesian methods for Tweedie
compound Poisson linear mixed
models](https://doi.org/10.1007/s11222-012-9343-7). *Statistics and
Computing*, **23**, 743–757.
