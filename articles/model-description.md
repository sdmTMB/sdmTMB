# sdmTMB model description

## Introduction

This vignette describes the statistical model underlying sdmTMB. Further
details are available in the sdmTMB paper (Anderson *et al.* 2024),
which should be referenced and cited when using sdmTMB in a publication.

## Other vignettes available

**If this vignette is being viewed on CRAN, note that many other
vignettes describing how to use sdmTMB are available on the
[documentation site](https://sdmTMB.github.io/sdmTMB/) under
[Articles](https://sdmTMB.github.io/sdmTMB/articles/).**

## Notation conventions

This appendix uses the following notation conventions, which generally
follow the guidance in Edwards & Auger-Méthé (2019):

- Greek symbols for parameters,

- the Latin/Roman alphabet for data except for cases where another
  symbol is used by convention,

- bold symbols for vectors or matrices (e.g., $\mathbf{ω}$ is a vector
  and $\omega_{\mathbf{s}}$ is the value of $\mathbf{ω}$ at point in
  space $\mathbf{s}$), with Latin/Roman symbols in non-italics (e.g.,
  $\mathbf{X}$, $\mathbf{Q}$),

- $\phi$ for all distribution dispersion parameters for consistency with
  the code,

- ${\mathbb{E}}\lbrack y\rbrack$ to define the expected value (mean) of
  variable $y$,

- ${Var}\lbrack y\rbrack$ to define the expected variance of the
  variable $y$,

- a $^{*}$ superscript represents interpolated or projected values as
  opposed to values at knot locations (e.g., $\mathbf{ω}$
  vs. ${\mathbf{ω}}^{*}$), and

- where possible, notation has been chosen to match VAST (Thorson 2019)
  to maintain consistency (e.g., $\mathbf{ω}$ for spatial fields and
  ${\mathbf{ϵ}}_{t}$ for spatiotemporal fields).

Tables of indices and symbols are summarized in [notation reference
tables](#notation-reference-tables) at the end of this vignette.

## sdmTMB model structure

The complete sdmTMB model can be written as

$$\begin{aligned}
\eta_{\mathbf{s},t} & {= \mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}} + O_{\mathbf{s},t} + \mathbf{Z}_{g}\mathbf{b}_{g} + \mathbf{X}_{\mathbf{s},t}^{tvc}{\mathbf{γ}}_{t} + \mathbf{X}_{\mathbf{s},t}^{svc}{\mathbf{ζ}}_{\mathbf{s}} + \omega_{\mathbf{s}} + \epsilon_{\mathbf{s},t},} \\
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \eta_{\mathbf{s},t} \right),} \\
{{\mathbb{E}}\left\lbrack y_{\mathbf{s},t} \right\rbrack} & {= \mu_{\mathbf{s},t},}
\end{aligned}$$

where

- $y_{\mathbf{s},t}$ represents the response data at point $\mathbf{s}$
  and time $t$;
- $\eta$ represents the linear predictor (i.e., in link space)
- $\mu$ represents the conditional mean on the response/data scale
  (after applying the inverse link);
- $f$ represents a link function (e.g., log or logit) and $f^{- 1}$
  represents its inverse, and $\eta$ represents the linear predictor
  before applying $f^{- 1}$;
- $\mathbf{X}^{main}$, $\mathbf{X}^{tvc}$, and $\mathbf{X}^{svc}$
  represent design matrices (the superscript identifiers ‘main’ = main
  effects, ‘tvc’ = time varying coefficients, and ‘svc’ = spatially
  varying coefficients);
- $\mathbf{β}$ represents a vector of fixed-effect coefficients;
- $O_{\mathbf{s},t}$ represents an offset: a covariate (usually log
  transformed) with a coefficient fixed at one that enters on the linear
  predictor scale (before the inverse link);
- $\mathbf{Z}_{g}$ represents the random-effect design row(s) for group
  $g$ (a subset of the matrix $\mathbf{Z}$) corresponding to
  $\mathbf{b}_{g}$; these random effects enter on the linear predictor
  scale; the scalar special case is a random intercept
  $\alpha_{g} \sim N\left( 0,\sigma_{\alpha}^{2} \right)$, and more
  generally random slopes and intercepts are vector-valued with full
  covariance as described below;
- ${\mathbf{γ}}_{t}$ represents a vector of time-varying coefficients,
  where each coefficient $\gamma_{p,t}$ can follow a random walk or
  AR(1) process;
- ${\mathbf{ζ}}_{\mathbf{s}}$ represents a vector of spatially varying
  coefficients (random fields), where each coefficient
  $\zeta_{l,\mathbf{s}} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\zeta,l} \right)$;
- $\omega_{\mathbf{s}}$ represents a spatial component (a Gaussian
  Markov random field, GMRF),
  $\omega_{\mathbf{s}} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\omega} \right)$
  with Matérn covariance built from the mesh; and
- $\epsilon_{\mathbf{s},t}$ represents a spatiotemporal component (a
  GMRF),
  $\epsilon_{\mathbf{s},t} \sim {MVN}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right)$
  with Matérn covariance over space and temporal evolution as specified.

A single sdmTMB model will rarely, if ever, contain all of the above
components. Next, we will split the model to describe the various parts
in more detail using ‘$\ldots$’ to represent the other optional
components.

### Main effects

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}} + \ldots \right)}
\end{aligned}$$

Within
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md),
$\mathbf{X}_{\mathbf{s},t}^{main}{\mathbf{β}}$ is defined by the
`formula` argument and represents the main-effect model matrix and a
corresponding vector of coefficients. This main effect formula can
contain optional penalized smoothers or non-linear functions as defined
below.

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

While mgcv can fit unpenalized (e.g., B-splines) or penalized splines
(P-splines), sdmTMB only implements penalized splines. The penalized
splines are constructed in sdmTMB using the function
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html),
which transforms splines into random effects (and associated design
matrices) that are estimable in a mixed-effects modelling framework.
This is the same approach as is implemented in the R packages gamm4
(Wood & Scheipl 2020) and brms (Bürkner 2017).

#### Linear break-point threshold models

The linear break-point or “hockey stick” model can be used to describe
threshold or asymptotic responses. This function consists of two pieces,
so that for $x < b_{1}$, $s(x) = x \cdot b_{0}$, and for $x > b_{1}$,
$s(x) = b_{1} \cdot b_{0}$. In both cases, $b_{0}$ represents the slope
of the function up to a threshold, and the product $b_{1} \cdot b_{0}$
represents the value at the asymptote. No constraints are placed on
parameters $b_{0}$ or $b_{1}$.

These models can be fit by including `+ breakpt(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
break-point covariate.

#### Logistic threshold models

Models with logistic threshold relationships between a predictor and the
response can be fit with the form

$$s(x) = \tau + \psi\ \left\lbrack 1 + e^{- \ln{(19)} \cdot {(x - s50)}/{(s95 - s50)}} \right\rbrack^{- 1},$$

where $s$ represents the logistic function, $\psi$ is a scaling
parameter (controlling the height of the y-axis for the response;
unconstrained), $\tau$ is an intercept, $s50$ is a parameter controlling
the point at which the function reaches 50% of the maximum ($\psi$), and
$s95$ is a parameter controlling the point at which the function reaches
95% of the maximum. The parameter $s50$ is unconstrained but $s95$ is
constrained to be larger than $s50$.

These models can be fit by including `+ logistic(x)` in the model
formula, where `x` is a covariate. The formula can contain a single
logistic covariate.

### Spatial random fields

Spatial random fields, $\omega_{\mathbf{s}}$, are included if
`spatial = 'on'` (or `TRUE`) and omitted if `spatial = 'off'` (or
`FALSE`).

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \omega_{\mathbf{s}} + \ldots \right),} \\
{\mathbf{ω}} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\omega} \right),} \\
 & 
\end{aligned}$$ The marginal standard deviation of $\mathbf{ω}$ is
indicated by `Spatial SD` in the printed model output or as `sigma_O` in
the output of `sdmTMB::tidy(fit, "ran_pars")`. The ‘O’ is for ‘omega’
($\omega$).

Internally, the random fields follow a Gaussian Markov random field
(GMRF)

$${\mathbf{ω}} \sim {MVNormal}\left( \mathbf{0},\sigma_{\omega}^{2}\mathbf{Q}_{\omega}^{- 1} \right),$$
where $\mathbf{Q}_{\omega}$ is a sparse precision matrix and
$\sigma_{\omega}^{2}$ is the marginal variance.

### Spatiotemporal random fields

Spatiotemporal random fields are included by default if there are
multiple time elements (`time` argument is not `NULL`) and can be set to
IID (independent and identically distributed, `'iid'`; default), AR(1)
(`'ar1'`), random walk (`'rw'`), or off (`'off'`) via the
`spatiotemporal` argument. These text values are case insensitive.

Spatiotemporal random fields are represented by ${\mathbf{ϵ}}_{t}$
within sdmTMB. This has been chosen to match the representation in VAST
(Thorson 2019). The marginal standard deviation of ${\mathbf{ϵ}}_{t}$ is
indicated by `Spatiotemporal SD` in the printed model output or as
`sigma_E` in the output of `sdmTMB::tidy(fit, "ran_pars")`. The ‘E’ is
for ‘epsilon’ ($\epsilon$).

#### IID spatiotemporal random fields

IID spatiotemporal random fields (`spatiotemporal = 'iid'`) can be
represented as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \epsilon_{\mathbf{s},t} + \ldots \right),} \\
{\mathbf{ϵ}}_{t} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right).}
\end{aligned}$$

where $\epsilon_{\mathbf{s},t}$ represent random field deviations at
point $\mathbf{s}$ and time $t$. The random fields are assumed
independent across time steps.

Similarly to the spatial random fields, these spatiotemporal random
fields (including all versions described below) are parameterized
internally with a sparse precision matrix ($\mathbf{Q}_{\epsilon}$)

$${\mathbf{ϵ}}_{t} \sim {MVNormal}\left( \mathbf{0},\sigma_{\epsilon}^{2}\mathbf{Q}_{\epsilon}^{- 1} \right).$$

#### AR(1) spatiotemporal random fields

First-order auto regressive, AR(1), spatiotemporal random fields
(`spatiotemporal = 'ar1'`) add a parameter defining the correlation
between random field deviations from one time step to the next. They are
defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \delta_{\mathbf{s},t} + \ldots \right),} \\
{\mathbf{δ}}_{t = 1} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),} \\
{\mathbf{δ}}_{t > 1} & {= \rho{\mathbf{δ}}_{t - 1} + \sqrt{1 - \rho^{2}}{\mathbf{ϵ}}_{t},\ {\mathbf{ϵ}}_{t} \sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),}
\end{aligned}$$ where $\rho$ is the correlation between subsequent
spatiotemporal random fields. The
$\rho{\mathbf{δ}}_{t - 1} + \sqrt{1 - \rho^{2}}$ term scales the
spatiotemporal variance by the correlation such that it represents the
steady-state marginal variance. The correlation $\rho$ allows for
mean-reverting spatiotemporal fields, and is constrained to be
$- 1 < \rho < 1$. Internally, the parameter is estimated as `ar1_phi`,
which is unconstrained. The parameter `ar1_phi` is transformed to $\rho$
with
$\rho = 2\left( {logit}^{- 1}\left( \texttt{𝚊𝚛𝟷\_𝚙𝚑𝚒} \right) \right) - 1$,
mapping `ar1_phi` to the range $( - 1,1)$.

#### Random walk spatiotemporal random fields (RW)

Random walk spatiotemporal random fields (`spatiotemporal = 'rw'`)
represent a model where the difference in spatiotemporal deviations from
one time step to the next are IID. They are defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \delta_{\mathbf{s},t} + \ldots \right),} \\
{\mathbf{δ}}_{t = 1} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),} \\
{\mathbf{δ}}_{t > 1} & {= {\mathbf{δ}}_{t - 1} + {\mathbf{ϵ}}_{t},\ {\mathbf{ϵ}}_{t} \sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\epsilon} \right),}
\end{aligned}$$

where the distribution of the spatiotemporal field in the initial time
step is the same as for the AR(1) model, but the absence of the $\rho$
parameter allows the spatiotemporal field to be non-stationary in time.
Note that, in contrast to the AR(1) parametrization, the variance is no
longer the steady-state marginal variance.

### Time-varying regression parameters

Parameters can be modelled as time-varying according to a random walk or
first-order autoregressive, AR(1), process. The time-series model is
defined by `time_varying_type`. For all types:

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \mathbf{X}_{\mathbf{s},t}^{tvc}{\mathbf{γ}}_{t} + \ldots \right),}
\end{aligned}$$ where ${\mathbf{γ}}_{t}$ is an optional vector of
time-varying regression parameters and $\mathbf{X}_{\mathbf{s},t}^{tvc}$
is the corresponding model matrix with covariate values. This is defined
via the `time_varying` argument, assuming that the `time` argument is
also supplied a column name. `time_varying` takes a *one-sided* formula.
`~ 1` implies a time-varying intercept.

The following equations describe the time-series process for a single
time-varying coefficient $\gamma_{p,t}$, where $p$ indexes the
time-varying coefficient. When multiple time-varying coefficients are
specified, each follows its time-series process independently with its
own variance $\sigma_{\gamma,p}^{2}$ and (for AR1) correlation parameter
$\rho_{\gamma,p}$.

For `time_varying_type = 'rw'`, the first time step is estimated
independently:

$$\begin{aligned}
\gamma_{p,t = 1} & {\sim {Uniform}( - \infty,\infty){\mspace{6mu}\text{(treated as an unconstrained parameter)}},} \\
\gamma_{p,t > 1} & {\sim {Normal}\left( \gamma_{p,t - 1},\sigma_{\gamma,p}^{2} \right).}
\end{aligned}$$

In this case, the first time-step value is given an implicit uniform
prior. I.e., the same variable should not appear in the fixed effect
formula since the initial value is estimated as part of the time-varying
formula. The formula `time_varying = ~ 1` implicitly represents a
time-varying intercept (assuming the `time` argument has been supplied)
and, this case, the intercept should be omitted from the main effects
(`formula ~ + 0 + ...` or `formula ~ -1 + ...`).

For `time_varying_type = 'rw0'`, the first time step is estimated from a
mean-zero prior:

$$\begin{aligned}
\gamma_{p,t = 1} & {\sim {Normal}\left( 0,\sigma_{\gamma,p}^{2} \right),} \\
\gamma_{p,t > 1} & {\sim {Normal}\left( \gamma_{p,t - 1},\sigma_{\gamma,p}^{2} \right).}
\end{aligned}$$ In this case, the time-varying variable (including the
intercept) *should* be included in the main effects.

For `time_varying_type = 'ar1'`:

$$\begin{aligned}
\gamma_{p,t = 1} & {\sim {Normal}\left( 0,\sigma_{\gamma,p}^{2} \right),} \\
\gamma_{p,t > 1} & {\sim {Normal}\left( \rho_{\gamma,p}\gamma_{p,t - 1},\left( 1 - \rho_{\gamma,p}^{2} \right)\sigma_{\gamma,p}^{2} \right),}
\end{aligned}$$ where $\rho_{\gamma,p}$ is the correlation between
subsequent time steps. The first time step is given a mean-zero prior.

### Spatially varying coefficients (SVC)

Spatially varying coefficient models are defined as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \mathbf{X}_{\mathbf{s},t}^{svc}{\mathbf{ζ}}_{\mathbf{s}} + \ldots \right),} \\
\zeta_{l,\mathbf{s}} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{\zeta,l} \right),}
\end{aligned}$$

where $\zeta_{l,\mathbf{s}}$ represents the $l$-th spatially varying
coefficient field at location $\mathbf{s}$. When multiple spatially
varying coefficients are specified, each is estimated as an independent
spatial random field with its own marginal variance
$\sigma_{\zeta,l}^{2}$. Usually, $\mathbf{X}_{\mathbf{s},t}^{svc}$ would
represent a prediction matrix that is constant spatially for a given
time $t$ as defined by a one-sided formula supplied to
`spatial_varying`. For example `spatial_varying = ~ 0 + x`, where `0`
omits the intercept.

The random fields are parameterized internally with a sparse precision
matrix ($\mathbf{Q}_{\zeta}$)

$$\zeta_{l,\mathbf{s}} \sim {MVNormal}\left( \mathbf{0},\sigma_{\zeta,l}^{2}\mathbf{Q}_{\zeta}^{- 1} \right).$$

### Random intercepts and slopes

Multilevel (hierarchical) intercepts and slopes follow the familiar
lme4/glmmTMB formulation. A single grouping factor can carry an
intercept-only term or a vector of random effects (intercept plus
slopes) with a full covariance:

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + \mathbf{Z}_{g}\mathbf{b}_{g} + \ldots \right),} \\
\mathbf{b}_{g} & {\sim {MVNormal}\left( \mathbf{0},\mathbf{\Sigma}_{g} \right),}
\end{aligned}$$

where $\mathbf{Z}_{g}$ is the design row(s) for group $g$ (e.g.,
intercept and covariate values) and $\mathbf{b}_{g}$ is the
corresponding vector of random effects for that group. The scalar
special case $\left( 1|g \right)$ corresponds to $n = 1$ with
$\alpha_{g} \sim {Normal}\left( 0,\sigma_{\alpha}^{2} \right)$.

Use lme4-style syntax in `formula`, e.g. `(1 | g)` for random
intercepts, `(1 + x | g)` for correlated intercepts and slopes, or
`(1 | g1) + (1 + x | g2)` to give each grouping factor its own
(co)variance.

Internally, standard deviations are estimated on the log scale, and
correlations are represented by unconstrained Cholesky parameters (via
`UNSTRUCTURED_CORR`), yielding a positive-definite $\mathbf{\Sigma}_{g}$
per grouping factor.

### Offset terms

Offset terms can be included through the `offset` argument in
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). These
are included in the linear predictor as

$$\begin{aligned}
\mu_{\mathbf{s},t} & {= f^{- 1}\left( \ldots + O_{\mathbf{s},t} + \ldots \right),}
\end{aligned}$$

where $O_{\mathbf{s},t}$ is an offset term—a **log transformed**
variable without a coefficient (assuming a log link). The offset is
*not* included in the prediction. Therefore, if `offset` represents a
measure of effort, for example, the prediction is for one unit of effort
(`log(1) = 0`).

## Observation model families

Here we describe the main observation families that are available in
sdmTMB and comment on their parametrization, statistical properties,
utility, and code representation in sdmTMB. Families are grouped by
outcome type to make it easier to locate an appropriate observation
model.

### Bounded or binary outcomes (0–1, proportions)

#### Binomial

$$\operatorname{Binomial}(N,\mu)$$ where $N$ is the size or number of
trials, and $\mu$ is the probability of success for each trial. If
$N = 1$, the distribution becomes the Bernoulli distribution.
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaecb5a18095a320b42e2d20c4b120f5f5)
in TMB, which is numerically stable when probabilities approach 0 or 1.
Following the structure of
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), lme4, and glmmTMB, a
binomial family can be specified in one of 4 ways:

1.  the response may be a factor (and the model classifies the first
    level versus all others)
2.  the response may be binomial (0/1)
3.  the response can be a matrix of form `cbind(success, failure)`, or
4.  the response may be the observed proportions, and the `weights`
    argument is used to specify the Binomial size ($N$) parameter
    (`probability ~ ..., weights = N`).

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaee11f805f02bc1febc6d7bf0487671be).

Example: `family = binomial(link = "logit")`

#### Beta-binomial

$$\operatorname{BetaBinomial}(N,\mu,\phi)$$ where $N$ is the number of
trials, $\mu$ is the mean success probability, and $\phi$ is a precision
parameter that controls overdispersion relative to a Binomial. The
implied Beta parameters are $\alpha = \mu\phi$ and
$\beta = (1 - \mu)\phi$, and the variance is \$ \[y\] = N (1 - ), , \$
which exceeds the Binomial variance unless
$\left. \phi\rightarrow\infty \right.$. Available links are logit and
cloglog. This family is useful for overdispersed counts of
successes/failures (e.g., aggregated Bernoulli data, proportions with
extra-binomial variation).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L505-L512)
and implemented in the TMB likelihood at `src/sdmTMB.cpp` (e.g., lines
974–982).

Example: `family = betabinomial(link = "logit")`

#### Beta

$$\operatorname{Beta}\left( \mu\phi,(1 - \mu)\phi \right)$$ where $\mu$
is the mean and $\phi$ is a precision parameter. This parametrization
follows Ferrari & Cribari-Neto (2004) and the betareg R package
(Cribari-Neto & Zeileis 2010). The variance is
$\mu(1 - \mu)/(\phi + 1)$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga5324c83759d5211c7c2fbbad37fa8e59).

Example: `family = Beta(link = "logit")`

### Count data

#### Poisson

$$\operatorname{Poisson}(\mu)$$ where $\mu$ represents the mean and
${Var}\lbrack y\rbrack = \mu$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa1ed15503e1441a381102a8c4c9baaf1).

Example: `family = poisson(link = "log")`

#### Negative Binomial 2 (NB2)

$$\operatorname{NB2}(\mu,\phi)$$

where $\mu$ is the mean and $\phi$ is the dispersion parameter. The
variance scales quadratically with the mean
${Var}\lbrack y\rbrack = \mu + \mu^{2}/\phi$(Hilbe 2011). The NB2
parametrization is more commonly seen in ecology than the NB1.
Internally, the distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga76266c19046e04b651fce93aa0810351).

Example: `family = nbinom2(link = "log")`

#### Negative Binomial 1 (NB1)

$$\operatorname{NB1}(\mu,\phi)$$

where $\mu$ is the mean and $\phi$ is the dispersion parameter. The
variance scales linearly with the mean
${Var}\lbrack y\rbrack = \mu + \mu/\phi$(Hilbe 2011). Internally, the
distribution is parameterized as the [robust
version](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gaa23e3ede4669d941b0b54314ed42a75c)
in TMB.

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/sdmTMB.cpp#L577-L582)
based on NB2 and borrowed from glmmTMB.

Example: `family = nbinom1(link = "log")`

#### Negative binomial 2 mixture

This is a 2 component mixture that extends the NB2 distribution,
following mixture-distribution approaches for shoaling/aggregation in
survey data (Thorson *et al.* 2011).

$$(1 - p) \cdot \operatorname{NB2}\left( \mu_{1},\phi \right) + p \cdot \operatorname{NB2}\left( \mu_{2},\phi \right)$$

where $\mu_{1}$ is the mean of the first (smaller component) of the
distribution, $\mu_{2}$ is the mean of the larger component, and $p$
controls the contribution of each component to the mixture.

Example: `family = nbinom2_mix(link = "log")`

### Positive continuous outcomes

#### Gamma

$$\operatorname{Gamma}\left( \phi,\frac{\mu}{\phi} \right)$$ where
$\phi$ represents the Gamma shape and $\mu/\phi$ represents the scale.
The mean is $\mu$ and variance is $\mu^{2}/\phi$.

Code defined [within
TMB](https://kaskr.github.io/adcomp/group__R__style__distribution.html#gab0e2205710a698ad6a0ed39e0652c9a3).

Example: `family = Gamma(link = "log")`

#### Lognormal

sdmTMB uses the “bias-corrected” lognormal distribution where $\phi$
represents the standard deviation in log-space:

$$\operatorname{Lognormal}\left( \log\mu - \frac{\phi^{2}}{2},\phi^{2} \right).$$
Because of the bias correction, ${\mathbb{E}}\lbrack y\rbrack = \mu$ and
${Var}\left\lbrack \log y \right\rbrack = \phi^{2}$.

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L47-L54)
based on the TMB [`dnorm()`](https://rdrr.io/r/stats/Normal.html) normal
density.

Example: `family = lognormal(link = "log")`

#### Generalized gamma

$$\operatorname{GenGamma}\left( \mu,\phi,Q^{gg} \right)$$

sdmTMB implements the Prentice (1974) parameterization introduced for
spatiotemporal models and index standardization by Dunic *et al.*
(2025), with parameters mean $\mu$, scale $\phi$ (standard deviation on
the log scale), and a shape parameter $Q^{gg}$ (reported as
`Generalized gamma Q`).

Here, $\mu$ is the mean on the data scale, $\phi$ is the log-scale
standard deviation, and $Q^{gg}$ controls the shape
($\left. Q^{gg}\rightarrow 0 \right.$ yields the lognormal;
$Q^{gg} = \phi$ yields the gamma). See Dunic *et al.* (2025) for the
full PDF in the Prentice formulation as implemented. Links available:
identity, log, inverse. This flexibility is useful for right-skewed
positive responses with tails heavier or lighter than gamma/lognormal,
and is often paired in a hurdle/mixture as
[`delta_gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
for zero-inflated biomass or catch data (see Dunic *et al.* (2025)).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L7-L43).

Example: `family = gengamma(link = "log")`; delta/hurdle:
`family = delta_gengamma(link1 = "logit", link2 = "log")`.

#### Tweedie

$$\operatorname{Tweedie}(\mu,p,\phi),\ 1 < p < 2$$

where $\mu$ is the mean, $p$ is the power parameter constrained between
1 and 2, and $\phi$ is the dispersion parameter. The Tweedie
distribution can be helpful for modelling data that are positive and
continuous but also contain zeros.

Internally, $p$ is transformed from
${logit}^{- 1}\left( \texttt{𝚝𝚑𝚎𝚝𝚊𝚏} \right) + 1$ to constrain it
between 1 and 2 and is estimated as an unconstrained variable.

The [source
code](https://kaskr.github.io/adcomp/tweedie_8cpp_source.html) is
implemented as in the [cplm](https://CRAN.R-project.org/package=cplm)
package (Zhang 2013) and is based on Dunn & Smyth (2005). The TMB
version is defined
[here](https://kaskr.github.io/adcomp/group__R__style__distribution.html#ga262f3c2d1cf36f322a62d902a608aae0).

Example: `family = tweedie(link = "log")`

#### Gamma mixture

This is a 2 component mixture that extends the Gamma distribution,
motivated by mixture-distribution treatments of aggregation in survey
data (Thorson *et al.* 2011),

$$(1 - p) \cdot \operatorname{Gamma}\left( \phi,\frac{\mu_{1}}{\phi} \right) + p \cdot \operatorname{Gamma}\left( \phi,\frac{\mu_{2}}{\phi} \right),$$
where $\phi$ represents the Gamma shape, $\mu_{1}/\phi$ represents the
scale for the first (smaller component) of the distribution,
$\mu_{2}/\phi$ represents the scale for the second (larger component) of
the distribution, and $p$ controls the contribution of each component to
the mixture (also interpreted as the probability of larger events).

The mean is $(1 - p) \cdot \mu_{1} + p \cdot \mu_{2}$. The variance
follows the usual mixture formula:
$(1 - p)\left( \mu_{1}^{2}/\phi + \mu_{1}^{2} \right) + p\left( \mu_{2}^{2}/\phi + \mu_{2}^{2} \right) - \left\lbrack (1 - p)\mu_{1} + p\mu_{2} \right\rbrack^{2}$.

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
again in the spirit of mixture approaches for aggregating/ shoaling data
(Thorson *et al.* 2011),

$$(1 - p) \cdot \operatorname{Lognormal}\left( \log\mu_{1} - \frac{\phi^{2}}{2},\phi^{2} \right) + p \cdot \operatorname{Lognormal}\left( \log\mu_{2} - \frac{\phi^{2}}{2},\phi^{2} \right).$$

Because of the bias correction,
${\mathbb{E}}\lbrack y\rbrack = (1 - p) \cdot \mu_{1} + p \cdot \mu_{2}$.
The log-scale variance of the mixture is not simply $\phi^{2}$; it can
be obtained with the standard mixture-variance formula using component
log-means $\log\mu_{i} - \phi^{2}/2$ and log-variance $\phi^{2}$.

As with the Gamma mixture, $p$ controls the contribution of each
component to the mixture (also interpreted as the probability of larger
events).

Example: `family = lognormal_mix(link = "log")`. See also
`family = delta_lognormal_mix()` for an extension incorporating this
distribution with delta models. Like with the gamma mixture, fixed
probabilities of extreme events ($p$ in notation above) can be passed
in, e.g.

``` r
sdmTMB(...,
  family = delta_lognormal_mix(p_extreme = 0.01)
)
```

### Continuous real-valued outcomes

#### Gaussian

$$\operatorname{Normal}\left( \mu,\phi^{2} \right)$$ where $\mu$ is the
mean and $\phi$ is the standard deviation. The variance is $\phi^{2}$.

Example: `family = Gaussian(link = "identity")`

Code defined [within
TMB](https://kaskr.github.io/adcomp/dnorm_8hpp.html).

#### Student-t

$$\operatorname{Student-t}(\mu,\phi,\nu)$$

where $\nu$, the degrees of freedom (`df`), is a user-supplied fixed
parameter. Lower values of $\nu$ result in heavier tails compared to the
Gaussian distribution. Above approximately `df = 20`, the distribution
becomes very similar to the Gaussian. The Student-t distribution with a
low degrees of freedom (e.g., $\nu \leq 7$) can be helpful for modelling
data that would otherwise be suitable for Gaussian but needs an approach
that is robust to outliers (e.g., Anderson *et al.* 2017).

Code defined [within
sdmTMB](https://github.com/sdmTMB/sdmTMB/blob/18a39eabc111e2179fa589f942c8820d87ad10df/src/utils.h#L37-L45)
based on the [`dt()`](https://rdrr.io/r/stats/TDist.html) distribution
in TMB.

Example: `family = student(link = "log", df = 7)`

### Zero-inflated and hurdle structures

#### Delta models

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
Default links associated with each family can be inspected with
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

Or for we may want separate main-effects formulas. For example:

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

## Gaussian random fields

### Matérn parameterization

The Matérn defines the covariance $\Phi\left( s_{j},s_{k} \right)$
between spatial locations $s_{j}$ and $s_{k}$ as

$$\Phi\left( s_{j},s_{k} \right) = \tau^{2}\left\lbrack 2^{\nu - 1}\Gamma(\nu) \right\rbrack^{- 1}\left( \kappa d_{jk} \right)^{\nu}K_{\nu}\left( \kappa d_{jk} \right),$$

where $\tau^{2}$ is a scale/precision factor (higher $\tau$ means lower
marginal variance), $\nu$ controls the smoothness, $\Gamma$ represents
the Gamma function, $d_{jk}$ represents the distance between locations
$s_{j}$ and $s_{k}$, $K_{\nu}$ represents the modified Bessel function
of the second kind, and $\kappa$ represents the decorrelation rate. The
parameter $\nu$ is set to 1 because the SPDE approach that gives us
sparse matrices (and therefore major computational speedups) requires
$\alpha = \nu + d/2$ to be an integer when working in 2D (so $\nu = 1$)
(Lindgren *et al.* 2011). Internally, the parameters $\kappa$ and $\tau$
are converted to range and marginal standard deviation $\sigma$ as
$\text{range} = \sqrt{8}/\kappa$ and
$\sigma = 1/\sqrt{4\pi\exp\left( 2\log(\tau) + 2\log(\kappa) \right)}$,
so increasing $\kappa$ or $\tau$ decreases $\sigma$.

In the case of a spatiotemporal model with both spatial and
spatiotemporal fields, if `share_range = TRUE` in
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) (the
default), then a single $\kappa$ and range are estimated with separate
$\sigma_{\omega}$ and $\sigma_{\epsilon}$. This often makes sense since
data are often only weakly informative about $\kappa$. If
`share_range = FALSE`, then separate $\kappa_{\omega}$ and
$\kappa_{\epsilon}$ are estimated. The spatially varying coefficient
field always shares $\kappa$ with the spatial random field.

### Projection $\mathbf{A}$ matrix

The values of the spatial variables at the knots are multiplied by a
projection matrix $\mathbf{A}$ that bilinearly interpolates from the
knot locations to the values at the locations of the observed or
predicted data (Lindgren & Rue 2015)

$${\mathbf{ω}}^{*} = \mathbf{A}{\mathbf{ω}},$$ where ${\mathbf{ω}}^{*}$
represents the values of the spatial random fields at the observed
locations or predicted data locations. The matrix $\mathbf{A}$ has a row
for each data point or prediction point and a column for each knot.
Three non-zero elements on each row define the weight of the
neighbouring 3 knot locations for location $\mathbf{s}$. The same
bilinear interpolation happens for any spatiotemporal random fields

$${\mathbf{ϵ}}_{t}^{*} = \mathbf{A}{\mathbf{ϵ}}_{t}.$$

### Anisotropy

TMB allows for anisotropy, where spatial covariance may be asymmetric
with respect to latitude and longitude ([full
details](https://kaskr.github.io/adcomp/namespaceR__inla.html)).
Anisotropy can be turned on or off with the logical `anisotropy`
argument to
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). There
are a number of ways to implement anisotropic covariance (Fuglstad *et
al.* 2015), and we adopt a 2-parameter rotation matrix $\mathbf{H}$. The
elements of $\mathbf{H}$ are defined by the parameter vector
$\mathbf{x}$ so that $H_{1,1} = x_{1}$, $H_{1,2} = H_{2,1} = x_{2}$ and
$H_{2,2} = \left( 1 + x_{2}^{2} \right)/x_{1}$.

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
constructs a mesh with a call to `sdmTMBextra::make_mesh()`, and (2)
passes the mesh to
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

[This website](https://haakonbakkagit.github.io/btopic128.html) by
Francesco Serafini and Haakon Bakka provides an illustration with INLA.
The implementation within TMB was borrowed from code written by Olav
Nikolai Breivik and Hans Skaug at the [TMB Case
Studies](https://github.com/skaug/tmb-case-studies) Github site.

## Optimization

### Optimization details

The sdmTMB model is fit by maximum marginal likelihood. Internally, a
TMB (Kristensen *et al.* 2016) model template calculates the marginal
log likelihood and its gradient, and the negative log likelihood is
minimized via the non-linear optimization routine
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) in R (Gay 1990;
R Core Team 2021). Random effects are represented by the values that
maximize the log likelihood conditional on the fixed effects (their
conditional modes), and the Laplace approximation integrates over them
when evaluating the likelihood (Kristensen *et al.* 2016).

Like AD Model Builder (Fournier *et al.* 2012), TMB allows for
parameters to be fit in phases and we include the `multiphase` argument
in
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
to allow this. For high-dimensional models (many fixed and random
effects), phased estimation may be faster, but this is not always the
case because it requires an extra call to
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html). In
sdmTMB, phased estimation proceeds by first estimating all fixed-effect
parameters contributing to the likelihood (holding random effects
constant at initial values). In the second phase, the random-effect
parameters (and their variances) are also estimated. Fixed-effect
parameters are also estimated in the second phase and are initialized at
their estimates from the first phase.

In some cases, a single call to
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) may not be
result in convergence (e.g., the maximum gradient of the marginal
likelihood with respect to fixed-effect parameters is not small enough
yet), and the algorithm may need to be run multiple times. In the
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
function, we include an argument `nlminb_loops` that will restart the
optimization at the previous best values. The number of `nlminb_loops`
should generally be small (e.g., 2 or 3), and defaults to 1. For some
sdmTMB models, the Hessian may also be unstable and need to be
re-evaluated. We do this optionally with the
[`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) routine after
the call to [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).
The [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) function
implements a Newton optimization routine to find the Hessian, and we
include the argument `newton_loops` in
[`sdmTMB::sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
to allow for multiple function evaluations (each starting at the
previous best value). By default, this is set to 1 evaluation. The
updated parameters are only accepted if they result in a lower negative
marginal log likelihood. If a model is already fit, the function
[`sdmTMB::run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.md)
can run additional optimization loops with either routine to further
reduce the maximum gradient.

### Assessing convergence

Much of the guidance around diagnostics and glmmTMB also applies to
sdmTMB, e.g. [the glmmTMB vignette on
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
the marginal likelihood with respect to fixed-effect parameters can be
checked by inspecting (`fit$gradients`). Guidance varies, but gradients
should be very close to zero (e.g., on the order of $10^{- 3}$ or
smaller after reasonable parameter scaling) before assuming the fitting
routine is consistent with convergence. If maximum gradients are already
relatively small, they can sometimes be reduced further with additional
optimization calls beginning at the previous best parameter vector as
described above with
[`sdmTMB::run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.md).

## Notation reference tables

Tables of all major indices (Table 1) and symbols (Table 2) are provided
here for quick reference.

| Symbol       | Description                                      |
|:-------------|:-------------------------------------------------|
| $\mathbf{s}$ | Index for space; a vector of x and y coordinates |
| $t$          | Index for time                                   |
| $g$          | Group                                            |
| $p$          | Index for time-varying coefficient               |
| $l$          | Index for spatially varying coefficient (SVC)    |

Table 1: Subscript notation

| Symbol                        | Code                 | Description                                                                                                 |
|:------------------------------|:---------------------|:------------------------------------------------------------------------------------------------------------|
| $y$                           | `y_i`                | Observed response data                                                                                      |
| $\mu$                         | `mu_i`               | Mean                                                                                                        |
| $\eta$                        | `eta_i`              | Linear predictor before applying the inverse link ($f^{- 1}$)                                               |
| $\phi$                        | `phi`                | A dispersion parameter for a distribution                                                                   |
| $f$                           | `fit$family$link`    | Link function                                                                                               |
| $f^{- 1}$                     | `fit$family$linkinv` | Inverse link function                                                                                       |
| $\mathbf{β}$                  | `b_j`                | Parameter vector                                                                                            |
| $\mathbf{X}$                  | `X_ij`               | A predictor model matrix                                                                                    |
| $\mathbf{Z}$                  | `Zt_list`            | Random effect design matrix (list of sparse blocks; rows correspond to observations; subset by group $g$)   |
| $O_{\mathbf{s},t}$            | `offset`             | An offset variable at point $\mathbf{s}$ and time $t$                                                       |
| $\omega_{\mathbf{s}}$         | `omega_s`            | Spatial random field at point $\mathbf{s}$ (knot)                                                           |
| $\omega_{\mathbf{s}}^{*}$     | `omega_s_A`          | Spatial random field at point $\mathbf{s}$ (interpolated)                                                   |
| $\zeta_{\mathbf{s}}$          | `zeta_s`             | Spatially varying coefficient random field at point $\mathbf{s}$ (knot)                                     |
| $\zeta_{\mathbf{s}}^{*}$      | `zeta_s_A`           | Spatially varying coefficient random field at point $\mathbf{s}$ (interpolated)                             |
| $\epsilon_{\mathbf{s},t}$     | `epsilon_st`         | Spatiotemporal random field at point $\mathbf{s}$ and time $t$ (knot)                                       |
| $\epsilon_{\mathbf{s},t}^{*}$ | `epsilon_st_A`       | Spatiotemporal random field at point $\mathbf{s}$ and time $t$ (interpolated)                               |
| $\delta_{\mathbf{s},t}$       | `b_t`                | AR(1) or random walk spatiotemporal deviations (knot)                                                       |
| $\alpha_{g}$                  | `RE`                 | IID random intercept deviation for group $g$                                                                |
| $\mathbf{b}_{g}$              | `re_b_pars`          | Vector of random intercepts/slopes for group $g$                                                            |
| $\mathbf{\Sigma}_{\omega}$    | `-`                  | Spatial random field covariance matrix                                                                      |
| $\mathbf{\Sigma}_{\zeta}$     | `-`                  | Spatially varying coefficient random field covariance matrix                                                |
| $\mathbf{\Sigma}_{\epsilon}$  | `-`                  | Spatiotemporal random field covariance matrix                                                               |
| $\mathbf{Q}_{\omega}$         | `Q_s`                | Spatial random field precision matrix                                                                       |
| $\mathbf{Q}_{\zeta}$          | `Q_s`                | Spatially varying coefficient random field precision matrix                                                 |
| $\mathbf{Q}_{\epsilon}$       | `Q_st`               | Spatiotemporal random field precision matrix                                                                |
| $\mathbf{\Sigma}_{g}$         | `re_cov_pars`        | Random effect covariance per grouping factor (Cholesky SD + correlation params)                             |
| $\sigma_{\alpha}^{2}$         | `sigma_G`            | IID random intercept variance                                                                               |
| $\sigma_{\epsilon}^{2}$       | `sigma_E`            | Spatiotemporal random field marginal variance                                                               |
| $\sigma_{\omega}^{2}$         | `sigma_O`            | Spatial random field marginal variance                                                                      |
| $\sigma_{\zeta}^{2}$          | `sigma_Z`            | Spatially varying coefficient random field marginal variance                                                |
| $\kappa_{\omega}$             | `kappa(0)`           | Spatial decorrelation rate                                                                                  |
| $\kappa_{\epsilon}$           | `kappa(1)`           | Spatiotemporal decorrelation rate                                                                           |
| $\rho$                        | `ar1_rho`            | Correlation between random fields in subsequent time steps                                                  |
| $\rho_{\gamma}$               | `rho_time`           | Correlation between time-varying coefficients in subsequent time steps                                      |
| $Q^{gg}$                      | `gengamma_Q`         | Generalized gamma shape parameter ($\left. Q\rightarrow 0 \right.$ gives lognormal; $Q = \phi$ gives gamma) |
| $\mathbf{A}$                  | `A`                  | Sparse projection matrix to interpolate between knot and data locations                                     |
| $\mathbf{H}$                  | `H`                  | 2-parameter rotation matrix used to define anisotropy                                                       |

Table 2: Symbol notation, code representation (in model output or in
model template code), and descriptions.

## References

Anderson, S.C., Branch, T.A., Cooper, A.B. & Dulvy, N.K. (2017).
Black-swan events in animal populations. *Proceedings of the National
Academy of Sciences*, **114**, 3252–3257.

Anderson, S.C., Ward, E.J., English, P.A., Barnett, L.A.K. & Thorson,
J.T. (2024). [sdmTMB: An R package for fast, flexible, and user-friendly
generalized linear mixed effects models with spatial and spatiotemporal
random fields](https://doi.org/10.1101/2022.03.24.485545). *bioRxiv*,
**2022.03.24.485545**.

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

Lindgren, F. & Rue, H. (2015). [Bayesian spatial modelling with
R-INLA](https://doi.org/10.18637/jss.v063.i19). *Journal of Statistical
Software*, **63**, 1–25.

Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between
Gaussian fields and Gaussian Markov random fields: The stochastic
partial differential equation approach. *J. R. Stat. Soc. Ser. B Stat.
Methodol.*, **73**, 423–498.

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

Thorson, J.T., Stewart, I.J. & Punt, A.E. (2011). [Accounting for fish
shoals in single- and multi-species survey data using mixture
distribution models](https://doi.org/10.1139/f2011-086). *Canadian
Journal of Fisheries and Aquatic Sciences*, **68**, 1681–1693.

Wood, S.N. (2017). *Generalized additive models: An introduction with
R*, 2nd edn. Chapman and Hall/CRC.

Wood, S. & Scheipl, F. (2020). *gamm4: Generalized Additive Mixed Models
using ’mgcv’ and ’lme4’*. Retrieved from
<https://CRAN.R-project.org/package=gamm4>

Zhang, Y. (2013). [Likelihood-based and Bayesian methods for Tweedie
compound Poisson linear mixed
models](https://doi.org/10.1007/s11222-012-9343-7). *Statistics and
Computing*, **23**, 743–757.
