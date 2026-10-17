# Additional families

Additional families compatible with
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md). In
addition to [`gaussian()`](https://rdrr.io/r/stats/family.html),
[`Gamma()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html), and
[`poisson()`](https://rdrr.io/r/stats/family.html), sdmTMB provides:

- Continuous: `student()`, `lognormal()`, `gengamma()`

- Proportions: `Beta()`, `ordbeta()`

- Non-negative with exact zeros: `tweedie()`

- Counts: `nbinom2()`, `nbinom1()`, `truncated_nbinom2()`,
  `truncated_nbinom1()`, `betabinomial()`

- Censored counts: `censored_poisson()`, `censored_nbinom2()`,
  `censored_nbinom1()`, `censored_binomial()`, `censored_betabinomial()`

- Two-component mixtures: `gamma_mix()`, `lognormal_mix()`,
  `nbinom2_mix()`

- Delta/hurdle: `delta_gamma()`, `delta_lognormal()`,
  `delta_gengamma()`, `delta_beta()`, `delta_truncated_nbinom2()`,
  `delta_truncated_nbinom1()`, `delta_gamma_mix()`,
  `delta_lognormal_mix()`

See 'Binomial families', 'Censored families', and 'Delta/hurdle models'
in the Details of
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

## Usage

``` r
student(link = "identity", df = NULL)

lognormal(link = "log")

gengamma(link = "log")

Beta(link = "logit")

ordbeta(link = "logit")

tweedie(link = "log")

nbinom2(link = "log")

nbinom1(link = "log")

truncated_nbinom2(link = "log")

truncated_nbinom1(link = "log")

betabinomial(link = "logit")

censored_poisson(link = "log")

censored_nbinom2(link = "log")

censored_nbinom1(link = "log")

censored_binomial(link = "logit")

censored_betabinomial(link = "logit")

gamma_mix(link = "log", p_extreme = NULL)

lognormal_mix(link = "log", p_extreme = NULL)

nbinom2_mix(link = "log", p_extreme = NULL)

delta_gamma(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_lognormal(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_gengamma(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_beta(link1 = "logit", link2 = "logit")

delta_truncated_nbinom2(link1 = "logit", link2 = "log")

delta_truncated_nbinom1(link1 = "logit", link2 = "log")

delta_gamma_mix(link1 = "logit", link2 = "log", p_extreme = NULL)

delta_lognormal_mix(
  link1,
  link2 = "log",
  type = c("standard", "poisson-link"),
  p_extreme = NULL
)
```

## Arguments

- link:

  Link.

- df:

  Student-t degrees of freedom parameter. Can be `NULL` to estimate
  (default) or a numeric value \> 1 to fix at a specific value.

- p_extreme:

  Optional fixed probability for the extreme component. If NULL
  (default), this is estimated. If specified, must be a proportion
  between 0 and 1.

- link1:

  Link for first part of delta/hurdle model. Defaults to `"logit"` for
  `type = "standard"` and `"log"` for `type = "poisson-link"`.

- link2:

  Link for second part of delta/hurdle model.

- type:

  Delta/hurdle family type. `"standard"` for a classic hurdle model.
  `"poisson-link"` for a Poisson-link delta model (Thorson 2018).

## Value

A list with elements common to standard R family objects including
`family`, `link`, `linkfun`, and `linkinv`. Delta/hurdle model families
also have elements `delta` (logical) and `type` (standard vs.
Poisson-link).

## Details

For `student()`, the degrees of freedom parameter is estimated by
default (`df = NULL`). You can fix it at a specific value by providing a
number \> 1 (e.g., `df = 3`).

The `gengamma()` family was implemented by J.T. Thorson and uses the
Prentice (1974) parameterization such that the lognormal occurs as the
internal parameter `gengamma_Q` (reported in
[`print()`](https://rdrr.io/r/base/print.html) or
[`summary()`](https://rdrr.io/r/base/summary.html) as "Generalized gamma
Q") approaches 0. If Q matches `phi` the distribution should be the
gamma.

The `ordbeta()` family implements the ordered-beta regression of Kubinec
(2023) for continuous responses on the closed unit interval `[0, 1]`
with point masses at the endpoints. It is a parsimonious alternative to
zero-one-inflated beta models because the three components (zeros, ones,
and continuous (0, 1) values) share a single linear predictor. Two
internal cutpoints are estimated and reported on the response scale
([`plogis()`](https://rdrr.io/r/stats/Logistic.html)).

The `nbinom2` negative binomial parameterization is the NB2 where the
variance grows quadratically with the mean (Hilbe 2011).

The `nbinom1` negative binomial parameterization lets the variance grow
linearly with the mean (Hilbe 2011).

The `censored_*()` families treat some counts as known only to lie
within a range, e.g., to account for hook competition in longline
surveys (Watson et al. 2023). Supply the bounds with
`sdmTMB(censored_upper = ...)`.

- `censored_poisson()`, `censored_nbinom1()`, and `censored_nbinom2()`
  have no upper limit on the count; use an offset of log hooks for
  effort.

- `censored_binomial()` and `censored_betabinomial()` take the number of
  trials (e.g., hooks) via `weights`, which caps the count.

All but `censored_poisson()` need the RTMB backend (the default). See
`censored_method` in
[`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
for how `censored_betabinomial()` computes its likelihood, and the [hook
competition
article](https://sdmTMB.github.io/sdmTMB/articles/hook-competition.html).

The families ending in `_mix()` are 2-component mixtures where each
distribution has its own mean but they share a scale parameter. (Thorson
et al. 2011). See the model-description vignette for details. The
parameter `p_extreme = plogis(logit_p_extreme)` is the probability of
the extreme (larger) mean and `exp(log_ratio_mix) + 1` is the ratio of
the larger extreme mean to the "regular" mean. You can see these
parameters in `model$sd_report`. The parameter `p_extreme` can be fixed
a priori and passed in as a proportion for these families.

The default first-component link (`link1`) for delta models of
`type = "standard"` is `"logit"`. For `type = "poisson-link"`, the
default `link1` is `"log"`.

## References

*Generalized gamma family*:

Prentice, R.L. 1974. A log gamma model and its maximum likelihood
estimation. Biometrika 61(3): 539–544.
[doi:10.1093/biomet/61.3.539](https://doi.org/10.1093/biomet/61.3.539)

Stacy, E.W. 1962. A Generalization of the Gamma Distribution. The Annals
of Mathematical Statistics 33(3): 1187–1192. Institute of Mathematical
Statistics.

Dunic, J.C., Conner, J., Anderson, S.C., and Thorson, J.T. 2025. The
generalized gamma is a flexible distribution that outperforms
alternatives when modelling catch rate data. ICES Journal of Marine
Science 82(4): fsaf040.
[doi:10.1093/icesjms/fsaf040](https://doi.org/10.1093/icesjms/fsaf040) .

*Ordered beta regression*:

Kubinec, R. 2023. Ordered beta regression: a parsimonious, well-fitting
model for continuous data with lower and upper bounds. Political
Analysis 31(4): 519–536.
[doi:10.1017/pan.2022.20](https://doi.org/10.1017/pan.2022.20)

*Negative binomial families*:

Hilbe, J. M. 2011. Negative binomial regression. Cambridge University
Press.

*Censored families*:

Watson, J., Edwards, A.M., and Auger-Méthé, M. 2023. A statistical
censoring approach accounts for hook competition in abundance indices
from longline surveys. Canadian Journal of Fisheries and Aquatic
Sciences. 80(3): 468–486.
[doi:10.1139/cjfas-2022-0159](https://doi.org/10.1139/cjfas-2022-0159)

*Families ending in `_mix()`*:

Thorson, J.T., Stewart, I.J., and Punt, A.E. 2011. Accounting for fish
shoals in single- and multi-species survey data using mixture
distribution models. Can. J. Fish. Aquat. Sci. 68(9): 1681–1693.
[doi:10.1139/f2011-086](https://doi.org/10.1139/f2011-086) .

*Poisson-link delta families*:

Thorson, J.T. 2018. Three problems with the conventional delta-model for
biomass sampling data, and a computationally efficient alternative.
Canadian Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
[doi:10.1139/cjfas-2017-0266](https://doi.org/10.1139/cjfas-2017-0266)

## Examples

``` r
student(link = "identity") # estimate df
#> Student-t degrees of freedom parameter will be estimated. This used to be fixed
#> at 3 by default. To fix it, supply a value to `df` (e.g., `df = 3`).
#> 
#> Family: student 
#> Link function: identity 
#> 
student(link = "identity", df = 3) # fix df at 3
#> Student-t degrees of freedom parameter fixed at 3. To estimate it, set `df =
#> NULL`.
#> 
#> Family: student 
#> Link function: identity 
#> 
lognormal(link = "log")
#> 
#> Family: lognormal 
#> Link function: log 
#> 
gengamma(link = "log")
#> 
#> Family: gengamma 
#> Link function: log 
#> 
Beta(link = "logit")
#> 
#> Family: Beta 
#> Link function: logit 
#> 
ordbeta(link = "logit")
#> 
#> Family: ordbeta 
#> Link function: logit 
#> 
tweedie(link = "log")
#> 
#> Family: tweedie 
#> Link function: log 
#> 
nbinom2(link = "log")
#> 
#> Family: nbinom2 
#> Link function: log 
#> 
nbinom1(link = "log")
#> 
#> Family: nbinom1 
#> Link function: log 
#> 
truncated_nbinom2(link = "log")
#> 
#> Family: truncated_nbinom2 
#> Link function: log 
#> 
truncated_nbinom1(link = "log")
#> 
#> Family: truncated_nbinom1 
#> Link function: log 
#> 
betabinomial(link = "logit")
#> 
#> Family: betabinomial 
#> Link function: logit 
#> 
censored_poisson(link = "log")
#> 
#> Family: censored_poisson 
#> Link function: log 
#> 
censored_nbinom2(link = "log")
#> 
#> Family: censored_nbinom2 
#> Link function: log 
#> 
censored_nbinom1(link = "log")
#> 
#> Family: censored_nbinom1 
#> Link function: log 
#> 
censored_binomial(link = "cloglog")
#> 
#> Family: censored_binomial 
#> Link function: cloglog 
#> 
censored_betabinomial(link = "cloglog")
#> 
#> Family: censored_betabinomial 
#> Link function: cloglog 
#> 
gamma_mix(link = "log")
#> 
#> Family: gamma_mix 
#> Link function: log 
#> 
lognormal_mix(link = "log")
#> 
#> Family: lognormal_mix 
#> Link function: log 
#> 
nbinom2_mix(link = "log")
#> 
#> Family: nbinom2_mix 
#> Link function: log 
#> 
delta_gamma()
#> 
#> Family: binomial Gamma 
#> Link function: logit log 
#> 
delta_lognormal()
#> 
#> Family: binomial lognormal 
#> Link function: logit log 
#> 
delta_gengamma()
#> 
#> Family: binomial gengamma 
#> Link function: logit log 
#> 
delta_beta()
#> 
#> Family: binomial Beta 
#> Link function: logit logit 
#> 
delta_truncated_nbinom2()
#> 
#> Family: binomial truncated_nbinom2 
#> Link function: logit log 
#> 
delta_truncated_nbinom1()
#> 
#> Family: binomial truncated_nbinom1 
#> Link function: logit log 
#> 
delta_gamma_mix()
#> 
#> Family: binomial gamma_mix 
#> Link function: logit log 
#> 
delta_lognormal_mix()
#> 
#> Family: binomial lognormal_mix 
#> Link function: logit log 
#> 
```
