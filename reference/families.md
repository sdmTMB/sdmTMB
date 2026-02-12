# Additional families

Additional families compatible with
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

## Usage

``` r
Beta(link = "logit")

lognormal(link = "log")

gengamma(link = "log")

gamma_mix(link = "log", p_extreme = NULL)

lognormal_mix(link = "log", p_extreme = NULL)

nbinom2_mix(link = "log", p_extreme = NULL)

nbinom2(link = "log")

nbinom1(link = "log")

truncated_nbinom2(link = "log")

truncated_nbinom1(link = "log")

student(link = "identity", df = NULL)

tweedie(link = "log")

censored_poisson(link = "log")

delta_gamma(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_gamma_mix(link1 = "logit", link2 = "log", p_extreme = NULL)

delta_gengamma(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_lognormal(link1, link2 = "log", type = c("standard", "poisson-link"))

delta_lognormal_mix(
  link1,
  link2 = "log",
  type = c("standard", "poisson-link"),
  p_extreme = NULL
)

delta_truncated_nbinom2(link1 = "logit", link2 = "log")

delta_truncated_nbinom1(link1 = "logit", link2 = "log")

delta_poisson_link_gamma(link1 = "log", link2 = "log")

delta_poisson_link_lognormal(link1 = "log", link2 = "log")

betabinomial(link = "logit")

delta_beta(link1 = "logit", link2 = "logit")
```

## Arguments

- link:

  Link.

- p_extreme:

  Optional fixed probability for the extreme component. If NULL
  (default), this is estimated. If specified, must be a proportion
  between 0 and 1.

- df:

  Student-t degrees of freedom parameter. Can be `NULL` to estimate
  (default) or a numeric value \> 1 to fix at a specific value.

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

The default `link1` for delta models of `type = "standard"` is
`"logit"`. The default `link1` for delta models of
`type = "poisson-link"` is `"log"`.

`delta_poisson_link_gamma()` and `delta_poisson_link_lognormal()` have
been deprecated in favour of `delta_gamma(type = "poisson-link")` and
`delta_lognormal(type = "poisson-link")`.

The `gengamma()` family was implemented by J.T. Thorson and uses the
Prentice (1974) parameterization such that the lognormal occurs as the
internal parameter `gengamma_Q` (reported in
[`print()`](https://rdrr.io/r/base/print.html) or
[`summary()`](https://rdrr.io/r/base/summary.html) as "Generalized gamma
Q") approaches 0. If Q matches `phi` the distribution should be the
gamma.

The families ending in `_mix()` are 2-component mixtures where each
distribution has its own mean but a shared scale parameter. (Thorson et
al. 2011). See the model-description vignette for details. The parameter
`p_extreme = plogis(logit_p_extreme)` is the probability of the extreme
(larger) mean and `exp(log_ratio_mix) + 1` is the ratio of the larger
extreme mean to the "regular" mean. You can see these parameters in
`model$sd_report`. The parameter `p_extreme` can be fixed a priori and
passed in as a proportion for these families.

The `nbinom2` negative binomial parameterization is the NB2 where the
variance grows quadratically with the mean (Hilbe 2011).

The `nbinom1` negative binomial parameterization lets the variance grow
linearly with the mean (Hilbe 2011).

For `student()`, the degrees of freedom parameter is estimated by
default (`df = NULL`). You can fix it at a specific value by providing a
number \> 1 (e.g., `df = 3`).

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

*Families ending in `_mix()`*:

Thorson, J.T., Stewart, I.J., and Punt, A.E. 2011. Accounting for fish
shoals in single- and multi-species survey data using mixture
distribution models. Can. J. Fish. Aquat. Sci. 68(9): 1681–1693.
[doi:10.1139/f2011-086](https://doi.org/10.1139/f2011-086) .

*Negative binomial families*:

Hilbe, J. M. 2011. Negative binomial regression. Cambridge University
Press.

*Poisson-link delta families*:

Thorson, J.T. 2018. Three problems with the conventional delta-model for
biomass sampling data, and a computationally efficient alternative.
Canadian Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
[doi:10.1139/cjfas-2017-0266](https://doi.org/10.1139/cjfas-2017-0266)

## Examples

``` r
Beta(link = "logit")
#> 
#> Family: Beta 
#> Link function: logit 
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
tweedie(link = "log")
#> 
#> Family: tweedie 
#> Link function: log 
#> 
censored_poisson(link = "log")
#> 
#> Family: censored_poisson 
#> Link function: log 
#> 
delta_gamma()
#> 
#> Family: binomial Gamma 
#> Link function: logit log 
#> 
delta_gamma_mix()
#> 
#> Family: binomial gamma_mix 
#> Link function: logit log 
#> 
delta_gengamma()
#> 
#> Family: binomial gengamma 
#> Link function: logit log 
#> 
delta_lognormal()
#> 
#> Family: binomial lognormal 
#> Link function: logit log 
#> 
delta_lognormal_mix()
#> 
#> Family: binomial lognormal_mix 
#> Link function: logit log 
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
betabinomial(link = "logit")
#> 
#> Family: betabinomial 
#> Link function: logit 
#> 
delta_beta()
#> 
#> Family: binomial Beta 
#> Link function: logit logit 
#> 
```
