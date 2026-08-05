# Extract a relative biomass/abundance index, center of gravity, effective area occupied, or weighted average

Extract a relative biomass/abundance index, center of gravity, effective
area occupied, or weighted average

## Usage

``` r
get_index(
  obj,
  bias_correct = TRUE,
  level = 0.95,
  area = 1,
  silent = TRUE,
  derived_link = NULL,
  ...
)

get_index_split(
  obj,
  newdata,
  bias_correct = FALSE,
  nsplit = 1,
  level = 0.95,
  area = 1,
  silent = FALSE,
  predict_args = list(),
  derived_link = NULL,
  ...
)

get_cog(
  obj,
  bias_correct = FALSE,
  level = 0.95,
  format = c("long", "wide"),
  area = 1,
  silent = TRUE,
  derived_link = NULL,
  ...
)

get_weighted_average(
  obj,
  vector,
  bias_correct = FALSE,
  level = 0.95,
  area = 1,
  silent = TRUE,
  derived_link = NULL,
  ...
)

get_eao(
  obj,
  bias_correct = FALSE,
  level = 0.95,
  area = 1,
  silent = TRUE,
  derived_link = NULL,
  ...
)
```

## Arguments

- obj:

  Output from
  [`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
  with `return_tmb_object = TRUE` (the most common case). Alternatively,
  if [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)
  was called with `do_index = TRUE`, or if using `get_index_split()`, an
  object from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

- bias_correct:

  Should bias correction be implemented via
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)? Bias
  correction accounts for the non-linear transformation of random
  effects when calculating the index. Recommended to be `TRUE` for final
  analyses, but can be set to `FALSE` for faster calculation while
  experimenting with models. See Thorson and Kristensen (2016) in the
  References.

- level:

  The confidence level.

- area:

  Grid cell area for area weighting the index. Can be: (1) a numeric
  vector of length `nrow(newdata)` with area for each grid cell, (2) a
  single numeric value to apply to all grid cells, or (3) a character
  value giving the column name in `newdata` containing areas. See
  Details for non-spatial uses of `area` as an integration multiplier.

- silent:

  Logical. Suppress progress messages?

- derived_link:

  Optional override for the inverse link used when calculating derived
  quantities such as the index. By default, the fitted family link is
  used. Currently supported for non-delta
  [`binomial()`](https://rdrr.io/r/stats/family.html) and
  [`betabinomial()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  models fit with `link = "cloglog"`.

- ...:

  Passed to
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- newdata:

  New data (e.g., a prediction grid by year) to pass to
  [`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
  in the case of `get_index_split()`.

- nsplit:

  The number of splits to do the calculation in. For memory intensive
  operations (large grids and/or models), it can be helpful to do the
  prediction, area integration, and bias correction on subsets of time
  slices (e.g., years) instead of all at once. If `nsplit > 1`, this
  will usually be slower but with reduced memory use.

- predict_args:

  A list of arguments to pass to
  [`predict.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
  in the case of `get_index_split()`.

- format:

  Long or wide.

- vector:

  A numeric vector of the same length as the prediction data, containing
  the values to be averaged (e.g., depth, temperature).

## Value

For `get_index()`: A data frame with columns for time, estimate
(area-weighted total abundance or biomass), lower and upper confidence
intervals, log estimate, and standard error of the log estimate.

For `get_cog()`: A data frame with columns for time, estimate (center of
gravity: the abundance-weighted mean x and y coordinates), lower and
upper confidence intervals, and standard error of center of gravity
coordinates.

For `get_eao()`: A data frame with columns for time, estimate (effective
area occupied: the area required if the population was spread evenly at
the arithmetic mean density), lower and upper confidence intervals, log
EAO, and standard error of the log EAO estimates.

For `get_weighted_average()`: A data frame with columns for time,
estimate (weighted average of the provided vector, weighted by predicted
density), lower and upper confidence intervals, and standard error of
the estimates.

## Details

More generally, `area` is the multiplier used to integrate
response-scale predictions. For binomial or beta-binomial models fit to
proportions with `weights` specifying the number of trials, predictions
are expected proportions per trial. In that case, `area` can be used as
a standard number of trials (e.g., hooks per longline set) to obtain an
expected-count index. The original fitting `weights` are not
automatically reused for index standardization; supply the desired
standardization multiplier through `area`.

## References

Geostatistical model-based indices of abundance (along with many newer
papers):

Shelton, A.O., Thorson, J.T., Ward, E.J., and Feist, B.E. 2014. Spatial
semiparametric models improve estimates of species abundance and
distribution. Canadian Journal of Fisheries and Aquatic Sciences 71(11):
1655–1666.
[doi:10.1139/cjfas-2013-0508](https://doi.org/10.1139/cjfas-2013-0508)

Thorson, J.T., Shelton, A.O., Ward, E.J., and Skaug, H.J. 2015.
Geostatistical delta-generalized linear mixed models improve precision
for estimated abundance indices for West Coast groundfishes. ICES J.
Mar. Sci. 72(5): 1297–1310.
[doi:10.1093/icesjms/fsu243](https://doi.org/10.1093/icesjms/fsu243)

Geostatistical model-based center of gravity:

Thorson, J.T., Pinsky, M.L., and Ward, E.J. 2016. Model-based inference
for estimating shifts in species distribution, area occupied and centre
of gravity. Methods Ecol Evol 7(8): 990–1002.
[doi:10.1111/2041-210X.12567](https://doi.org/10.1111/2041-210X.12567)

Geostatistical model-based effective area occupied:

Thorson, J.T., Rindorf, A., Gao, J., Hanselman, D.H., and Winker, H.
2016. Density-dependent changes in effective area occupied for
sea-bottom-associated marine fishes. Proceedings of the Royal Society B:
Biological Sciences 283(1840): 20161853.
[doi:10.1098/rspb.2016.1853](https://doi.org/10.1098/rspb.2016.1853)

Bias correction:

Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method
for bias correction in statistical models using random effects, with
spatial and population dynamics examples. Fisheries Research 175: 66–74.
[doi:10.1016/j.fishres.2015.11.016](https://doi.org/10.1016/j.fishres.2015.11.016)

## See also

[`get_index_sims()`](https://sdmTMB.github.io/sdmTMB/reference/get_index_sims.md)

## Examples

``` r
# \donttest{
if (ggplot2_installed()) {
library(ggplot2)

# use a small number of knots for this example to make it fast:
mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 60)

# fit a spatiotemporal model:
m <- sdmTMB(
 data = pcod,
 formula = density ~ 0 + as.factor(year),
 time = "year", mesh = mesh, family = tweedie(link = "log")
)

# prepare a prediction grid:
nd <- replicate_df(qcs_grid, "year", unique(pcod$year))

# Note `return_tmb_object = TRUE` and the prediction grid:
predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)

# biomass index:
ind <- get_index(predictions, bias_correct = TRUE)
ind
ggplot(ind, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  ylim(0, NA)

# do that in 2 chunks
# only necessary for very large grids to save memory
# will be slower but save memory
# note the first argument is the model fit object:
ind <- get_index_split(m, newdata = nd, nsplit = 2, bias_correct = TRUE)

# center of gravity:
cog <- get_cog(predictions, format = "wide")
cog
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_point() +
  geom_linerange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_linerange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()

# effective area occupied:
eao <- get_eao(predictions)
eao
ggplot(eao, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  ylim(0, NA)

# weighted average (e.g., depth-weighted by biomass):
wa <- get_weighted_average(predictions, vector = nd$depth)
wa
ggplot(wa, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4)
}
#> Calculating index in 2 chunks ■■■■■■■■■■■■■■■■                  50% | ETA:  0s
#> Calculating index in 2 chunks ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
#> Bias correction is turned off.
#> It is recommended to turn this on for final inference.
#> Bias correction is turned off.
#> It is recommended to turn this on for final inference.
#> Bias correction is turned off.
#> It is recommended to turn this on for final inference.
#> Bias correction is turned off.
#> It is recommended to turn this on for final inference.

# }
```
