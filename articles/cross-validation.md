# Cross-validation for model evaluation and comparison

**If the code in this vignette has not been evaluated, a rendered
version is available on the [documentation
site](https://sdmTMB.github.io/sdmTMB/index.html) under ‘Articles’.**

``` r
library(ggplot2)
library(dplyr)
library(sdmTMB)
```

## Overview

Cross-validation is one of the best approaches that can be used to
quantify model performance and compare sdmTMB models with different
structures (unlike AIC, this approach will also factor in uncertainty in
random effects). Arguably the most challenging decision in implementing
cross-validation is how to specify the folds (each fold representing a
subset of data that is in turn held out and used as a test set). Folds
may vary in number and how data are partitioned, and will likely be
slightly different for each application.

The goals of some sdmTMB applications may be focused on spatial
prediction; these include making prediction to new spatial regions
(e.g. unsampled areas or areas not sampled in every year). For these
types of models we recommend exploring folds using the `blockCV` or
`spatialsample` packages (Valavi *et al.* 2019; Silge 2021). In general,
these spatial sampling approaches assign observations that are spatially
autocorrelated to the same fold. Accounting for the spatial correlation
can lead to better estimates of covariate effects, as well as prediction
errors.

Alternatively, the goals of an analysis with sdmTMB may be to evaluate
the predictive accuracy of a model in time (e.g. a missing survey year,
or prediction to future years). For retrospective analyses, all points
within a year may be assigned to a fold (or groups of years to the same
fold). In contrast, models that are forward looking would use Leave
Future Out Cross-Validation (LFOCV). In LFOCV, data up to year $t$ are
used to predict observations at $t + 1$, etc.

## Cross validation in sdmTMB

Cross validation in sdmTMB is implemented using the
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
function, with the `k_folds` argument specifying the number of folds
(defaults to 8). The function uses parallelization by default a
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
is set, but this can be turned off with the `parallel` argument.

``` r
data(pcod)
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
pcod$fyear <- as.factor(pcod$year)
```

``` r
# Set parallel processing if desired:
library(future)
plan(multisession, workers = 2)
m_cv <- sdmTMB::sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  k_folds = 4
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
```

In the above example, folds are assigned randomly—but these can be
modified to specific spatial or temporal applications. Without getting
into the complexities of the `blockCV` or `spatialsample` packages, we
could simply use `kmeans` to generate spatial clusters, e.g.

``` r
clust <- kmeans(pcod[, c("X", "Y")], 20)$cluster

m_cv <- sdmTMB::sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  family = tweedie(link = "log")
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
```

Or similarly, these clusters could be assigned in time—here, each year
to a unique fold. Note that year is not included as a factor and
spatiotemporal fields are turned off because they cannot be estimated in
missing years.

``` r
clust <- as.numeric(as.factor(pcod$year))

m_cv <- sdmTMB::sdmTMB_cv(
  density ~ 0 + s(depth_scaled),
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  spatiotemporal = "off",
  family = tweedie(link = "log")
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.
```

## Measuring model performance

Lots of measures of predictive accuracy can be used to evaluate model
performance. By default,
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
returns a list that contains the sum of the log likelihoods for each
left-out fold and the total summed across the left-out folds. This is
roughly equivalent to the expected log predictive density (ELPD) in the
Bayesian literature and can be interpreted as the predictive ability of
the model for new observations. These can be accessed as below, and
inspecting the quantities across folds may help elucidate whether there
are particular folds that are difficult to predict.

``` r
m_cv <- sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  k_folds = 4
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

m_cv$fold_loglik # fold log-likelihood
#> [1] -1640.111 -1732.549 -1722.904 -1507.211
m_cv$sum_loglik # total log-likelihood
#> [1] -6602.776
```

The `fold_loglik` values represent the sum of predictive log-likelihood
values across observations for each fold, and the `sum_loglik` value
represents the sum of `fold_loglik` values, representing the sum for the
entire dataset. We can also calculate derived statistics based on the
predictive values – here, using examples with Root Mean Square Error
(RMSE) and Mean Absolute Error (MAE).

``` r
m <- sdmTMB_cv(
  density ~ depth_scaled + depth_scaled2,
  data = pcod, mesh = make_mesh(pcod, c("X", "Y"), cutoff = 25),
  family = tweedie(link = "log"), k_folds = 3
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# RMSE across entire dataset:
sqrt(mean((m$data$density - m$data$cv_predicted)^2)) 
#> [1] 194.9289
# MAE across entire dataset:
mean(abs(m$data$density - m$data$cv_predicted))
#> [1] 50.53269
```

Alternatively, we might be interested in calculating RMSE and MAE by
fold,

``` r
# RMSE and MAE by fold:
group_by(m$data, cv_fold) |> 
  summarize(
    rmse = sqrt(mean((density - cv_predicted)^2)),
    mae = mean(abs(density - cv_predicted))
  )
#> # A tibble: 3 × 3
#>   cv_fold  rmse   mae
#>     <int> <dbl> <dbl>
#> 1       1  281.  58.2
#> 2       2  150.  48.8
#> 3       3  113.  44.6
```

## Single splits

In cases where only a single test set is evaluated (e.g., 10% of the
data), using the
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
function may be overkill because two
[`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) models
will be fit, but using this function may be worthwhile to reduce coding
errors (in the log-likelihood calculations). For example, here we assign
two folds, randomly holding out 10% of the observations as a test set
(the test set is given ID = 1, and the training set is given ID = 2).

``` r
clust <- sample(1:2, size = nrow(pcod), replace = TRUE, prob = c(0.1, 0.9))

m_cv <- sdmTMB_cv(
  density ~ 0 + s(depth_scaled) + fyear,
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  family = tweedie(link = "log"),
  k_folds = length(unique(clust))
)
```

We can ignore the total log-likelihood, and just focus on the first
element of list list:

``` r
m_cv$fold_loglik[[1]]
#> [1] -771.5139
```

## Comparing two or more models

We can use the output of
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
to compare two or more models. For example, if we wanted to evaluate the
support for a depth effect or not, we could do 10-fold cross validation
(it’s important that the folds be the same across the two models). In
this example, using either the predictive log-likelihood or ELPD would
lead one to conclude that including depth improves the predictive
accuracy of the model.

``` r
clust <- sample(seq_len(10), size = nrow(pcod), replace = TRUE)

m1 <- sdmTMB::sdmTMB_cv(
  density ~ 0 + fyear,
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  family = tweedie(link = "log")
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

m2 <- sdmTMB::sdmTMB_cv(
  density ~ 0 + fyear + s(depth_scaled),
  data = pcod,
  mesh = mesh,
  fold_ids = clust,
  family = tweedie(link = "log")
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# Compare log-likelihoods -- higher is better!
m1$sum_loglik
#> [1] -6728.762
m2$sum_loglik
#> [1] -6554.709
```

## Model ensembling

Finally, instead of identifying single “best” models, we may be
interested in doing model averaging. In the sdmTMB package, we’ve
implemented the model stacking procedure described by (Yao *et al.*
2018) in the
[`sdmTMB_stacking()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_stacking.md)
function. This procedure uses optimization to find the normalized
weights that maximize the total log-likelihood across models (other
metrics may also be used). Inputs to the function are a list of models
(a fictitious `model_list`), where each list element is the output of a
call to
[`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md):

``` r
weights <- sdmTMB_stacking(model_list)
```

By default this calculation uses data from each fold. If instead, we had
split the data into the 10/90 split (as in the example above), we
wouldn’t want to use the 2nd model fit to generate these weights. If we
had just wanted to use the predictions from the first fold onto the 10%
test set, we could specify that using the `include_folds` argument.

``` r
weights <- sdmTMB_stacking(model_list, include_folds = 1)
```

## Calculating measures of predictive skill for binary data

For delta models, or models of presence-absence data, several measures
of predictive ability are available. These are applicable to cross
validation, although we demonstrate them here first in a non-cross
validation context for simplicity.

A first commonly used diagnostic is the AUC (Area Under the Curve),
which quantifies the ability of a model to discriminate between the two
classes; this is done from the Receiver Operating Characteristic (ROC)
curve, which plots the true positive rate vs. false positive rate. There
are several packages to calculate AUC in R, but this can be done with
the `pROC` package, where inputs are a vector of 0s and 1s (or factor
equivalents) in the raw data, and a vector of estimated probabilities
(generated from a call to
[`predict()`](https://rdrr.io/r/stats/predict.html), as shown below).
The [`plogis()`](https://rdrr.io/r/stats/Logistic.html) function is
needed to convert estimated values in logit space to probabilities in
natural (zero to one) space.

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(present ~ s(depth), data = pcod, mesh = mesh)
pred <- predict(fit) # presence-absence model
roc <- pROC::roc(pcod$present, plogis(pred$est))
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
auc <- pROC::auc(roc)
auc
#> Area under the curve: 0.8831
```

With a delta model, two estimated values are returned, so only the first
would be used. E.g.,

``` r
fit <- sdmTMB(density ~ 1, data = pcod, 
  mesh = mesh, family = delta_gamma())
pred <- predict(fit)

# the first linear predictor is the binomial component (est1):
roc <- pROC::roc(pcod$present, plogis(pred$est1))
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases
auc <- pROC::auc(roc)
auc
#> Area under the curve: 0.8604
```

If we wanted to apply this in the context of cross validation, we could
do it like this:

``` r
x <- sdmTMB_cv(
  present ~ s(depth), data = pcod, spatial = "off",
  mesh = mesh, family = binomial(), k_folds = 2
)
roc <- pROC::roc(x$data$present, plogis(x$data$cv_predicted))
auc <- pROC::auc(roc)
auc
```

AUC may be sensitive to imbalances in the data, however, and alternative
metrics may better approximate skill. Here we highlight an example of
using true skill score (implemented in packages such as SDMtune):

``` r
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
fit <- sdmTMB(present ~ 1, data = pcod, 
  mesh = mesh, family = binomial())
```

Next, we can generate predicted probabilities and classes using a
threshold of 0.5 as an example:

``` r
pred <- predict(fit)
pred$p <- plogis(pred$est)
pred$pred_01 <- ifelse(pred$p < 0.5, 0, 1)
```

Next we create a confusion matrix and calculate the true skill score:

``` r
conmat <- table(pred$pred_01, pred$present)
true_neg <- conmat[1, 1]
false_neg <- conmat[1, 2]
false_pos <- conmat[2, 1]
true_pos <- conmat[2, 2]

# Calculate TSS:
true_pos_rate <- true_pos / (true_pos + false_neg)
true_neg_rate <- true_neg / (true_neg + false_pos)
TSS <- true_pos_rate + true_neg_rate - 1
TSS
#> [1] 0.5238745
```

In some cases, reporting the true negative or true positive rate might
be of interest in addition to TSS.

## Multi-scale spatial assessment with waywiser

The waywiser package (Mahoney 2023) provides tools for assessing spatial
models, including methods to evaluate model performance at multiple
spatial scales. A common pattern with spatial models is that predictions
for individual observation units (e.g., points, pixels) may be
aggregated to arbitrary scales (e.g., management units, grid cells of
varying sizes). Because prediction errors can be spatially distributed
and may compound or counteract each other when aggregated, it’s useful
to assess model performance across multiple scales.

The
[`cv_to_waywiser()`](https://sdmTMB.github.io/sdmTMB/reference/cv_to_waywiser.md)
function converts `sdmTMB_cv` objects to sf (simple features) format,
which can then be passed to waywiser functions like `ww_multi_scale()`.

``` r
# First run cross-validation
set.seed(42)
mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
m_cv <- sdmTMB_cv(
  density ~ s(depth_scaled),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  k_folds = 2
)
#> Running fits with `future.apply()`.
#> Set a parallel `future::plan()` to use parallel processing.

# Convert to sf format for spatial assessment
cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"))
#> Detected UTM zone 9N; CRS = 32609.
#> Visit https://epsg.io/32609 to verify.
print(cv_sf)
#> Simple feature collection with 2143 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 343061.7 ymin: 5635893 xmax: 579368.1 ymax: 5839019
#> Projected CRS: WGS 84 / UTM zone 9N
#> First 10 features:
#>        truth     estimate                 geometry
#> 1  113.13848 37.532270279 POINT (446475.2 5793426)
#> 2   41.70492 11.173070839 POINT (446459.4 5800136)
#> 3    0.00000  9.561890006 POINT (448598.7 5801687)
#> 4   15.70614 98.612458614 POINT (436915.7 5802305)
#> 5    0.00000  7.037695083 POINT (420610.1 5771055)
#> 6    0.00000  0.891090606   POINT (417713 5772205)
#> 7    0.00000  0.033984196 POINT (408208.8 5771287)
#> 8    0.00000  0.007256825 POINT (408356.2 5766882)
#> 9    0.00000  0.427827876 POINT (414365.6 5760894)
#> 10  13.94762  6.558124623 POINT (400474.5 5749034)
```

This sf object can now be used with waywiser functions. For example, to
assess model performance at multiple spatial scales:

``` r
# Assess performance at different grid resolutions
# n controls the number of grid cells in x and y directions
multi_scale_results <- waywiser::ww_multi_scale(
  cv_sf,
  truth,    # column name (unquoted)
  estimate, # column name (unquoted)
  # 10x10, 5x5, and 2x2 grids:
  n = list(c(10, 10), c(5, 5), c(2, 2))
)

multi_scale_results
#> # A tibble: 6 × 6
#>   .metric .estimator .estimate .grid_args       .grid          .notes          
#>   <chr>   <chr>          <dbl> <list>           <list>         <list>          
#> 1 rmse    standard       19.4  <tibble [1 × 1]> <sf [100 × 5]> <tibble [0 × 2]>
#> 2 mae     standard       12.2  <tibble [1 × 1]> <sf [100 × 5]> <tibble [0 × 2]>
#> 3 rmse    standard        7.23 <tibble [1 × 1]> <sf [25 × 5]>  <tibble [0 × 2]>
#> 4 mae     standard        4.80 <tibble [1 × 1]> <sf [25 × 5]>  <tibble [0 × 2]>
#> 5 rmse    standard        3.45 <tibble [1 × 1]> <sf [4 × 5]>   <tibble [0 × 2]>
#> 6 mae     standard        3.40 <tibble [1 × 1]> <sf [4 × 5]>   <tibble [0 × 2]>
```

The results show how metrics like RMSE and MAE change as predictions are
aggregated to coarser spatial scales, which can reveal whether your
model’s performance is scale-dependent.

See the waywiser documentation
([`?waywiser::ww_multi_scale`](https://docs.ropensci.org/waywiser/reference/ww_multi_scale.html)
and
[`vignette("multi-scale-assessment", package = "waywiser")`](https://docs.ropensci.org/waywiser/articles/multi-scale-assessment.html))
for more details on interpreting multi-scale assessment results and
additional functions for spatial model evaluation.

## References

Mahoney, M. (2023). *Waywiser: Ergonomic methods for assessing spatial
models*. Retrieved from <https://CRAN.R-project.org/package=waywiser>

Silge, J. (2021). *Spatialsample: Spatial resampling infrastructure*.
Retrieved from <https://CRAN.R-project.org/package=spatialsample>

Valavi, R., Elith, J., Lahoz-Monfort, J.J. & Guillera-Arroita, G.
(2019). blockCV: An r package for generating spatially or
environmentally separated folds for k-fold cross-validation of species
distribution models. *Methods in Ecology and Evolution*, **10**,
225–232.

Yao, Y., Vehtari, A., Simpson, D. & Gelman, A. (2018). Using Stacking to
Average Bayesian Predictive Distributions (with Discussion). *Bayesian
Analysis*, **13**, 917–1007.
