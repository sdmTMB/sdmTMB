# Package index

## Mesh construction

Tools for converting coordinates into UTMs, constructing SPDE meshes
prior to model fitting, and adding correlation barriers.

- [`add_utm_columns()`](https://sdmTMB.github.io/sdmTMB/reference/add_utm_columns.md)
  [`get_crs()`](https://sdmTMB.github.io/sdmTMB/reference/add_utm_columns.md)
  : Add UTM coordinates to a data frame
- [`make_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md)
  [`plot(`*`<sdmTMBmesh>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/make_mesh.md)
  : Construct an SPDE mesh for sdmTMB
- [`add_barrier_mesh()`](https://sdmTMB.github.io/sdmTMB/reference/add_barrier_mesh.md)
  : Transform a mesh object into a mesh with correlation barriers

## Fitting and predicting

Core tools for model fitting, prediction, and inspection.

- [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md) :
  Fit a spatial or spatiotemporal GLMM with TMB

- [`sanity()`](https://sdmTMB.github.io/sdmTMB/reference/sanity.md) :
  Sanity check of an sdmTMB model

- [`tidy(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/tidy.sdmTMB.md)
  [`tidy(`*`<sdmTMB_cv>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/tidy.sdmTMB.md)
  : Turn sdmTMB model output into a tidy data frame

- [`predict(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/predict.sdmTMB.md)
  : Predict from an sdmTMB model

- [`residuals(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/residuals.sdmTMB.md)
  : Residuals method for sdmTMB models

- [`update(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/update.sdmTMB.md)
  : Update an sdmTMB model

- [`dharma_residuals()`](https://sdmTMB.github.io/sdmTMB/reference/dharma_residuals.md)
  : DHARMa residuals

- [`sdmTMBcontrol()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMBcontrol.md)
  : Optimization control options

- [`run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.md)
  **\[experimental\]** : Run extra optimization on an already fitted
  object

- [`replicate_df()`](https://sdmTMB.github.io/sdmTMB/reference/replicate_df.md)
  : Replicate a prediction data frame over time

- [`set_delta_model()`](https://sdmTMB.github.io/sdmTMB/reference/set_delta_model.md)
  :

  Set delta model for
  [`ggeffects::ggpredict()`](https://strengejacke.github.io/ggeffects/reference/ggpredict.html)

- [`make_category_svc()`](https://sdmTMB.github.io/sdmTMB/reference/make_category_svc.md)
  : Set up spatially varying coefficients for category composition
  models

- [`coef(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/coef.sdmTMB.md)
  : Get fixed-effect coefficients

- [`sigma(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/sigma.sdmTMB.md)
  : Extract residual standard deviation or dispersion parameter

- [`cAIC()`](https://sdmTMB.github.io/sdmTMB/reference/cAIC.md) :
  Calculate conditional AIC

## Families

Additional families beyond the standard R families.

- [`Beta()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`gamma_mix()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`lognormal_mix()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`nbinom2_mix()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`nbinom1()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`truncated_nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`truncated_nbinom1()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`student()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`tweedie()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`censored_poisson()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_gamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_gamma_mix()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_gengamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_lognormal_mix()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_truncated_nbinom2()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_truncated_nbinom1()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_poisson_link_gamma()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_poisson_link_lognormal()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`betabinomial()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  [`delta_beta()`](https://sdmTMB.github.io/sdmTMB/reference/families.md)
  : Additional families

## Priors

Optional priors or penalties on parameters.

- [`sdmTMBpriors()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  [`normal()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  [`halfnormal()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  [`gamma_cv()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  [`mvnormal()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  [`pc_matern()`](https://sdmTMB.github.io/sdmTMB/reference/priors.md)
  **\[experimental\]** : Prior distributions
- [`plot_pc_matern()`](https://sdmTMB.github.io/sdmTMB/reference/plot_pc_matern.md)
  : Plot PC Matérn priors

## Simulation

Simulating new data with an sdmTMB model.

- [`simulate(`*`<sdmTMB>`*`)`](https://sdmTMB.github.io/sdmTMB/reference/simulate.sdmTMB.md)
  : Simulate from a fitted sdmTMB model

- [`sdmTMB_simulate()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_simulate.md)
  : Simulate from a spatial/spatiotemporal model

- [`project()`](https://sdmTMB.github.io/sdmTMB/reference/project.md)
  **\[experimental\]** :

  Project from an sdmTMB model using simulation

## Plotting

Functions for plotting.

- [`plot_anisotropy()`](https://sdmTMB.github.io/sdmTMB/reference/plot_anisotropy.md)
  [`plot_anisotropy2()`](https://sdmTMB.github.io/sdmTMB/reference/plot_anisotropy.md)
  : Plot anisotropy from an sdmTMB model

- [`plot_smooth()`](https://sdmTMB.github.io/sdmTMB/reference/plot_smooth.md)
  : Plot a smooth term from an sdmTMB model

- [`plot_pc_matern()`](https://sdmTMB.github.io/sdmTMB/reference/plot_pc_matern.md)
  : Plot PC Matérn priors

- [`visreg_delta()`](https://sdmTMB.github.io/sdmTMB/reference/visreg_delta.md)
  [`visreg2d_delta()`](https://sdmTMB.github.io/sdmTMB/reference/visreg_delta.md)
  :

  Plot sdmTMB models with the visreg package

## Cross validation

Functions related to cross validation.

- [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
  : Cross validation with sdmTMB models

- [`sdmTMB_stacking()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_stacking.md)
  **\[experimental\]** :

  Perform stacking with log scores on
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
  output

- [`cv_to_waywiser()`](https://sdmTMB.github.io/sdmTMB/reference/cv_to_waywiser.md)
  **\[experimental\]** :

  Convert
  [`sdmTMB_cv()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB_cv.md)
  objects to sf format for spatial assessment with waywiser

## Derived quantities

Derived quantities that can be calculated from sdmTMB models.

- [`get_index()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  [`get_index_split()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  [`get_cog()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  [`get_weighted_average()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  [`get_eao()`](https://sdmTMB.github.io/sdmTMB/reference/get_index.md)
  : Extract a relative biomass/abundance index, center of gravity,
  effective area occupied, or weighted average
- [`get_index_sims()`](https://sdmTMB.github.io/sdmTMB/reference/get_index_sims.md)
  **\[experimental\]** : Calculate a population index via simulation
  from the joint precision matrix
- [`get_range_edge()`](https://sdmTMB.github.io/sdmTMB/reference/get_range_edge.md)
  **\[experimental\]** : Calculate range edges via simulation from the
  joint precision matrix

## Miscellaneous parameter extraction

Functions for calculating effect sizes or extracting parameter samples.

- [`Effect.sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/Effect.sdmTMB.md)
  : Calculate effects

- [`emmeans.sdmTMB`](https://sdmTMB.github.io/sdmTMB/reference/emmeans.sdmTMB.md)
  :

  Estimated marginal means with the emmeans package with sdmTMB

- [`spread_sims()`](https://sdmTMB.github.io/sdmTMB/reference/gather_sims.md)
  [`gather_sims()`](https://sdmTMB.github.io/sdmTMB/reference/gather_sims.md)
  : Extract parameter simulations from the joint precision matrix
