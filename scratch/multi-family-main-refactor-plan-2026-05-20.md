# Multi-Family Refactor Plan on `main`

## Summary

Implement multi-family support on top of `main` by introducing one canonical internal family specification that covers both ordinary single-family fits and row-wise multi-family fits. `n_m` keeps its existing meaning: number of linear predictors/components, capped at 2. Multi-family adds a separate row-to-family mapping axis; it does not change component semantics.

Public API stays aligned with the `dev` vignette: multi-family fits use a named `family` list plus `distribution_column`. Existing single-family behavior on `main` must remain numerically unchanged. Saved `dev`-branch multi-family fits will not be supported; they must be refit. First pass preserves current unsupported combinations and method limits.

## Implementation Changes

### 1. Add one canonical internal family spec in R

Create a single R-side normalization layer, used for **all** fits, that converts the user input into `family_spec` before any formula/TMB setup:

- `n_f`: number of observation families in the fit.
- `n_m`: number of linear predictors/components, equal to `max(components_per_family)`, with valid values `1` or `2`.
- `family_list`: normalized family entries, even for ordinary single-family fits.
- `distribution_column`: stored for multi-family fits, `NULL` otherwise.
- `family_id_i`: 1-based row-to-family mapping in R.
- `active[f, m]`: logical matrix; whether family `f` contributes to component `m`.
- `family_code[f, m]` and `link_code[f, m]`: TMB enums for each active component.
- `combine_kind[f]`: one of `single`, `delta`, or `poisson_link_delta`.
- `param_slot` vectors for `ln_phi`, `thetaf`, `ln_student_df`, and `gengamma_Q`.

Rules:

- An ordinary single-family fit becomes `n_f = 1`, `family_id_i = 1` for all rows.
- A model has `n_m = 2` if any family entry is a 2-component family.
- Reject family entries with more than 2 components.
- Keep formulas, spatial settings, spatiotemporal settings, and `share_range` **component-based**, not family-based.
- If `n_m = 2` and the user supplies a single formula or scalar component settings, replicate current delta behavior.
- LP1-only families ignore LP2; they do not create family-specific formulas.

Response construction must also be unified here:

- LP1-only rows: `y_i[,1] = processed response`, `y_i[,2] = NA`.
- Standard delta rows: `y_i[,1] = encounter`, `y_i[,2] = positive response or NA`.
- Poisson-link delta rows keep the same two-component response shape as current delta code.
- Binomial and betabinomial preprocessing must be row-family aware, including `weights`/`size`.

Store `object$family_spec` on the fitted object for all fits. Keep `object$family` user-facing and backwards-compatible with current `main` behavior.

### 2. Unify TMB family metadata and parameter slots

Change TMB data assembly in [R/fit.R](/Users/seananderson/src/sdmTMB/R/fit.R) so TMB always receives one per-family table, for both single-family and multi-family fits. Do not keep a separate single-family-vs-multi-family representation inside TMB.

Use one 0-based row-family id at the TMB boundary and one set of per-family vectors/matrices:

- `obs_family_id[i]`: 0-based family id for row `i`.
- `component_active[f, m]`: integer active flag.
- `family_code[f, m]`, `link_code[f, m]`.
- `combine_kind[f]`.
- `ln_phi_slot[f]`, `thetaf_slot[f]`, `ln_student_df_slot[f]`, `gengamma_Q_slot[f]`, with `-1` meaning “not used”.

Parameter storage should also be unified:

- Replace the current dual parameter-source approach with per-family slot vectors for all fits.
- Single-family fits simply have length-0 or length-1 slot vectors as appropriate.
- Keep single-family `dispformula` behavior working by treating observation-level dispersion as a row-level override on top of the resolved family slot, but do not support `dispformula` for multi-family in this pass.

Explicit offset rules:

- LP1-only families: apply offset to LP1.
- Standard 2-component families: apply offset to LP2 only, matching current delta behavior.
- Poisson-link delta families: preserve current special handling and do not add the raw offset through the standard branch.

### 3. Replace repeated C++ family branching with one resolver

Refactor [src/sdmTMB.cpp](/Users/seananderson/src/sdmTMB/src/sdmTMB.cpp) to resolve row/component family behavior through one stack-only helper that is safe inside `PARALLEL_REGION`. The resolver should return, for row `i` and component `m`:

- `active`
- `family_code`
- `link_code`
- `combine_kind`
- `offset_applies`
- resolved parameter values for `phi`, Tweedie `p`, student `df`, and generalized-gamma `Q`

This resolver must be the only place that decides whether a row participates in LP1 only or LP1+LP2. All row/component loops should then use `if (!resolved.active) continue;` rather than checking “regular delta” vs “multi-family”.

Migrate all family-sensitive C++ paths together to the resolver:

- LP and `mu` construction
- simulated-observation replacement path
- observation likelihood evaluation
- projection combined-SE reporting
- combined index / weighted-average / EAO aggregation

Combined prediction semantics must be generic and row-family driven:

- `single`: combined response is LP1 inverse-link; combined link is LP1 link space.
- `delta`: combined response is LP1 response times LP2 response; combined link is LP2 link of the product.
- `poisson_link_delta`: preserve current additive/log-space combined logic.

Also replace delta-specific combined report names with generic combined report variables used by all two-component fits, so downstream code extracts one combined prediction path rather than branching by model type.

### 4. Generalize downstream R methods through `family_spec`

Refactor [R/predict.R](/Users/seananderson/src/sdmTMB/R/predict.R) and the small helper layer around it so ordinary single-family fits and multi-family fits both use the same family-spec-driven code.

Prediction rules:

- `newdata` for multi-family must include `distribution_column`, mapped through `family_spec`.
- `predict(..., model = NA)` returns rowwise combined predictions.
- If `n_m = 2`, ordinary predictions also return `est1` and `est2`; `est2` is `NA` for LP1-only rows.
- `predict(..., model = 1)` returns LP1 for all rows.
- `predict(..., model = 2)` returns LP2 where active and `NA` for LP1-only rows.
- `type = "response", se_fit = TRUE` remains unsupported for multi-family in this pass.
- `nsim`/`mcmc_samples` prediction paths must use the same rowwise combiner as ordinary predictions.

Other downstream behavior:

- `print()` and `tidy()` should report per-family summaries and per-family extra parameters using `family_spec`.
- `simulate()` should support combined-response multi-family simulations and use the same rowwise combination logic.
- Standard `residuals()` remain explicitly unsupported for multi-family fits in this pass; keep a clear error recommending simulation-based residual workflows.
- `project()` remains unchanged and unsupported for multi-family fits in this pass.
- Keep single-family fit objects from `main` working normally.
- Do not add a reconstruction layer for saved `dev`-branch multi-family objects; emit a clear “refit required” error if such objects are encountered.

### 5. Cleanup and guardrails

After the unified path is in place:

- Remove redundant `regular_delta` vs `multi_family` routing that is no longer needed.
- Keep only small helpers that answer:
  - how many components the fit has,
  - whether a row-family map is required,
  - how to combine LP1/LP2 for a given family kind.
- Ensure every unsupported combination fails early in R setup, before TMB object construction.

## Test Plan

Add focused equivalence and behavior tests that make the refactor safe:

1. Single-family ordinary fit equals current `main` behavior for Gaussian, Poisson, and one delta family.
2. `family = list(one = gaussian())` matches `family = gaussian()` for log-likelihood, coefficients, and link/response predictions.
3. `family = list(one = delta_gamma())` matches `family = delta_gamma()` for `est`, `est1`, `est2`, and combined link-space behavior.
4. Mixed LP1-only families route row likelihoods, links, and response transforms correctly.
5. Mixed LP1-only plus LP2-active families route rows correctly in fit, prediction, simulation, and index calculations.
6. Offset semantics are correct for LP1-only rows, standard two-component rows, and Poisson-link delta rows.
7. Per-family parameter slotting works for `phi`, Tweedie, student, and generalized-gamma families.
8. `predict(..., model = 1/2/NA)` on mixed fits returns `NA` only where the requested component is inactive.
9. `predict(..., se_fit = TRUE, type = "link")` works for multi-family; `type = "response"` remains blocked with a stable message.
10. `simulate()` returns combined responses for mixed fits and preserves integer-like behavior for count families.
11. `residuals()` and `project()` error clearly for multi-family fits.
12. Unsupported combinations (`dispformula`, `covariate_diffusion`, `_mix()` families) fail early with stable errors.

Add numerical fingerprint tests and run them after each implementation phase:

- single Gaussian
- single standard delta
- single Poisson-link delta
- mixed non-delta
- mixed delta/non-delta
- single-family `dispformula`

Compare objective value, selected fixed effects, selected distribution parameters, a small prediction vector, and index totals before and after the refactor.

## Assumptions and Defaults

- Public multi-family API matches the `dev` vignette: named `family` list plus `distribution_column`.
- `n_m` continues to mean “number of linear predictors/components”, not number of families.
- Multi-family introduces row-wise observation likelihood selection only; it does not introduce family-specific formulas.
- Saved `dev`-branch multi-family fit objects are not supported; they must be refit.
- First pass preserves unsupported combinations for multi-family: `dispformula`, `covariate_diffusion`, `_mix()` families, standard residuals, and projection helpers.
- The implementation target is maintainability over minimal diff size: TMB should have one canonical family resolution path for all fits, not parallel single-family and multi-family paths.
