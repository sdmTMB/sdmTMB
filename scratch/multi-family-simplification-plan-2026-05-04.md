# Multi-family maintainability refactor plan

Date: 2026-05-04

## Goal

Make multi-family support easier to maintain by replacing scattered `multi_family`
branches with one normalized internal observation-family specification that can
represent:

- ordinary single-family models,
- ordinary delta models,
- multi-family non-delta models,
- multi-family models with delta families.

The public API should stay unchanged. The first implementation target is internal
clarity and test-equivalent behavior, not new user-facing functionality.

## Current problem

The current implementation has two internal representations:

- Single-family and regular delta models use `family`, `link`, `ln_phi`, `thetaf`,
  `ln_student_df`, `gengamma_Q`, and `poisson_link_delta`.
- Multi-family models add separate row-level and family-level metadata:
  `multi_family`, `e_i`, `e_g`, `family_e1`, `family_e2`, `link_e1`, `link_e2`,
  `delta_family_e`, `poisson_link_delta_e`, and parameter offset vectors.

This creates repeated routing logic in several places:

- `R/fit.R`: TMB data, parameter, and map construction.
- `R/predict.R`: prediction assembly, component transforms, combined delta output,
  population prediction, and standard-error handling.
- `src/sdmTMB.cpp`: repeated per-observation family/link/delta/parameter resolution.
- `R/tmb-sim.R`, `R/residuals.R`, and index helpers: feature-specific guards and
  special cases.

The maintainability issue is not the likelihood switch itself. The issue is the
repeated logic that decides which family, component, link, and parameter slot
applies to a row.

## Desired internal model

Introduce a normalized `family_spec` used by all observation models.

Conceptually:

```text
family_table
  family_id                 integer, 1-based in R; 0-based in TMB
  name                      user-facing family-list name, or default
  family1                   family enum for component 1
  link1                     link enum for component 1
  family2                   family enum for component 2, or -1
  link2                     link enum for component 2, or -1
  is_delta                  logical/integer
  is_poisson_link_delta     logical/integer
  ln_phi_slot               parameter slot, or -1
  thetaf_slot               parameter slot, or -1
  ln_student_df_slot        parameter slot, or -1
  gengamma_Q_slot           parameter slot, or -1
  family_object             R family object, for R-side transforms only

obs_family_id
  integer vector of length n_obs

pred_family_id
  integer vector of length n_pred
```

Representations:

- Single-family model: one `family_table` row; all `obs_family_id = 1`.
- Regular delta model: one `family_table` row with `is_delta = TRUE`.
- Multi-family model: one `family_table` row per named family; `obs_family_id`
  comes from `distribution_column`.
- Multi-family with delta rows: rows that correspond to delta family objects have
  `is_delta = TRUE` and valid component-2 metadata.

## Non-goals

- Do not change the user-facing `family = list(...)` API.
- Do not add support for `dispformula` in multi-family models as part of this
  refactor.
- Do not add support for distributed lags plus multi-family models as part of this
  refactor.
- Do not rewrite the full C++ likelihood switch unless a later step proves it is
  necessary.
- Do not remove existing unsupported-case guards until equivalent behavior is
  covered by tests.

## Step 1: Add a normalized R family spec builder

Create or extend helpers in `R/multi-likelihood.R`.

Proposed helpers:

```r
.build_family_spec <- function(family, data = NULL, distribution_column = NULL)
.family_spec_for_prediction <- function(object, newdata)
.family_spec_param_offsets <- function(family_table)
.validate_family_spec <- function(spec, where)
```

Behavior:

1. Detect whether `family` is a named list or a standard family object.
2. For non-list families, create a one-row `family_table`.
3. For regular delta families, create a one-row `family_table` with both
   components populated.
4. For named-list families, preserve the existing validation rules:
   - names required,
   - no duplicate names,
   - all elements inherit from `"family"`,
   - no `_mix` families,
   - max two components,
   - second component only for delta families,
   - delta first component must be binomial,
   - current allowed-family list unchanged.
5. Build `obs_family_id` from `distribution_column` when needed.
6. Store both 1-based R ids and 0-based TMB ids deliberately:
   - `obs_family_id` for R,
   - `obs_family_id0` for TMB.
7. Replace or wrap `.parse_multi_family()` so existing call sites can migrate
   gradually.

Acceptance criteria:

- Existing multi-family parser tests still pass.
- Add new tests proving `.build_family_spec()` produces equivalent metadata for:
  - single-family Gaussian,
  - regular delta Gamma,
  - multi-family Gaussian/Gamma/Tweedie,
  - multi-family with one delta family and one non-delta family.

## Step 2: Move parameter-slot logic onto `family_table`

Replace `.multi_family_param_offsets()` with a generalized helper that operates
on the normalized table.

Rules:

1. Determine the target distribution family for parameter slots:
   - non-delta: component 1,
   - delta: component 2.
2. Assign slots for:
   - `ln_phi`,
   - `thetaf`,
   - `ln_student_df`,
   - `gengamma_Q`.
3. Store slot ids directly on `family_table`.
4. Preserve the current behavior that fixed student df is unsupported for
   multi-family models.

Important design choice:

- Keep single-family TMB parameters unchanged at first.
- For multi-family, continue using compact vectors such as `ln_phi_e`.
- The normalized table should hide this difference from R-side prediction and
  C++ routing.

Acceptance criteria:

- Existing parameter-map tests pass.
- Add tests that inspect slot assignment for families with and without dispersion
  parameters.
- Add a test where only one of several multi-family rows uses Tweedie and confirm
  exactly one `thetaf_e` slot is estimated.

## Step 3: Use `family_spec` in TMB data construction

Update `R/fit.R` to build `family_spec` early and use it when constructing
`tmb_data`.

Initial compatibility mapping:

```r
tmb_data$family_table_family1
tmb_data$family_table_family2
tmb_data$family_table_link1
tmb_data$family_table_link2
tmb_data$family_table_is_delta
tmb_data$family_table_is_poisson_link_delta
tmb_data$family_table_ln_phi_slot
tmb_data$family_table_thetaf_slot
tmb_data$family_table_ln_student_df_slot
tmb_data$family_table_gengamma_Q_slot
tmb_data$obs_family_id
tmb_data$pred_family_id
```

Keep the old fields temporarily:

- `multi_family`
- `e_i`
- `e_g`
- `family_e1`
- `family_e2`
- `link_e1`
- `link_e2`
- `delta_family_e`
- `poisson_link_delta_e`
- parameter start/length vectors

This makes the first pass non-disruptive. The old fields can be removed after
C++ and prediction use the new fields.

Acceptance criteria:

- Fitted objects contain `family_spec`.
- `tmb_data` contains both legacy and new fields.
- Existing tests still pass.
- New tests compare legacy fields to the new table for representative models.

## Step 4: Centralize R-side family transforms

Add helpers in `R/multi-likelihood.R` for prediction-scale transformations.

Proposed helpers:

```r
.family_spec_response_transform <- function(eta1, eta2 = NULL, family_id, spec)
.family_spec_link_transform <- function(eta1, eta2 = NULL, family_id, spec)
.family_spec_component_estimates <- function(eta1, eta2 = NULL, family_id, spec, type)
```

Responsibilities:

1. Transform component predictions to response scale.
2. Combine standard delta predictions.
3. Combine Poisson-link delta predictions.
4. Return a consistent list:

```r
list(
  est = ...,
  est1 = ...,
  est2 = ...,
  component = ...
)
```

Migration:

- Reimplement `.multi_family_predict_est()` as a wrapper around these helpers.
- Add equivalent helpers for regular delta models so `predict()` does not need
  separate regular-delta and multi-family transform logic.

Acceptance criteria:

- Unit tests compare old `.multi_family_predict_est()` behavior to the new helper.
- Add direct tests for:
  - non-delta response transform,
  - standard delta link and response transforms,
  - Poisson-link delta link and response transforms,
  - mixed family ids in the same prediction vector.

## Step 5: Simplify prediction assembly

Refactor `R/predict.R` in a narrow pass.

Add a helper:

```r
.assemble_prediction_output <- function(
  nd,
  report,
  family_spec,
  pred_family_id,
  type,
  model,
  pop_pred,
  se_fit,
  sdreport_est = NULL,
  sdreport_se = NULL,
  object
)
```

Responsibilities:

1. Select `proj_eta` or `proj_fe` depending on `pop_pred`.
2. Apply component selection when `model` is `1` or `2`.
3. Apply combined prediction logic when `model = NA`.
4. Add `est`, `est1`, `est2`, and `est_se` consistently.
5. Leave random-field columns (`omega_s`, `epsilon_st`, `zeta_s_*`) in a smaller
   helper so the statistical prediction logic is separate from output decoration.

Migration order:

1. Route only multi-family non-`se_fit` predictions through the new helper.
2. Route regular non-delta predictions through it.
3. Route regular delta predictions through it.
4. Route `se_fit` predictions last.

Acceptance criteria:

- Existing prediction tests pass after each migration step.
- Add snapshot-like tests for column names in:
  - single-family prediction,
  - regular delta prediction,
  - multi-family non-delta prediction,
  - multi-family mixed delta/non-delta prediction.
- Confirm current unsupported cases still error with the same or clearer messages:
  - multi-family `type = "response", se_fit = TRUE`,
  - simulated combined multi-family delta predictions,
  - simulated response-scale multi-family delta predictions.

## Step 6: Add C++ family-resolution helpers

Create a small C++ helper, ideally in a new header such as
`src/observation-family.h`.

Proposed struct:

```cpp
template <class Type>
struct ObservationFamilyResolved {
  bool active;
  int family;
  int link;
  bool is_delta;
  bool is_poisson_link_delta;
  Type ln_phi;
  Type phi;
  Type thetaf;
  Type ln_student_df;
  Type gengamma_Q;
};
```

Proposed resolver:

```cpp
template <class Type>
ObservationFamilyResolved<Type> resolve_observation_family(
  int obs_i,
  int model_m,
  const ObservationFamilyContext<Type>& ctx
);
```

Responsibilities:

1. Determine `family_id` for observation `i`.
2. Determine whether component `m` is active.
3. Resolve family and link for component `m`.
4. Resolve delta and Poisson-link-delta flags.
5. Resolve distribution parameter values:
   - use single-family parameters for single-family specs,
   - use table slots for multi-family specs.

Migration:

- First implement resolver using both old and new TMB fields if needed.
- Replace repeated C++ blocks in:
  - linear predictor/mu construction,
  - simulated observation mean replacement,
  - likelihood evaluation.

Acceptance criteria:

- The C++ code no longer repeats the same `if (multi_family) { ... }` family/link
  selection block in multiple loops.
- Existing fit, simulation, and prediction tests pass.
- Add a low-level test if feasible through fitted reports: same likelihood before
  and after the resolver refactor for representative models.

## Step 7: Remove legacy multi-family TMB fields

Only do this after Steps 3-6 are stable.

Remove or stop populating:

- `e_i`
- `e_g`
- `family_e1`
- `family_e2`
- `link_e1`
- `link_e2`
- `delta_family_e`
- `poisson_link_delta_e`
- `ln_phi_start`
- `ln_phi_len`
- `thetaf_start`
- `thetaf_len`
- `ln_student_df_start`
- `ln_student_df_len`
- `gengamma_Q_start`
- `gengamma_Q_len`

Keep `multi_family` only if it remains useful as a user-facing model-type flag.
Prefer deriving behavior from `family_table` and `family_spec`.

Acceptance criteria:

- No C++ code depends on legacy fields.
- No R code depends on legacy fields except compatibility methods for old fitted
  objects, if needed.
- Tests inspect `family_spec` instead of legacy vectors.

## Step 8: Preserve compatibility with old fitted objects

Add a compatibility helper:

```r
.family_spec_from_object <- function(object)
```

Behavior:

1. If `object$family_spec` exists, validate and return it.
2. Else, reconstruct from legacy fields:
   - `object$family`,
   - `object$tmb_data$multi_family`,
   - `object$tmb_data$e_i`,
   - `object$tmb_data$family_e1`,
   - `object$tmb_data$family_e2`,
   - `object$tmb_data$link_e1`,
   - `object$tmb_data$link_e2`,
   - `object$tmb_data$delta_family_e`,
   - `object$tmb_data$poisson_link_delta_e`.
3. Emit no warning for normal prediction from old objects unless reconstruction is
   incomplete.

Acceptance criteria:

- Existing saved-object tests, if any, continue to work.
- `predict()` on old multi-family objects either works or errors with a clear
  "please refit" message only when essential metadata is missing.

## Step 9: Revisit unsupported combinations deliberately

After the representation is normalized, review unsupported cases one at a time.

Keep unsupported unless there is a clear implementation path:

- `dispformula` with multi-family models,
- distributed lags with multi-family models,
- response-scale `se_fit` for multi-family models,
- simulated combined predictions for multi-family delta models,
- index calculations for mixed delta/non-delta family rows if semantics are
  unclear.

For each unsupported case, add a test that asserts the exact guard. This prevents
future refactors from accidentally exposing partially implemented behavior.

## Step 10: Documentation and developer notes

Update internal documentation after the refactor lands.

Recommended files:

- `scratch/multi-family-maintainability-assessment-2026-02-05.md`, if still
  relevant.
- A new short developer note in `R/multi-likelihood.R` near the spec builder.
- `NEWS.md`, only if user-visible behavior changes.

Document:

- the normalized family table,
- how observation rows map to family rows,
- component semantics for delta rows,
- which unsupported cases are intentional.

## Suggested implementation sequence

Use small PR-sized commits:

1. Add `family_spec` builder and tests; no call-site migration.
2. Add parameter-slot table and tests; no C++ migration.
3. Store `family_spec` on fitted objects and duplicate legacy TMB fields.
4. Add R transform helpers and rewrap `.multi_family_predict_est()`.
5. Refactor `predict()` multi-family non-`se_fit` output assembly.
6. Refactor regular prediction paths to use the same assembly helper.
7. Add C++ resolver helper and migrate one loop at a time.
8. Remove legacy TMB fields after tests are green.
9. Add compatibility fallback for old fitted objects.
10. Run full reverse-style checks that exercise predictions, simulations,
    residuals, and index calculations.

## Test checklist

Run focused tests frequently:

```r
testthat::test_local(filter = "multi-family-core")
testthat::test_local(filter = "multi-family-fit")
testthat::test_local(filter = "multi-family-output")
testthat::test_local(filter = "8-delta2")
testthat::test_local(filter = "prediction")
testthat::test_local(filter = "sdmTMB-simulate|tmb-simulation")
```

Run broader checks before removing legacy fields:

```r
testthat::test_local()
```

Also compare representative fitted likelihoods before and after the C++ resolver
migration. The refactor should be behavior-preserving unless a test explicitly
documents a corrected bug.

## Risk controls

- Keep old and new metadata side by side during migration.
- Do not change public API or error semantics until the internal representation is
  stable.
- Migrate prediction assembly before likelihood code; it is easier to test and
  carries less numerical risk.
- Avoid changing the likelihood switch and the family-routing abstraction in the
  same commit.
- Preserve explicit guards for unsupported combinations.

## Success criteria

The refactor is successful when:

- `predict.R` no longer has separate large branches for multi-family, regular
  delta, and ordinary models where a shared family-spec helper can decide the
  behavior.
- `src/sdmTMB.cpp` resolves active family/link/parameter metadata through one
  helper rather than repeating the same branch in multiple loops.
- `R/fit.R` builds one normalized family representation for all observation model
  types.
- Existing tests pass.
- New tests cover family-table construction, parameter-slot assignment, prediction
  transforms, and unsupported-case guards.
