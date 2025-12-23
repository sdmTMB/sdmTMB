---
name: multi-family-sdmtmb
description: Plan for adding multi-observation likelihoods to sdmTMB
---

# Plan

Implement multi-observation likelihoods in sdmTMB using a per-row family selector while preserving the existing single-family and delta-family paths. Keep the implementation simple by reusing the existing likelihood code paths and only adding minimal per-family parameter plumbing.

## Requirements
- Support named list of families plus a `distribution_column` that maps each row to a family.
- Include delta families; exclude *_mix families from multi-family mode.
- Support families: gaussian, poisson, binomial, betabinomial, Beta, nbinom1/2, Gamma, lognormal, tweedie, student, generalized gamma, and delta variants.
- For multi-family fits, disable regular residuals and require simulation/DHARMa residuals.
- Preserve existing behavior for single-family and delta-family fits without `distribution_column`.
- Allow mixing delta and non-delta families in the same fit.

## Scope
- In: per-row family/link selection, per-family extra parameters where required, predictions that use `distribution_column`.
- Out: *_mix families in multi-family mode; regular residuals for multi-family mode.

## Delta families in multi-family mode

Delta families are treated as another option in the named family list. Key design:

- Each row has one `e_i` value mapping to a family name in the list
- Rows assigned to a **delta family** (e.g., `delta_gamma`) contribute **two** likelihood terms:
  1. Binomial likelihood for presence/absence (using formula1)
  2. Conditional positive likelihood (e.g., Gamma) for positive values only (using formula2; zeroes skipped)
- Rows assigned to **non-delta families** contribute **one** likelihood term
- Delta models retain their two-part formula structure (formula1 for binomial, formula2 for positive)
- Non-delta families use only the main formula (formula1); formula2 is ignored for non-delta families in multi-family mode
- In C++: when `e_i(i)` maps to a delta family, evaluate both binomial and positive likelihoods; otherwise evaluate single likelihood
- Predictions: delta families return predictions from both parts (as currently); non-delta families return single predictions

**Implications:**
- Model fitting: some observations contribute 1 likelihood term, some contribute 2
- Parameter counts: delta families have parameters for both parts
- Residuals: delta families have residuals for both parts (though standard residuals blocked in multi-family mode)
- Predictions require specifying which family via `distribution_column` in `newdata`

## Files and entry points
- `R/fit.R` (family parsing, validation, TMB data/params)
- `R/families.R`, `R/enum.R` (family/link codes and validation helpers)
- `src/sdmTMB.cpp` (likelihood dispatch)
- `R/predict.R`, `R/tmb-sim.R`, `R/residuals.R`, `R/dharma.R`, `R/methods.R`, `R/print.R` (family assumptions, residuals, reporting)
- Reference tinyVAST implementation for the offsets approach: `~/src/tinyVAST/R/fit.R`, `~/src/tinyVAST/src/tinyVAST.cpp`, and vignette `~/src/tinyVAST/vignettes/multiple_data.Rmd`.

## Data model / API changes
- Add optional `distribution_column` argument (default `NULL`).
- If `family` is a named list, require `distribution_column` in `data` and map rows to `e_i`.
- Add TMB data: `e_i` (per-row family index), `family_e`, `link_e`, `multi_family` flag.
- Use single parameter vectors with per-family start/length offsets (tinyVAST-style) for `ln_phi`, `thetaf`, `ln_student_df`, and `gengamma_Q` instead of TMB mapping.
- Store per-family offset metadata in TMB data (e.g., integer vectors for start/length per parameter group) to avoid ambiguity.

## Action items
[ ] Define a small R helper to validate multi-family inputs (named family list, no *_mix families, `distribution_column` values match names) and compute `e_i`.
[ ] Extend `R/fit.R` to branch: current path for single/delta family; new path builds `family_e/link_e/e_i` and per-family parameter offsets for multi-family mode.
[ ] Update `src/sdmTMB.cpp` to dispatch likelihoods using `family_e(e_i(i))` and `link_e(e_i(i))` when `multi_family == 1`, reusing existing switch logic with minimal duplication.
[ ] Handle delta families in multi-family: when `e_i(i)` maps to a delta family, evaluate both binomial and positive likelihoods (2 terms); when it maps to non-delta family, evaluate single likelihood (1 term).
[ ] Enforce residual behavior: in multi-family mode, block standard residuals and direct users to simulation/DHARMa residuals.
[ ] Update predictions to require `distribution_column` in `newdata` for multi-family and pass `e_g` into TMB.
[ ] Add documentation and a vignette or example showing mixed data types with `distribution_column`, ideally mirroring the tinyVAST multiple-data vignette example.

## Implementation stages

Implement incrementally with testing at each stage to avoid getting overwhelmed.

### Stage 0: Infrastructure (Foundation)
**Goal:** Set up basic data structures without changing likelihood evaluation

**Tasks:**
- [x] Add family/link enum helpers in `R/enum.R` for multi-family
- [x] Add validation function for named family lists (no *_mix, check names)
- [x] Add `distribution_column` parameter to `sdmTMB()` (but error if used initially)
- [x] Add basic tests for validation logic

**Test:** Validation catches bad inputs; existing single-family fits still work

**Why first:** No C++ changes, pure R infrastructure. Easy to test and debug.

**Implementation notes (Stage 0):**
- Added `.enum_family()` / `.enum_link()` helpers in `R/enum.R` to reuse enum lookups.
- Added `.validate_multi_family_list()` in `R/multi-likelihood.R` to enforce:
  - named list with unique, non-empty names
  - all elements inherit from `family`
  - excludes any `_mix` families
  - verifies family/link enums are supported
  - optional `distribution_column` mapping returns `e_i` (for later stages)
- Added `distribution_column` arg to `sdmTMB()` with an explicit error guard.
- Added unit tests in `tests/testthat/test-multi-family-validation.R`.

### Stage 1: Simple multi-family (non-delta, no extra parameters)
**Goal:** Get basic multi-family working with simplest families

**Scope:** Only families that don't need extra parameters: **gaussian, poisson, binomial**

**Tasks:**
- [x] Extend `R/fit.R` to build `e_i`, `family_e`, `link_e` arrays for multi-family mode
- [x] Add `multi_family` flag to TMB data
- [x] Update `src/sdmTMB.cpp` to dispatch using `e_i` for these three families only
- [x] Create simple test: mix poisson + binomial + gaussian data

**Test:** Fit a model with 3 different families, check:
- `e_i` maps correctly
- Likelihoods are evaluated correctly (check nll values make sense)
- Can recover known parameters from simulated data
- Existing single-family fits still work
**Add:** a minimal unit test that asserts `e_i` mapping from `distribution_column` to family index.

**Constraint (Stage 1):** Binomial rows must be numeric 0/1 responses, or proportions with `weights` supplying the binomial size (no cbind successes/failures yet).

**Implemented (Stage 1):**
- Added multi-family parsing in `sdmTMB()` for named `family` lists, requiring `distribution_column`.
- Built `multi_family`, `e_i`, `family_e`, and `link_e` in `tmb_data` and stored the original family list on the fit object.
- Added per-row validation for binomial proportions (with `weights` as size) and log-link non-negative response checks.
- Updated C++ to dispatch link/family per-row using `e_i` for gaussian/poisson/binomial.
- Added tests for `distribution_column` mapping and Stage 1 TMB data creation.

**Why second:** Gets core architecture working with minimal complexity. No parameter offsets yet.

### Stage 2: Add extra parameter families (non-delta)
**Goal:** Implement parameter offset logic for families with phi, theta, df, Q

**Scope:** Add **nbinom1, nbinom2, Gamma, lognormal, tweedie, student, gengamma, Beta, betabinomial**

**Tasks:**
- [x] Implement tinyVAST-style parameter vectors with offsets for `ln_phi`, `thetaf`, `ln_student_df`, `gengamma_Q`
- [x] Update TMB data structure with offset indices
- [x] Extend C++ likelihood dispatch to handle these families
- [x] Create test mixing families with different parameter needs (e.g., poisson + nbinom2 + gamma)

**Test:**
- Fit mixed model with nbinom2 + gamma (both need phi)
- Check parameters are indexed correctly
- Verify parameter estimates are reasonable
- Test that families sharing phi get separate estimates

**Notes (Stage 2):**
- Fixed student df is not supported for multi-family yet (error if supplied).
- Gamma/lognormal rows must have response values > 0 in multi-family mode.
- Stage 2 test uses `newton_loops = 0` and `getsd = FALSE` to avoid tiny-dataset Hessian issues.

**Why third:** Builds on Stage 1 infrastructure, adds complexity of parameter management.

### Stage 3: Add delta families
**Goal:** Support delta families in multi-family mode

**Scope:** **delta_gamma, delta_lognormal, delta_truncated_nbinom2**, etc.

**Tasks:**
- [x] Update C++ to check if `e_i(i)` is delta family and evaluate 2 likelihoods
- [x] Handle formula1/formula2 for delta parts
- [x] Ensure parameter indexing works for both delta parts
- [x] Create test mixing delta_gamma + poisson + binomial

**Test:**
- Fit model with delta_gamma + non-delta families
- Check observations with delta family contribute 2 likelihood terms
- Verify binomial and positive parts both fit correctly
- Test with multiple delta families in same model

**Notes (Stage 3):**
- Delta families now allowed in multi-family; non-delta rows only contribute to component 1.
- Per-family `poisson_link_delta` is supported via `poisson_link_delta_e`.
- Student-t `student_df` SEs are derived from `ln_student_df` when `student_df` SE is unnamed in `sdreport`.
- Avoid `$` when accessing `tmb_map$ln_student_df` because partial matching can hit `ln_student_df_e`; use `[[` instead.

**Why fourth:** Most complex part; better to have Stages 1-2 solid first. Delta logic is self-contained.

### Stage 3.5: tinyVAST parity check
**Goal:** Compare multi-family fits against the tinyVAST multiple-data vignette before adding prediction logic. See the file: ~/src/tinyVAST/vignettes/multiple_data.Rmd

**Tasks:**
- Recreate the tinyVAST vignette example data and model spec in sdmTMB.
- Fit both tinyVAST and sdmTMB with the same data and compare key parameter estimates and likelihoods.
- Document any differences and adjust multi-family implementation if needed.

**Why here:** Validates likelihood behavior before predictions/residuals add extra complexity.

**Completed (Stage 3.5):**
- Script: `scratch/multi-family-stage3.5.R` (tinyVAST installed + sdmTMB load_all; random fields off).
- Matched within 1e-3 for fixed effects, Tweedie phi/p, and logLik; max abs diffs: fixed effects ~1.3e-05, phi ~7.8e-04, tweedie_p ~1.5e-06, logLik ~7.4e-07.
- Added prediction parity checks against tinyVAST (`predict(..., what = "p_g")` and `what = "mu_g"`). Link-scale `p_g` matches within 1e-3 absolute; response-scale `mu_g` uses relative tolerance (max rel diff ~2.1e-05). Script now reports `rel_diff` and uses absolute vs relative tolerances appropriately.
- Added index parity checks (tinyVAST `integrate_output()` vs sdmTMB `get_index()`), looping over years with the vignette-style grid if `sf` is available. Uses relative tolerance for index totals (max rel diff ~2.1e-05). Script now orders `Family` with `Biomass_KG` first to ensure index link consistency for the multi-family path.
- Added a single-family delta check for `delta_gamma(type = "poisson-link")` (no multi-family) using only `Biomass_KG` data:
  - tinyVAST delta formula uses `delta_options = list(formula = ~ .)` so the second predictor matches the base formula (including offset).
  - sdmTMB link predictions for delta don't return combined `est`, so the link-scale `p_g` comparison uses `est1 + est2` (poisson-link definition).
  - tinyVAST gamma dispersion is parameterized as CV (`log_sigma`), so phi comparisons convert to shape via `exp(-2 * log_sigma)` to match sdmTMB's `ln_phi`.
  - Fixed effects (both components), phi, logLik, and `p_g`/`mu_g` predictions now match within tolerance.
- Added a gamma-specific phi conversion in `scratch/multi-family-stage3.5.R` so `extract_tv_params()` reports Gamma phi in the same parameterization as sdmTMB.

**Open issue:** the mixed multi-family case with `delta_gamma(type = "poisson-link")` + `binomial(cloglog)` still mismatches on the Encounter component (fixed effects + `p_g`). This looks consistent with offset handling: sdmTMB applies the global `offset` only to the delta positive component when any delta family is present, so non-delta families (Encounter) do not get the offset in component 1. tinyVAST includes the offset in the binomial predictor via the formula.

**Resolution:**
- Updated `src/sdmTMB.cpp` to apply offsets per-row in multi-family predictions:
  - non-delta rows now add the offset to component 1
  - delta rows add the offset to component 2 unless poisson-link delta, which keeps the offset in component 1 (matching tinyVAST’s `p1_g` definition)
- Updated `scratch/multi-family-stage3.5.R` to compare `pred_p_g_pl` using tinyVAST `p_g` for delta rows and `p1_g` for Encounter rows (per the paper).
- Re-ran `devtools::load_all()` + `scratch/multi-family-stage3.5.R`: all comparisons within tolerance, including the multi-family delta + cloglog case.

### Stage 4: Predictions
**Goal:** Make predictions work with multi-family models

**Tasks:**
- Update `predict.sdmTMB()` to require `distribution_column` in newdata
- Error if `newdata` has family names not present in the fitted multi-family model
- Add `e_g` (family index for prediction grid) to TMB
- Handle prediction type (link, response, etc.) per family
- Update prediction tests

**Test:**
- Predict from multi-family model with newdata containing distribution_column
- Test prediction fails gracefully without distribution_column
- Test predictions for delta vs non-delta families
- Verify predictions match single-family model when all rows use same family

**Why fifth:** Predictions depend on fitting working correctly. Can't test until models fit.

**Completed (Stage 4):**
- Multi-family prediction support implemented in `R/predict.R` with minimal branching:
  - Requires `distribution_column` in `newdata`; defaults to original data when `newdata = NULL`.
  - Validates family names and maps `e_g` via `.multi_family_predict_e_g()` in `R/multi-likelihood.R`.
  - Per-row link/response transforms handled by `.multi_family_predict_est()`; mixed delta/non-delta returns `est` plus `est1`/`est2` with `NA` for non-delta rows (choice A).
  - `e_g` added to `tmb_data` for prediction runs; `distribution_column` stored on fit object.
  - `src/sdmTMB.cpp` now declares `DATA_IVECTOR(e_g)` for prediction compatibility.
- Multi-family prediction tests added in `tests/testthat/test-multi-family-predict.R`.
- Simulation predictions:
  - `type = "response"` works for non-delta multi-family via per-row link inversion.
  - Multi-family delta sims with `type = "response"` or `model = NA` are blocked with clear errors (pending future support).

**Implementation notes (Stage 4):**
- New helpers in `R/multi-likelihood.R` keep `predict.sdmTMB()` maintainable:
  - `.multi_family_predict_e_g()` reuses `.validate_multi_family_list()` for mapping and validation.
  - `.multi_family_predict_est()` centralizes per-row transform logic and poisson-link delta handling.
- `predict.sdmTMB()` now sets `multi_family`, `has_delta_multi`, and `dist_col` once near the top and uses these to keep branching localized.
- Tests run: `tests/testthat/test-multi-family-predict.R` and `tests/testthat/test-multi-family-stage3.R` (both pass).


### Stage 5: Residuals and diagnostics
**Goal:** Block standard residuals, ensure simulation residuals work

**Tasks:**
- Block `residuals.sdmTMB()` for multi-family (error with helpful message)
- Ensure `simulate.sdmTMB()` and DHARMa residuals work using per-row family simulation
- Update `print.sdmTMB()` to show multi-family info
- Update `summary.sdmTMB()` for multi-family

**Test:**
- Verify standard residuals give clear error
- Test DHARMa residuals work for multi-family
- Check print/summary show correct family information

**Completed (Stage 5):**
- `residuals.sdmTMB()` now errors for multi-family fits with guidance to use `simulate(..., type = "mle-mvn")` + `dharma_residuals()`.
- `simulate.sdmTMB()` combines multi-family delta simulations into a single response per row (non-delta: component 1; delta: component1 * component2), so DHARMa gets the expected `n_obs x n_sim` matrix.
- `dharma_residuals()` now treats multi-family delta responses like delta models when choosing observed values.
- `fitted.sdmTMB()` now uses `predict(..., type = "response")` for multi-family fits (needed for DHARMa).
- `print.sdmTMB()` now shows multi-family family info in the header and exits early with a note that detailed fixed-effect summaries will be added later (Stage 6).
- Tests added: `tests/testthat/test-multi-family-residuals.R` (residuals error + multi-family delta simulate shape).
- Tests run: `devtools::load_all(); testthat::test_file("tests/testthat/test-multi-family-residuals.R")`.

### Stage 6: Print/tidy support
**Goal:** Handle multi-family family lists in `print()`/`summary()`/`tidy()`

**Tasks:**
- Update `print.sdmTMB()` and `tidy.sdmTMB()` to tolerate `family` as a named list
- Ensure summary output lists per-family names and links without assuming a single family

**Completed (Stage 6):**
- `tidy.sdmTMB()` now supports multi-family fits:
  - Handles delta vs non-delta via `tmb_data$delta_family_e`
- Skips single-family-only checks and adds per-family `phi_*/tweedie_p_*/student_df_*/gengamma_Q_*` rows (filtered by component for delta fits)
  - Guards `rho_time_unscaled` when absent
- Per-family parameter rows use the family list names (e.g., `phi_gaussian`, `tweedie_p_Biomass_KG`, `gengamma_Q_bio`); use `model = 1` (non-delta families) or `model = 2` (delta families) when both are present.
- `print.sdmTMB()` now prints multi-family fixed effects and random parameters using the updated `tidy()` output.
- `print_other_parameters()` now formats per-family dispersion parameters for multi-family fits.
- `print_anisotropy()` detects multi-family delta fits.
- Tests added: `tests/testthat/test-multi-family-tidy.R` (fixed + ran_pars for multi-family).
- Tests run: `devtools::load_all(); testthat::test_file("tests/testthat/test-multi-family-tidy.R")`.

### Stage 6.5: Add beta and betabinomial

**Completed (Stage 6.5):**
- Allowed `Beta` and `betabinomial` families in multi-family models and include them in per-family `phi` offsets.
- Added shared binomial-like validation/transform logic for binomial (proportions only) and betabinomial (counts/proportions) in `R/multi-likelihood.R`.
- Added tests for binomial proportions with `weights` and betabinomial size handling in `tests/testthat/test-multi-family-stage1.R` and `tests/testthat/test-multi-family-stage2.R`.

### Stage 7: Documentation and examples
**Goal:** Make it usable and documented

**Tasks:**
- Add vignette with example (mirror tinyVAST vignette)
- Document `distribution_column` parameter
- Add examples to help files
- Create comprehensive integration tests

**Progress (Stage 7):**
- Drafted `vignettes/articles/multi-family-dogfish.Rmd` using the `dogfish` dataset with binomial (cloglog), nbinom2, and delta_lognormal (poisson-link) families and `log(area_swept)` offset.

### Handoff notes for Stage 2+
- Binomial rows allow 0/1 or proportions with `weights` supplying size; cbind responses and binomial counts are not supported in multi-family mode.
- `print()`/`tidy()` currently break for multi-family objects; update these once core fitting is stable.
- The fit object stores the original named family list; C++ uses `family_e`/`link_e` for multi-family and `family`/`link` for single-family.
- `phi` reporting now checks `family_e` when `multi_family == 1`; preserve this behavior when adding parameter offsets.
- Multi-component and delta families are explicitly blocked in Stage 1; revisit those checks in Stage 2/3 as support expands.

### Key principles
1. **Each stage is independently testable** - can commit and verify correctness
2. **Backward compatibility preserved** - existing code keeps working
3. **Complexity added gradually** - simple families → parameters → delta → predictions
4. **Early stages unlock later stages** - dependencies are clear
5. **Can stop at any stage** - even partial implementation is useful

## Testing and validation
- Add a minimal test for mixed families (e.g., binomial + poisson + tweedie) checking `e_i` mapping and model fit.
- Add a delta-family mixed test (e.g., delta_gamma + binomial) if supported by the chosen delta design.
- Added `tests/testthat/test-multi-family-recovery.R` to simulate a shared linear predictor across several families (including a delta family) and verify fixed-effect recovery.
- Add a prediction test that fails when `distribution_column` is missing in `newdata` for multi-family fits.
- Create a test/example that matches the tinyVAST vignette dataset and model setup as closely as possible.

## Risks and edge cases
- Per-family parameter vectors must align with families that do or don’t use `phi`/`thetaf`/`df`/`gengamma_Q`.
- Performance: per-row branching in likelihood loops; keep branching tight and reuse existing code paths.
