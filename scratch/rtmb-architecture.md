# How the RTMB backend fits together

The RTMB backend reuses everything `sdmTMB()` already builds (the TMB `data`,
`parameters`, `map`, and `random` lists). Only the objective function differs:
instead of `src/sdmTMB.cpp`, an R function computes the joint negative log
likelihood, and RTMB tapes it.

## Choosing a backend

- `sdmTMBcontrol(backend = getOption("sdmTMB.backend", "tmb"))` (`R/utils.R`).
- `sdmTMB()` reads `control$backend` and stores it on the fit as `fit$backend`.
- `backend_sdmTMB(fit)` (`R/backend.R`) returns that value, or `"tmb"` for old
  fits that don't have it. `update()` keeps the fitted backend unless you pass
  a new `control`.

## The one dispatch point: `R/backend.R`

```
make_sdmTMB_adfun(data, parameters, map, random, backend)
  tmb:  TMB::MakeADFun(..., DLL = "sdmTMB")
  rtmb: prepared <- rtmb_prepare(data)                  # rtmb-prepare.R
        rtmb_validate(data, prepared, parameters, random)
        RTMB::MakeADFun(rtmb_make_objective(prepared), parameters, map, random)
sdreport_sdmTMB(obj)   # TMB::sdreport or RTMB::sdreport
```

Every caller that builds an objective goes through this function: `fit.R`
(the fit, plus the multiphase first stage), `predict.R`, `index.R`,
`tmb-sim.R`, `cross-val.R`, `caic.R`, `project.R`, and `utils.R`. Post-fit
methods rebuild the objective with new `data` (such as the `proj_*` slots)
and read reports by name, so they don't need to know which backend ran.

## Tracing one objective evaluation

```
rtmb_objective(par, prepared)                          rtmb-objective.R
├─ rtmb_value()                strip simulation refs    rtmb-effects.R
├─ rtmb_transform()        -> theta (natural scale)     rtmb-objective.R
│    rtmb_log_field_sd(), rtmb_aniso_H()                rtmb-precision.R
├─ rtmb_latent_effects()   -> effects (+ nll)           rtmb-effects.R
│    rtmb_precision()          Q for SPDE/aniso/barrier/SAR/CAR   rtmb-precision.R
│    rtmb_gmrf()               omega_s, zeta_s
│    rtmb_spatiotemporal_field()  epsilon_st (IID/AR1/RW)
│    rtmb_iid_effects(), rtmb_time_varying(), rtmb_smooth_effects()
├─ rtmb_linear_predictors(prepared$fit) -> fitted       rtmb-predictors.R
│    rtmb_linear_predictor() per component; rtmb_threshold();
│    rtmb_diffusion_terms()                             rtmb-diffusion.R
├─ rtmb_observations()     -> obs (jnll_obs, devresid, y_i)   rtmb-observations.R
│    rtmb_obs_state()          per-family mu, phi, etc.
│    rtmb_obs_log_density()    likelihoods (rtmb_dcenspois(), rtmb_dordbeta(), ...)
│    rtmb_obs_simulate()       when simulating the response
│    rtmb_obs_deviance()       deviance residuals
├─ rtmb_prior_nll()                                     rtmb-priors.R
├─ if projecting (prepared$proj):
│    rtmb_linear_predictors(prepared$proj) -> projected rtmb-predictors.R
│    rtmb_mixture_eta(), rtmb_combined_projection()
│    rtmb_derived_indices()   totals, COG-style weighted averages, EAO
└─ rtmb_report()           REPORT/ADREPORT under C++ names   rtmb-report.R
```

The returned `jnll` = latent effects nll + observation nll + prior nll
(+ the `eps_index` term used for bias correction).

## The objects passed between steps

| Object     | Built by                   | Contents |
| ---------- | -------------------------- | -------- |
| `prepared` | `rtmb_prepare()`           | Named flags, one-based indices, `families`, `precision` templates, `fit`/`proj` rows. Doesn't depend on parameters. |
| `par`      | RTMB                       | Parameters with the C++ names and shapes. |
| `theta`    | `rtmb_transform()`         | Natural-scale parameters (`sigma_O`, `kappa`, `range`, `phi`, `H`, ...). |
| `effects`  | `rtmb_latent_effects()`    | Latent arrays (after any simulation), `log_sigma_E`, `nll`. |
| `rows`     | `rtmb_row_inputs()`        | `prepared$fit` or `prepared$proj`: design matrices, `A` matrices, time/station/family indices. |
| predictors | `rtmb_linear_predictors()` | `fixed`, `smooth`, `rw`, `iid`, `omega`, `epsilon`, `svc`, `fe`, `rf`, `eta` (n by component); `zeta` (n by SVC by component). |

`rtmb_prepare()` and its `rtmb_*_inputs()` helpers are the **only** code that
reads the C++ data list. Everything downstream reads `prepared`.

## Where to look when changing something

| Change | RTMB | Also touch (TMB side) |
| ------ | ---- | --------------------- |
| New data flag / option from `fit.R` | `rtmb_prepare()` (add a named flag), then use it downstream | `fit.R` data list, `src/sdmTMB.cpp` |
| New or modified family | `rtmb_obs_state()`, `rtmb_obs_log_density()`, `rtmb_obs_simulate()`, `rtmb_obs_deviance()`; `rtmb_component_mean()` if its mean isn't the inverse link | `R/families.R`, `R/enum.R`, cpp |
| New parameter transform | `rtmb_transform()`; report it in `rtmb_report()` | cpp REPORT/ADREPORT |
| Spatial field / precision type | `rtmb_precision_inputs()` + `rtmb_precision()`; scale in `rtmb_gmrf_scale()` | cpp |
| Spatiotemporal structure | `rtmb_spatiotemporal_field()` | cpp |
| New random effect type | `rtmb_latent_effects()` + a helper in `rtmb-effects.R`; add its name to `rtmb_validate()`'s allowed random list | cpp |
| New predictor term | `rtmb_linear_predictor()`; add the covariates to both branches of `rtmb_row_inputs()` | `fit.R` fit + `proj_*` data |
| Index / derived quantity | `rtmb_derived_indices()`, then `rtmb_report()` | `index.R`, cpp |
| Prior | `rtmb_prior_inputs()`, `rtmb_prior_nll()` | cpp |
| Report name or shape | `rtmb_report()` (names must match C++; post-fit code indexes them) | cpp |
| Unsupported feature guard | `rtmb_validate()` | — |

## RTMB-specific points to watch for

- `"[<-" <- RTMB::ADoverload("[<-")` at the top of any function that assigns
  AD values into plain vectors/matrices.
- Any branch on a parameter value is fixed at taping time. Branch only on
  `prepared`, or use AD-safe expressions (`logspace_add`, `ifelse`-style
  arithmetic). That is why the censored Poisson log probability
  (`rtmb_censpois_logprob()`) is an `RTMB::ADjoint` with an analytic
  derivative.
- RTMB differentiates linear combinations of constant sparse matrices but not
  products of AD sparse matrices, so `rtmb_precision_inputs()` pre-expands
  every product.
- Simulation: random parameters arrive as `simref`s; `rtmb_latent_effects()`
  draws those listed in `prepared$simulate_re`, and `RTMB::OBS()` marks the
  response for `rtmb_obs_simulate()`.
- `rtmb_make_objective()` builds the closure in its own frame so the tape
  doesn't capture `data`/`parameters`.

## Tests

`tests/testthat/test-rtmb-*.R`, grouped by topic (families, fields, effects,
delta, diffusion, index, postfit, backend). Most compare objectives, gradients,
and reports against TMB with the helpers in `helper-rtmb.R`. Run the whole
suite with `SDMTMB_TEST_BACKEND=rtmb` to use RTMB as the default backend.
