## sdmTMB benchmark: old commit vs current TMB vs current RTMB

Compares sdmTMB at commit `891d8e7f914759338ef91d078032ee85d2ed778b` (TMB only) with the current version using the TMB backend and the default RTMB backend. Each configuration ran in a fresh R process from a non-debug install.

**Models**

- `fit`, `predict`, `get_index()`, `simulate()`: `pcod` data, mesh with `cutoff = 8`, 
  `density ~ 0 + as.factor(year)`, delta-gamma family (`delta_gamma()`), `time = "year"`, 
  spatial and IID spatiotemporal random fields. `predict()` and `get_index()` use 
  `qcs_grid` replicated across all years; `get_index()` is run with `bias_correct` 
  `FALSE` and `TRUE`; `simulate()` uses `nsim = 400`.
- `project()`: `dogfish` data (as in `?project`), mesh with `cutoff = 10`, Tweedie 
  with an area-swept offset, spatial field and AR1 spatiotemporal field, fitted to 
  2004-2022 and projected 30 further years on `wcvi_grid` with `nsim = 400`.

Times are elapsed seconds (median over reps; `fit` is run once). Speedup is 
old TMB time divided by current RTMB time.

### Time (s)

| step | old TMB | current TMB | current RTMB | RTMB speedup vs old |
|:--|--:|--:|--:|--:|
| `fit` | 15.05 | 16.32 | 12.65 | 1.2x |
| `predict` | 0.85 | 0.63 | 0.06 | 13.5x |
| `get_index(bias_correct = FALSE)` | 5.21 | 0.87 | 0.50 | 10.4x |
| `get_index(bias_correct = TRUE)` | 13.59 | 9.24 | 9.21 | 1.5x |
| `simulate(nsim = 400)` | 3.04 | 3.12 | 1.62 | 1.9x |
| `project(nsim = 400, 30 yr AR1)` | 35.62 | 9.07 | 3.81 | 9.3x |

### Peak memory (MB)

| step | old TMB | current TMB | current RTMB |
|:--|--:|--:|--:|
| `fit` | 888.4 | 888.8 | 1063.4 |
| `predict` | 1022.3 | 1056.7 | 886.7 |
| `get_index(bias_correct = FALSE)` | 1389.1 | 1132.1 | 1397.8 |
| `get_index(bias_correct = TRUE)` | 6557.5 | 6299.0 | 7827.4 |
| `simulate(nsim = 400)` | 6263.9 | 6130.6 | 7415.3 |
| `project(nsim = 400, 30 yr AR1)` | 7117.4 | 4159.1 | 4391.7 |

Peak memory is the peak process RSS during the step (includes baseline session memory).
