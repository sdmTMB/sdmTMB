# Reference fits

A regression suite of cached model results. Every case is refit and compared
with `ref/reference.csv` so that changes to sdmTMB can't silently change old
model results. It runs in CI (`.github/workflows/reference-fits.yaml`) for
both the TMB and RTMB backends against the same reference file.

This folder is excluded from the package build (`.Rbuildignore`).

## Running

From the package root:

```sh
make reference-check     # both backends, optimized DLL, 6 cores (REF_CORES=N to change)
make reference-record    # re-record after an accepted change

Rscript reference-fits/run.R check --backend tmb|rtmb|both [--filter REGEX] [--cores N]
Rscript reference-fits/run.R validate [--backend ...] [--filter REGEX] [--cores N]
Rscript reference-fits/run.R record [--filter REGEX] [--cores N]
```

- `check` refits each case (no standard errors) and compares its metrics with
  the reference. Failures are printed and written to `check-<backend>.csv`
  (git-ignored; uploaded as an artifact in CI).
- `validate` fits with standard errors, runs `sanity()`, reports timing, and
  compares TMB with RTMB. It writes nothing.
- `record` is `validate` for both backends. It writes the TMB values to
  `ref/reference.csv` (and versions to `ref/session.txt`) only if every case
  passes `sanity()` and the two backends agree within tolerance. Cases limited
  to one backend (`backends` in `cases.R`, e.g., RTMB-only families) record
  that backend's values and skip the comparison. With
  `--filter`, only the matching cases are replaced.

## Accepting a change

If a change is supposed to alter results, run `make reference-record`, review
`git diff reference-fits/ref/`, and commit the new reference with an
explanation of why the values moved.

## What's here

- `make-data.R`: one-time generator for the fixtures in `data/`. Responses are
  simulated in base R from the model each case fits, with parameters chosen
  to be clearly identifiable. The CSV files are the fixtures; rerunning the
  script is only for deliberately replacing them (then re-record).
- `data/`: frozen fixtures, including the mesh as vertices and triangles
  (`mesh-loc.csv`, `mesh-tv.csv`) so results don't depend on fmesher's mesh
  construction, and the Ohio county adjacency (`ohio-edges.csv`, used with
  `ohio_df`) so sf isn't needed.
- `cases.R`: the case table.
- `run.R`: fitting, metrics, comparison, and recording.

## Cases

- `family/`: every family and documented link, non-spatial (fast), including
  binomial with `cbind()` trials and with proportions plus `weights`, and the
  censored families (all but `censored_poisson()` are RTMB-only), including
  `censored_method = "direct"` and an observation-level random intercept.
- `offset/`: one family per link with an offset.
- `delta/`: every delta family, non-default `link1`/`link2`,
  `type = "poisson-link"`, and offsets.
- `fields/`: spatial and spatiotemporal configurations (IID, AR1, RW),
  `share_range`, anisotropy, SVCs, time-varying (RW, RW0, AR1),
  `extra_time`, REML, IID random intercepts and correlated slopes, smoothers,
  threshold models, `dispformula`, priors, `bayesian = TRUE`, count and
  positive-continuous families with fields, multi-family
  (`distribution_column`), nonlocal (covariate diffusion), and delta models
  with shared and separate component structure.
- `areal/`: SAR and CAR.

All fits use `multiphase = FALSE`. Identity and inverse links get starting
values near the generating intercept, as a user would supply.

## Metrics

Each case records `nll` (objective at the optimum), `loglik`, `n_fixed`, and
the mean and SD of link-scale predictions on the fitted data and on a fixed
`newdata` (both components for delta models), plus the mean response-scale
prediction. Some cases add population-level predictions (`re_form = NA`),
`get_index()` estimates and standard errors, and `get_cog()`.

`record` also stores `seconds_tmb` and `seconds_rtmb`: each case's elapsed
time (fit with standard errors, metrics, and `sanity()`) to 0.01 s. These are
for information only and are not compared by `check`.

Every `check` and `validate` run also appends to a local, git-ignored
`reference-fits/timings-local.csv` (latest timing per mode/backend/case, to
0.01 s) and flags cases at least 10% (and 0.02 s) slower than the previous run.
Thresholds are `slow_rel` and `slow_abs` in `run.R`.

A metric passes if `|new - ref| <= abs + rel * |ref|`. Defaults are in
`run.R` (`default_tol`), and cases can override them with `tol` in `cases.R`.
Summaries at the optimum are used because they are insensitive to the path
the optimizer took; each case's data contain the structure the case
estimates, so no variance sits on a boundary where tiny numeric changes could
move it.

## Adding a case

Add it to `cases.R`. If it needs new data, add the column at the end of its
block in `make-data.R` under its own `set.seed()`, so the random streams of
existing columns don't shift. Rerun the script and confirm that `check` still
passes for every existing case. Check the new case with
`Rscript reference-fits/run.R validate --filter <name>`, then
`Rscript reference-fits/run.R record --filter <name>`.
