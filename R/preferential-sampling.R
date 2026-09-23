# Preferential-sampling grid-cell presence/absence sub-model.
#
# Mirrors the shape of R/covariate-diffusion.R: a self-contained set of
# internal helpers, threaded into sdmTMB() in R/fit.R at the same points
# nonlocal_formula/nonlocal_data are threaded in. See the project design
# doc "preferential-sampling-integration-design.md" for the full rationale.
#
# First-pass scope (deliberately excluded, matching the RTMB reference this
# was ported from): smooths, time-varying effects, IID random effects/
# slopes, threshold effects, and spatially varying coefficients do not
# contribute to the grid-level linear predictor. Only the plain fixed-effect
# design matrix (X_pref_ij %*% b_j, model column 0 only -- not delta-aware)
# plus the spatial (omega_s) and spatiotemporal (epsilon_st) fields,
# reprojected onto the grid via A_pref, are reused from the catch model.

#' @noRd
.validate_preferential_args <- function(preferential_grid, preferential_response,
                                         preferential_b_type, mesh_missing) {
  if (!is.null(preferential_grid) && is.null(preferential_response)) {
    cli_abort("`preferential_grid` was supplied but `preferential_response` was not.")
  }
  if (!is.null(preferential_response)) {
    if (!is.character(preferential_response) || length(preferential_response) != 1L) {
      cli_abort("`preferential_response` must be a single column name (character string).")
    }
    if (mesh_missing) {
      cli_abort("`mesh` must be supplied when using `preferential_response`.")
    }
  }
  match.arg(preferential_b_type[1L], c("constant", "rw", "iid"))
}

#' @noRd
.default_preferential_grid <- function(preferential_grid, preferential_response, data) {
  if (is.null(preferential_response)) {
    return(NULL)
  }
  if (is.null(preferential_grid)) data else preferential_grid
}

#' Validate and prepare the preferential-sampling grid
#'
#' Checks `grid` has usable coordinates, a usable response column, and full
#' time coverage, then builds the mesh-projection matrix and per-row year
#' index used by the sampling sub-model. Mirrors
#' `.prepare_nonlocal_grid_inputs()`.
#' @noRd
.prepare_preferential_grid_inputs <- function(grid,
                                               xy_cols,
                                               time,
                                               time_df,
                                               full_time_vec,
                                               preferential_response,
                                               mesh) {
  if (!inherits(grid, "data.frame")) {
    cli_abort("`preferential_grid` must be `NULL` or a data frame.")
  }
  if (is.null(xy_cols) || length(xy_cols) != 2L) {
    cli_abort("The preferential-sampling grid requires a mesh built with known `xy_cols` (e.g., from `make_mesh()`).")
  }

  missing_xy <- setdiff(xy_cols, names(grid))
  if (length(missing_xy)) {
    cli_abort(c(
      "`preferential_grid` is missing required coordinate column(s).",
      "x" = "Missing: {.code {paste(missing_xy, collapse = ', ')}}"
    ))
  }
  non_numeric_xy <- xy_cols[!vapply(xy_cols, function(col) is.numeric(grid[[col]]), logical(1L))]
  if (length(non_numeric_xy)) {
    cli_abort(c(
      "`preferential_grid` coordinates must be numeric.",
      "x" = "Non-numeric coordinate column(s): {.code {paste(non_numeric_xy, collapse = ', ')}}"
    ))
  }
  invalid_xy <- xy_cols[!vapply(xy_cols, function(col) all(is.finite(grid[[col]])), logical(1L))]
  if (length(invalid_xy)) {
    cli_abort(c(
      "`preferential_grid` coordinates must be finite and cannot contain `NA` values.",
      "x" = "Invalid coordinate column(s): {.code {paste(invalid_xy, collapse = ', ')}}"
    ))
  }

  if (!preferential_response %in% names(grid)) {
    cli_abort("`preferential_grid` is missing the `preferential_response` column {.code {preferential_response}}.")
  }
  r_raw <- grid[[preferential_response]]
  if (is.logical(r_raw)) r_raw <- as.numeric(r_raw)
  if (!is.numeric(r_raw) || anyNA(r_raw) || !all(r_raw %in% c(0, 1))) {
    cli_abort(c(
      "`preferential_grid${preferential_response}` must be a 0/1 (or logical) indicator with no missing values.",
      "i" = "It represents whether each grid cell/row was sampled (1) or not (0) in that time slice."
    ))
  }

  if (!time %in% names(grid)) {
    if (identical(time, "_sdmTMB_time")) {
      grid[[time]] <- 0L # internal placeholder time column; not user-facing
    } else {
      cli_abort("`preferential_grid` is missing the time column {.code {time}}.")
    }
  }
  missing_slices <- setdiff(full_time_vec, grid[[time]])
  if (length(missing_slices)) {
    cli_abort(c(
      "`preferential_grid` does not cover all fitted (+ `extra_time`) time slices.",
      "x" = "Missing time slice(s): {.code {paste(missing_slices, collapse = ', ')}}",
      "i" = "`preferential_grid` must be pre-expanded by the user across every modeled time slice (one row per grid cell per time slice)."
    ))
  }
  year_i <- time_df$year_i[match(grid[[time]], time_df$time_from_data)]
  if (anyNA(year_i)) {
    cli_abort("`preferential_grid` contains time value(s) not present in the fitted (+ `extra_time`) time slices.")
  }

  A_pref <- fmesher::fm_basis(mesh, loc = as.matrix(grid[, xy_cols, drop = FALSE]))

  list(
    data = grid,
    A_pref = A_pref,
    year_i = year_i,
    R_i = as.numeric(r_raw),
    n_pref = nrow(grid)
  )
}

#' Build the preferential-sampling fixed-effect design matrix
#'
#' Reuses the fitted model's own `Terms`/`xlev`/`contrasts` (the same
#' mechanism `predict.sdmTMB()` uses for `newdata`) so that `X_pref_ij` is
#' guaranteed either to have identical columns to `X_ij[[1]]`, or to fail
#' loudly. A factor level in `preferential_grid` that was never observed in
#' `data` triggers R's standard "factor ... has new levels" error from
#' `model.matrix()`, which is caught and re-raised as a `cli_abort()`; an
#' exact column-name/order comparison is done afterward as defense in depth
#' against any subtler mismatch. No attempt is
#' made to pad, reorder, or otherwise reconcile a mismatch -- fitting stops.
#' @noRd
.build_preferential_X <- function(grid, formula_terms, xlev, contrasts, X_main) {
  Terms_pref <- stats::delete.response(formula_terms)
  required_vars <- all.vars(Terms_pref)
  missing_vars <- setdiff(required_vars, names(grid))
  if (length(missing_vars)) {
    cli_abort(c(
      "`preferential_grid` is missing fixed-effect predictor column(s) required by `formula`.",
      "x" = "Missing: {.code {paste(missing_vars, collapse = ', ')}}"
    ))
  }
  mf_pref <- tryCatch(
    stats::model.frame(Terms_pref, grid, xlev = xlev, na.action = stats::na.pass),
    error = function(e) {
      cli_abort(c(
        "Failed to build the preferential-sampling fixed-effect design matrix from `preferential_grid`.",
        "x" = conditionMessage(e)
      ))
    }
  )
  X_pref <- tryCatch(
    stats::model.matrix(Terms_pref, mf_pref, contrasts.arg = contrasts),
    error = function(e) {
      cli_abort(c(
        "Failed to build the preferential-sampling fixed-effect design matrix from `preferential_grid`.",
        "i" = "This usually means `preferential_grid` has a factor level for a fixed-effect predictor that was not present in `data` when the model was fit.",
        "x" = conditionMessage(e)
      ))
    }
  )
  if (!identical(colnames(X_pref), colnames(X_main))) {
    missing_cols <- setdiff(colnames(X_main), colnames(X_pref))
    extra_cols <- setdiff(colnames(X_pref), colnames(X_main))
    msgs <- c(
      "The fixed-effect design matrix built from `preferential_grid` does not match the fitted model's design matrix (`formula`).",
      "i" = "Every fixed-effect predictor in `formula` must be present in `preferential_grid`, producing the same columns in the same order."
    )
    if (length(missing_cols)) {
      msgs <- c(msgs, "x" = paste0("Column(s) missing from `preferential_grid`: ", paste(missing_cols, collapse = ", ")))
    }
    if (length(extra_cols)) {
      msgs <- c(msgs, "x" = paste0("Extra/unexpected column(s) from `preferential_grid`: ", paste(extra_cols, collapse = ", ")))
    }
    cli_abort(msgs)
  }
  X_pref
}

#' Assemble the TMB data list for the preferential-sampling sub-model
#'
#' Takes the grid inputs already validated/projected by
#' `.prepare_preferential_grid_inputs()` (computed earlier in `sdmTMB()`,
#' right after `time_df` is built -- mirroring where nonlocal's grid inputs
#' are prepared) and combines them with the fixed-effect design matrix built
#' from the fitted model's own `formula` terms (available only later in
#' `sdmTMB()`, once `X_ij`/`mf`/`mt` exist). Returns `NULL` when
#' preferential sampling is off (`grid_inputs` is `NULL`), so callers can
#' gate all downstream wiring on `is.null(preferential_tmb)`.
#' @noRd
.build_preferential_tmb_data <- function(grid_inputs,
                                          preferential_b_type,
                                          formula_terms,
                                          xlev,
                                          contrasts,
                                          X_main) {
  if (is.null(grid_inputs)) {
    return(NULL)
  }
  X_pref <- .build_preferential_X(
    grid = grid_inputs$data,
    formula_terms = formula_terms,
    xlev = xlev,
    contrasts = contrasts,
    X_main = X_main
  )
  b_type_int <- switch(preferential_b_type,
    constant = 0L,
    iid       = 1L,
    rw        = 2L
  )
  list(
    n_pref = as.integer(grid_inputs$n_pref),
    R_i = grid_inputs$R_i,
    X_pref_ij = X_pref,
    A_pref = grid_inputs$A_pref,
    year_i_pref = as.integer(grid_inputs$year_i),
    b_pref_type = b_type_int
  )
}

#' Placeholder preferential-sampling TMB data when the feature is off
#'
#' Always-present, zero-length fields, mirroring how `covariate_diffusion`
#' is always passed to TMB with `n_terms = 0` when `nonlocal_formula` is
#' unused. Keeps the `tmb_data` list shape constant regardless of whether
#' preferential sampling is requested.
#' @noRd
.default_preferential_tmb <- function(n_b_j) {
  list(
    n_pref = 0L,
    R_i = numeric(0),
    X_pref_ij = matrix(0, nrow = 0, ncol = n_b_j),
    A_pref = Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, 0L)),
    year_i_pref = integer(0),
    b_pref_type = 0L
  )
}

#' Initial parameter values for the preferential-sampling sub-model
#'
#' Zero-length for every element when the feature is off, matching how
#' `tmb_data$preferential` is always present but empty. When on: `gamma_0`,
#' `ln_tau_xi`, `ln_kappa_xi` start at 0, matching how the main model's own
#' `ln_tau_O`/`ln_kappa` start at 0 (no domain-aware heuristic is used
#' elsewhere in \pkg{sdmTMB}, so none is added here either); `xi_s` starts
#' at all zeros, length `n_s` (mesh vertex count); `b_pref` is length 1 for
#' `preferential_b_type == "constant"` or length `n_t` otherwise;
#' `log_sigma_b_pref` is length 0 for `"constant"` (no rw/iid penalty
#' applies) or length 1 (starting at `log(0.3)`) otherwise.
#' @noRd
.preferential_init_params <- function(is_on, preferential_b_type, n_s, n_t) {
  if (!is_on) {
    return(list(
      gamma_0 = numeric(0),
      b_pref = numeric(0),
      log_sigma_b_pref = numeric(0),
      ln_tau_xi = numeric(0),
      ln_kappa_xi = numeric(0),
      xi_s = numeric(0)
    ))
  }
  list(
    gamma_0 = 0,
    b_pref = if (identical(preferential_b_type, "constant")) 0 else rep(0, n_t),
    log_sigma_b_pref = if (identical(preferential_b_type, "constant")) numeric(0) else log(0.3),
    ln_tau_xi = 0,
    ln_kappa_xi = 0,
    xi_s = rep(0, n_s)
  )
}

#' Names of preferential-sampling parameters to unmap (freely estimate)
#'
#' `map_all_params()` in `R/fit.R` starts every parameter mapped off
#' (`factor(NA)`); callers explicitly unmap whatever should actually be
#' estimated. Unmapping a zero-length parameter is a no-op, so this can be
#' called unconditionally even when the feature is off.
#' @noRd
.preferential_map_names <- function(is_on) {
  if (!is_on) {
    return(character(0))
  }
  c("gamma_0", "b_pref", "log_sigma_b_pref", "ln_tau_xi", "ln_kappa_xi", "xi_s")
}

#' Names of preferential-sampling parameters that are random effects
#'
#' `xi_s` is always a random effect (a GMRF field) whenever the feature is
#' on. `b_pref` is a random effect only for `"rw"`/`"iid"` -- for
#' `"constant"` it's a plain length-1 fixed effect, same treatment as `b_j`.
#' @noRd
.preferential_random_names <- function(is_on, preferential_b_type) {
  if (!is_on) {
    return(character(0))
  }
  nms <- "xi_s"
  if (!identical(preferential_b_type, "constant")) {
    nms <- c(nms, "b_pref")
  }
  nms
}
