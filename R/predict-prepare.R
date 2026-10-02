# Build the prediction `tmb_data` for `predict.sdmTMB()` with `newdata`.
#
# Returns `tmb_data` with the `proj_*` slots filled and the internal frame
# `nd` (the user's rows plus any dummy columns needed for model matrices,
# including a fake response for `model.frame()` and smoothers).
predict_prepare <- function(object, req, tmb_data, nonlocal_newdata) {
  newdata <- req$newdata
  has_nonlocal <- !is.null(object$nonlocal_parsed)
  nonlocal_uses_external_grid <- .nonlocal_uses_external_grid(object, nonlocal_newdata)

  newdata <- predict_check_newdata(object, req, newdata, nonlocal_newdata)
  spatial <- predict_spatial_index(object, req, newdata)
  newdata <- spatial$newdata

  if (length(object$formula) == 1L) {
    # this formula has breakpt() etc. in it:
    thresh <- list(check_and_parse_thresh_params(object$formula[[1]], newdata))
    formula <- list(thresh[[1]]$formula) # this one does not
  } else {
    thresh <- list(check_and_parse_thresh_params(object$formula[[1]], newdata),
      check_and_parse_thresh_params(object$formula[[2]], newdata))
    formula <- list(thresh[[1]]$formula, thresh[[2]]$formula)
  }
  threshold_columns <- unique(unlist(
    lapply(thresh, `[[`, "threshold_parameter"),
    use.names = FALSE
  ))

  response <- get_response(object$formula[[1]])
  if (!response %in% names(newdata)) newdata[[response]] <- 0 # fake for model.matrix

  .check_no_missing_covariates(
    data = newdata,
    formulas = lapply(formula, reformulas::nobars),
    shared_formulas = c(
      .formula_list(object$spatial_varying_formula),
      .formula_list(object$time_varying),
      list(object$dispformula)
    ),
    # A supplied (or stored) nonlocal grid provides its own covariates. The
    # prediction locations therefore need not duplicate those columns.
    # `.prepare_nonlocal_grid_inputs()` validates an overriding grid below.
    required_columns = if (is.null(object$nonlocal_formula_parsed) || nonlocal_uses_external_grid) {
      threshold_columns
    } else {
      c(object$nonlocal_formula_parsed$covariates, threshold_columns)
    },
    stage = "prediction"
  )

  design <- .newdata_design(object, newdata,
    include_iid = sum(object$tmb_data$n_re_groups) > 0 && isFALSE(req$pop_pred_iid),
    warn_new_levels = isFALSE(req$allow_new_levels))
  Zt_list <- design$Zt_list
  proj_X_ij <- predict_append_nonlocal(object, req, design$X_ij)
  proj_Xdisp_ij <- predict_disp_matrix(object, newdata)
  proj_X_rw_ik <- predict_tv_matrix(object, newdata)

  n <- nrow(proj_X_ij[[1]])
  tmb_data$proj_offset_i <- predict_offset(req, object, n)
  tmb_data$proj_X_threshold <- thresh[[1]]$X_threshold # TODO DELTA HARDCODED TO 1
  tmb_data$area_i <- rep(1, n)
  tmb_data$proj_mesh <- spatial$proj_mesh
  tmb_data$proj_X_ij <- proj_X_ij
  tmb_data$proj_Xdisp_ij <- proj_Xdisp_ij
  tmb_data$proj_X_rw_ik <- proj_X_rw_ik
  tmb_data$Zt_list_proj <- Zt_list
  time_lu <- object$time_lu
  tmb_data$proj_year <- time_lu$year_i[match(newdata[[object$time]], time_lu$time_from_data)]
  tmb_data$proj_time_include <- as.integer(time_lu$time_from_data %in% newdata[[object$time]])
  tmb_data$proj_family_id <- .family_spec_row_family_id(req$family_spec, newdata) - 1L
  if (req$is_areal) {
    tmb_data$proj_lon <- rep(0, nrow(newdata))
    tmb_data$proj_lat <- rep(0, nrow(newdata))
  } else {
    tmb_data$proj_lon <- if (req$xy_cols[[1]] %in% names(newdata)) newdata[[req$xy_cols[[1]]]] else rep(0, nrow(newdata))
    tmb_data$proj_lat <- if (req$xy_cols[[2]] %in% names(newdata)) newdata[[req$xy_cols[[2]]]] else rep(0, nrow(newdata))
  }
  tmb_data$calc_se <- as.integer(req$se_fit)
  tmb_data$pop_pred <- as.integer(req$pop_pred)
  tmb_data$exclude_RE <- req$exclude_RE
  tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
  tmb_data$covariate_diffusion$proj_covariate_vertex_time <- if (has_nonlocal) {
    predict_nonlocal_field(object, tmb_data, newdata, spatial$proj_mesh,
      nonlocal_newdata, nonlocal_uses_external_grid)
  } else {
    array(0, dim = c(1L, 1L, 1L))
  }
  tmb_data$proj_Zs <- design$Zs
  tmb_data$proj_Xs <- design$Xs
  tmb_data$proj_z_i <- predict_svc_matrix(object, newdata)

  list(tmb_data = tmb_data, nd = newdata)
}

# Validate `newdata` coordinates and time; fill in a time column where the
# model has none and dummy coordinates for population predictions.
predict_check_newdata <- function(object, req, newdata, nonlocal_newdata) {
  no_spatial <- as.logical(object$tmb_data$no_spatial)
  has_nonlocal <- !is.null(object$nonlocal_parsed)

  needs_xy <- if (has_nonlocal) TRUE else isFALSE(req$pop_pred) && !no_spatial && !req$is_areal
  if (any(!req$xy_cols %in% names(newdata)) && needs_xy)
    cli_abort(c("`xy_cols` (the column names for the x and y coordinates) are not in `newdata`.",
      "Did you miss specifying the argument `xy_cols` to match your data?",
      "The newer `make_mesh()` (vs. `make_spde()`) takes care of this for you."))

  if (isFALSE(req$pop_pred) && !req$is_areal && (!no_spatial || has_nonlocal)) {
    xy_orig <- object$data[,req$xy_cols]
    xy_nd <- newdata[,req$xy_cols]
    all_outside <- function(x1, x2) {
      min(x1) > max(x2) || max(x1) < min(x2)
    }
    if (all_outside(xy_orig[,1], xy_nd[,1]) || all_outside(xy_orig[,2], xy_nd[,2])) {
      cli_warn(c("`newdata` prediction coordinates appear to be outside the fitted coordinates.",
        "This will likely cause all your random field values to be returned as 0.",
        "Check your coordinates including any conversions between projections.",
        "If working with UTMs, are both in km or m?"))
    }
  }

  if (object$time == "_sdmTMB_time") newdata[[object$time]] <- 0L

  check_time_class(object, newdata)
  original_time <- object$time_lu$time_from_data
  new_data_time <- unique(newdata[[object$time]])

  if (!all(new_data_time %in% original_time))
    cli_abort(c("Some new time values were found in `newdata`. ",
      "If you would like to predict on new time values,",
      "see the `extra_time` argument in `?sdmTMB`.")
    )
  if (.nonlocal_prediction_requires_full_time(object, nonlocal_newdata) &&
    !setequal(new_data_time, original_time)) {
    cli_abort(c(
      "Temporal nonlocal prediction currently requires full time coverage in `newdata`.",
      "i" = "Include exactly the same time values used in the fitted model."
    ))
  }

  # If making population predictions (with standard errors), we don't need
  # to worry about space, so fill in dummy values if the user hasn't made any:
  if (req$pop_pred && !req$is_areal) {
    for (i in c(1, 2)) {
      if (!req$xy_cols[[i]] %in% names(newdata)) {
        suppressWarnings({
          newdata[[req$xy_cols[[i]]]] <- mean(object$data[[req$xy_cols[[i]]]], na.rm = TRUE)
        })
      }
    }
  }

  if (sum(is.na(new_data_time)) > 0)
    cli_abort(c("There is at least one NA value in the time column.",
      "Please remove it."))

  newdata
}

# Add `sdm_orig_id` and `sdm_spatial_id` to `newdata` and build the
# projection matrix from spatial vertices (or areal units) to unique
# prediction locations.
predict_spatial_index <- function(object, req, newdata) {
  no_spatial <- as.logical(object$tmb_data$no_spatial)
  has_nonlocal <- !is.null(object$nonlocal_parsed)
  newdata$sdm_orig_id <- seq(1L, nrow(newdata))

  if (req$is_areal) {
    newdata[["sdm_spatial_id"]] <- seq_len(nrow(newdata)) - 1L
    if (isFALSE(req$pop_pred) && !no_spatial) {
      if (!object$spde$space_column %in% names(newdata)) {
        cli_abort("Areal space column {.field {object$spde$space_column}} was not found in `newdata`.")
      }
      if (anyNA(newdata[[object$spde$space_column]])) {
        cli_abort("Areal space column {.field {object$spde$space_column}} contains missing values in `newdata`.")
      }
      proj_mesh <- areal_projection_matrix(object$spde, newdata)
    } else {
      proj_mesh <- Matrix::sparseMatrix(
        i = integer(0L),
        j = integer(0L),
        x = numeric(0L),
        dims = c(nrow(newdata), object$spde$n_s)
      )
    }
  } else if (!no_spatial || has_nonlocal) {
    locations <- .project_unique_locations(object$spde$mesh, newdata, req$xy_cols)
    newdata[["sdm_spatial_id"]] <- locations$index
    proj_mesh <- locations$A
  } else {
    proj_mesh <- object$spde$A_st # fake
    if (!all(object$spde$xy_cols %in% names(newdata))) {
      newdata[[req$xy_cols[1]]] <- NA_real_ # fake
      newdata[[req$xy_cols[2]]] <- NA_real_ # fake
    }
    newdata[["sdm_spatial_id"]] <- rep(0L, nrow(newdata)) # fake
  }
  list(newdata = newdata, proj_mesh = proj_mesh)
}

# Append the nonlocal coefficient columns to each linear predictor's fixed-effect
# model matrix.
predict_append_nonlocal <- function(object, req, proj_X_ij) {
  if (!is.null(object$nonlocal_parsed)) {
    proj_X_ij[[1]] <- .append_nonlocal_coef_columns(
      X = proj_X_ij[[1]],
      coef_names = object$nonlocal_parsed$term_coef_name
    )
    if (req$has_two_components) {
      proj_X_ij[[2]] <- .append_nonlocal_coef_columns(
        X = proj_X_ij[[2]],
        coef_names = object$nonlocal_parsed$term_coef_name
      )
    }
  }
  proj_X_ij
}

predict_disp_matrix <- function(object, newdata) {
  if (!isTRUE(object$has_dispformula)) return(NULL)
  predict_design_matrix(object, newdata, "dispformula")
}

predict_tv_matrix <- function(object, nd) {
  if (is.null(object$time_varying)) {
    return(matrix(0, ncol = 1, nrow = 1)) # dummy
  }
  predict_design_matrix(object, nd, "time_varying")
}

# Apply a saved auxiliary design (`spatial_varying`, `time_varying`, or
# `dispformula`) to the prediction rows. New fits always save the design;
# only fits from earlier versions reconstruct it.
predict_design_matrix <- function(object, newdata, which) {
  if (is.null(object$design_specs)) {
    return(predict_design_matrix_legacy(object, newdata, which))
  }
  spec <- object$design_specs[[which]]
  if (is.null(spec)) {
    cli_abort("Internal error: no saved `{which}` design specification.")
  }
  .apply_design(spec, newdata)
}

# Fallback for fits saved before `design_specs` existed. This rebuilds the
# design from the stored fitted data under the current session's options, so
# the original encoding (e.g., contrasts set by `options()` at fit time)
# cannot always be recovered exactly.
predict_design_matrix_legacy <- function(object, newdata, which) {
  formula <- switch(which,
    spatial_varying = object$spatial_varying_formula,
    time_varying = object$time_varying,
    dispformula = object$dispformula
  )
  mf_fit <- stats::model.frame(stats::terms(formula), object$data,
    na.action = stats::na.pass)
  tt <- attr(mf_fit, "terms") # with `predvars` from the fitted data
  X_fit <- stats::model.matrix(tt, mf_fit)
  mf_new <- stats::model.frame(tt, newdata,
    xlev = stats::.getXlevels(tt, mf_fit), na.action = stats::na.pass)
  if (anyNA(mf_new)) {
    cli_abort("NAs are not allowed in variables used by `{which}` in `newdata`.")
  }
  X <- stats::model.matrix(tt, mf_new, contrasts.arg = attr(X_fit, "contrasts"))
  fit_cols <- switch(which,
    spatial_varying = object$spatial_varying,
    time_varying = colnames(X_fit),
    dispformula = colnames(object$tmb_data$Xdisp_ij)
  )
  if (!all(fit_cols %in% colnames(X))) {
    cli_abort(c(
      "The `{which}` prediction matrix has different columns than the fitted model.",
      "i" = "Check factor levels, contrasts, and transformed covariates in `newdata`."
    ))
  }
  X[, fit_cols, drop = FALSE]
}

# Covariate field (vertex x time) for nonlocal (covariate diffusion) terms.
predict_nonlocal_field <- function(object, tmb_data, nd, proj_mesh,
                                   nonlocal_newdata, nonlocal_uses_external_grid) {
  if (!is.null(nonlocal_newdata)) {
    # override grid: rebuild the field from the supplied nonlocal_newdata
    override_grid_inputs <- .prepare_nonlocal_grid_inputs(
      grid = nonlocal_newdata,
      nonlocal_formula = object$nonlocal_formula_parsed,
      mesh = object$spde$mesh,
      xy_cols = object$spde$xy_cols,
      time = object$time,
      time_df = object$time_lu,
      full_time_vec = object$time_lu$time_from_data,
      time_indexed = .nonlocal_time_indexed_from_object(object)
    )
    .build_nonlocal_tmb_data(
      nonlocal_formula = object$nonlocal_formula_parsed,
      data = override_grid_inputs$data,
      A_st = override_grid_inputs$A_st,
      A_spatial_index = override_grid_inputs$A_spatial_index,
      year_i = override_grid_inputs$year_i,
      n_t = tmb_data$n_t,
      time_values = object$time_lu$time_from_data
    )$covariate_vertex_time
  } else if (nonlocal_uses_external_grid) {
    # reuse the fitted field: same mesh vertices, all time slices already present
    object$nonlocal_parsed$covariate_vertex_time
  } else {
    # no grid was used at fit: rebuild the field from newdata
    .build_nonlocal_tmb_data(
      nonlocal_formula = object$nonlocal_formula_parsed,
      data = nd,
      A_st = proj_mesh,
      A_spatial_index = nd$sdm_spatial_id,
      year_i = tmb_data$proj_year,
      n_t = tmb_data$n_t,
      time_values = object$time_lu$time_from_data
    )$covariate_vertex_time
  }
}

# Spatially varying coefficient covariate matrix.
predict_svc_matrix <- function(object, newdata) {
  if (is.null(object$spatial_varying)) {
    return(matrix(0, nrow(newdata), 0L))
  }
  predict_design_matrix(object, newdata, "spatial_varying")
}
