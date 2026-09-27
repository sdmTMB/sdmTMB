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

  Zt_list <- predict_iid_re_matrices(object, req, newdata, length(formula))
  proj_X_ij <- predict_fixed_matrices(object, req, newdata)
  proj_Xdisp_ij <- predict_disp_matrix(object, newdata)
  # TODO DELTA hardcoded to 1:
  sm <- parse_smoothers(object$smoothers$formula_no_bars, data = object$data,
    newdata = newdata, basis_prev = object$smoothers$basis_out)
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
  tmb_data$proj_Zs <- sm$Zs
  tmb_data$proj_Xs <- sm$Xs
  tmb_data$proj_z_i <- predict_svc_matrix(object, newdata)
  tmb_data$epsilon_predictor <- predict_epsilon_covariate(object, tmb_data, newdata)

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
    if (requireNamespace("dplyr", quietly = TRUE)) { # faster
      unique_newdata <- dplyr::distinct(newdata[, req$xy_cols, drop = FALSE])
    } else {
      unique_newdata <- unique(newdata[, req$xy_cols, drop = FALSE])
    }
    unique_newdata[["sdm_spatial_id"]] <- seq(1, nrow(unique_newdata)) - 1L

    if (requireNamespace("dplyr", quietly = TRUE)) { # much faster
      newdata <- dplyr::left_join(newdata, unique_newdata, by = req$xy_cols)
    } else {
      newdata <- base::merge(newdata, unique_newdata, by = req$xy_cols,
        all.x = TRUE, all.y = FALSE)
      newdata <- newdata[order(newdata$sdm_orig_id),, drop = FALSE]
    }
    proj_mesh <- fmesher::fm_basis(object$spde$mesh,
      loc = as.matrix(unique_newdata[, req$xy_cols, drop = FALSE]))
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

# IID random intercept/slope sparse model matrices, one per linear predictor.
# All linear predictors share the same random effect structure.
predict_iid_re_matrices <- function(object, req, newdata, n_formula) {
  if (sum(object$tmb_data$n_re_groups) == 0 || !isFALSE(req$pop_pred_iid)) {
    return(list())
  }
  re_formula_no_response <- stats::formula(
    stats::delete.response(
      stats::terms(remove_s_and_t2(object$smoothers$formula_no_sm))
    )
  )
  # factor level checks:
  RE_names <- barnames(reformulas::findbars(re_formula_no_response))
  missing_RE_names <- setdiff(RE_names, names(newdata))
  if (length(missing_RE_names) > 0) {
    cli_abort(c(
      "Random effect group column(s) missing from `newdata`: {.val {missing_RE_names}}.",
      "i" = "Use `re_form_iid = NA` or `re_form_iid = ~0` to exclude random effects in prediction."
    ))
  }
  new_level_rows <- integer(0)
  for (i in seq_along(RE_names)) {
    assert_that(is.factor(newdata[[RE_names[i]]]),
      msg = sprintf("Random effect group column `%s` in newdata is not a factor.", RE_names[i]))
    levels_fit <- levels(object$data[[RE_names[i]]])
    values_nd <- as.character(newdata[[RE_names[i]]])
    is_new_level <- !is.na(values_nd) & !values_nd %in% levels_fit
    if (any(is_new_level)) {
      new_level_rows <- union(new_level_rows, which(is_new_level))
      if (isFALSE(req$allow_new_levels)) {
        cli_warn(c(
          "Found new levels in random effect grouping variable {.field {RE_names[i]}}.",
          "i" = "These rows will use population-level IID random effect predictions (`re_form_iid = NA`).",
          "i" = "Set `allow_new_levels = TRUE` to suppress this warning."
        ))
      }
    }
  }

  # now do with a joint data frame to ensure factor levels match
  common_cols <- intersect(colnames(object$data), colnames(newdata))
  nd_aligned <- newdata[, common_cols, drop = FALSE]
  for (col_name in common_cols) {
    if (is.factor(object$data[[col_name]]) && is.factor(nd_aligned[[col_name]])) {
      nd_aligned[[col_name]] <- factor(
        as.character(nd_aligned[[col_name]]),
        levels = levels(object$data[[col_name]])
      )
      if (anyNA(nd_aligned[[col_name]])) {
        nd_aligned[[col_name]][is.na(nd_aligned[[col_name]])] <-
          levels(object$data[[col_name]])[1]
      }
    }
  }
  joint_df <- rbind(object$data[, common_cols, drop = FALSE], nd_aligned)
  xx <- parse_formula(re_formula_no_response, joint_df)
  # drop the original data:
  Zt <- xx$re_cov_terms$Zt[, seq(nrow(object$data) + 1, nrow(object$data) + nrow(newdata)), drop = FALSE]
  if (length(new_level_rows) > 0) {
    Zt[, new_level_rows] <- 0
  }
  rep(list(Zt), n_formula)
}

# Fixed-effect model matrices, one per linear predictor.
predict_fixed_matrices <- function(object, req, newdata) {
  proj_X_ij <- list()
  for (i in seq_along(object$formula)) {
    f2 <- remove_s_and_t2(object$split_formula[[i]]$form_no_bars)
    tt <- stats::terms(f2)
    attr(tt, "predvars") <- attr(object$terms[[i]], "predvars")
    Terms <- stats::delete.response(tt)
    mf <- model.frame(Terms, newdata, xlev = object$xlevels[[i]])
    proj_X_ij[[i]] <- model.matrix(Terms, mf, contrasts.arg = object$contrasts[[i]])
  }
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
  tt_disp <- stats::terms(object$dispformula)
  mf_disp_fit <- model.frame(tt_disp, object$data)
  xlevels_disp <- stats::.getXlevels(attr(mf_disp_fit, "terms"), mf_disp_fit)
  mf_disp <- model.frame(tt_disp, newdata, xlev = xlevels_disp, na.action = stats::na.pass)
  if (sum(is.na(mf_disp)) > 0) {
    cli_abort("NAs are not allowed in variables used by `dispformula` in `newdata`.")
  }
  proj_Xdisp_ij <- model.matrix(tt_disp, mf_disp)
  fit_disp_cols <- colnames(object$tmb_data$Xdisp_ij)
  pred_disp_cols <- colnames(proj_Xdisp_ij)
  missing_cols <- setdiff(fit_disp_cols, pred_disp_cols)
  extra_cols <- setdiff(pred_disp_cols, fit_disp_cols)
  if (length(missing_cols) > 0 || length(extra_cols) > 0) {
    cli_abort(c(
      "Dispersion model matrix in `newdata` does not match the fitted `dispformula` terms.",
      if (length(missing_cols) > 0) paste0("x Missing terms: ", paste(missing_cols, collapse = ", ")),
      if (length(extra_cols) > 0) paste0("x New terms: ", paste(extra_cols, collapse = ", ")),
      "i" = "Check factor levels and columns used in `dispformula`."
    ))
  }
  proj_Xdisp_ij[, fit_disp_cols, drop = FALSE]
}

predict_tv_matrix <- function(object, nd) {
  if (is.null(object$time_varying)) {
    return(matrix(0, ncol = 1, nrow = 1)) # dummy
  }
  tv_terms <- stats::terms(object$time_varying)
  mf_tv_orig <- stats::model.frame(
    tv_terms,
    object$data,
    na.action = stats::na.pass
  )
  tv_xlevels <- stats::.getXlevels(tv_terms, mf_tv_orig)
  X_tv_orig <- stats::model.matrix(tv_terms, mf_tv_orig)
  tv_contrasts <- attr(X_tv_orig, "contrasts")
  mf_tv_new <- stats::model.frame(
    tv_terms,
    nd,
    xlev = tv_xlevels,
    na.action = stats::na.pass
  )
  proj_X_rw_ik <- stats::model.matrix(
    tv_terms,
    mf_tv_new,
    contrasts.arg = tv_contrasts
  )
  if (!identical(colnames(proj_X_rw_ik), colnames(X_tv_orig))) {
    cli::cli_abort(c(
      "The time-varying prediction matrix has different columns than the fitted model.",
      "This may be caused by changed factor levels, contrasts, or transformed covariates in `newdata`."
    ))
  }
  proj_X_rw_ik
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
    # no grid was used at fit: rebuild the field from newdata, as before
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
  # recreate original data SVC formula stuff:
  z_i_orig <- model.matrix(object$spatial_varying_formula, object$data)
  svc_contrasts <- attr(z_i_orig, which = "contrasts")
  ttsv <- stats::terms(object$spatial_varying_formula)
  mfsv <- model.frame(ttsv, object$data)
  mtsv <- attr(mfsv, "terms")
  xlevelssv <- stats::.getXlevels(mtsv, mfsv)
  # apply it to prediction data:
  mfsv_new <- model.frame(ttsv, newdata, xlev = xlevelssv)
  z_i <- model.matrix(ttsv, mfsv_new, contrasts.arg = svc_contrasts)
  .int <- grep("(Intercept)", colnames(z_i))
  if (length(.int) > 0L && isTRUE(object$svc_omega_is_intercept)) {
    z_i <- z_i[, -.int, drop = FALSE]
  }
  z_i
}

# Epsilon-model covariate, one value per time step in `newdata`.
predict_epsilon_covariate <- function(object, tmb_data, newdata) {
  time_steps <- unique(newdata[[object$time]])
  epsilon_covariate <- rep(0, length(time_steps))
  if (tmb_data$est_epsilon_model) {
    for (i in seq_along(time_steps)) {
      epsilon_covariate[i] <- newdata[newdata[[object$time]] == time_steps[i],
        object$epsilon_predictor, drop = TRUE][[1]]
    }
  }
  epsilon_covariate
}
