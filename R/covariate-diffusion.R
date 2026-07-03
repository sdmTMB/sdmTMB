.as_dgC <- function(x) {
  if (!inherits(x, "sparseMatrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  methods::as(x, "dgCMatrix")
}

.extract_nonlocal_term_exprs <- function(expr) {
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    if (identical(fn, "+")) {
      if (length(expr) == 2L) {
        return(.extract_nonlocal_term_exprs(expr[[2]]))
      }
      if (length(expr) == 3L) {
        return(c(
          .extract_nonlocal_term_exprs(expr[[2]]),
          .extract_nonlocal_term_exprs(expr[[3]])
        ))
      }
    }
    if (identical(fn, "(") && length(expr) == 2L) {
      return(.extract_nonlocal_term_exprs(expr[[2]]))
    }
  }
  list(expr)
}

.append_nonlocal_coef_columns <- function(X, coef_names) {
  if (!length(coef_names)) {
    return(X)
  }
  existing <- intersect(colnames(X), coef_names)
  if (length(existing)) {
    cli_abort(c(
      "Covariate diffusion coefficient names collide with existing fixed-effect columns.",
      "x" = "Conflicting name(s): {.code {paste(existing, collapse = ', ')}}"
    ))
  }
  lag_cols <- matrix(
    0,
    nrow = nrow(X),
    ncol = length(coef_names),
    dimnames = list(NULL, coef_names)
  )
  cbind(X, lag_cols)
}

.parse_nonlocal_formula <- function(nonlocal_formula) {
  if (is.null(nonlocal_formula)) {
    return(NULL)
  }

  if (!inherits(nonlocal_formula, "formula")) {
    cli_abort("`nonlocal_formula` must be `NULL` or a one-sided formula.")
  }

  if (length(nonlocal_formula) != 2L) {
    cli_abort("`nonlocal_formula` must be a one-sided formula such as `~ diffusion(x) + time_lag(x)`.")
  }

  term_exprs <- .extract_nonlocal_term_exprs(nonlocal_formula[[2]])
  if (!length(term_exprs)) {
    cli_abort("`nonlocal_formula` must include at least one lag term.")
  }

  allowed_wrappers <- c("diffusion", "time_lag")

  parsed_terms <- lapply(term_exprs, function(expr) {
    term_label <- paste(deparse(expr), collapse = "")

    if (!is.call(expr)) {
      cli_abort(c(
        "Unsupported term in `nonlocal_formula`.",
        "i" = "Terms must be wrapped in `diffusion()` or `time_lag()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    wrapper <- as.character(expr[[1]])
    if (!wrapper %in% allowed_wrappers) {
      cli_abort(c(
        "Unsupported wrapper in `nonlocal_formula`.",
        "i" = "Allowed wrappers are `diffusion()` and `time_lag()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    if (length(expr) != 2L || !is.symbol(expr[[2]])) {
      cli_abort(c(
        "Unsupported `nonlocal_formula` term structure.",
        "i" = "Use a bare variable name inside each wrapper, e.g. `diffusion(x)`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    variable <- as.character(expr[[2]])
    list(component = wrapper, variable = variable)
  })

  terms_df <- do.call(rbind, lapply(parsed_terms, function(x) {
    data.frame(component = x$component, variable = x$variable, stringsAsFactors = FALSE)
  }))

  duplicated_terms <- duplicated(paste(terms_df$component, terms_df$variable, sep = "::"))
  if (any(duplicated_terms)) {
    dup_labels <- paste0(terms_df$component[duplicated_terms], "(", terms_df$variable[duplicated_terms], ")")
    cli_abort(c(
      "Duplicate `nonlocal_formula` terms are not supported.",
      "x" = "Duplicated term(s): {.code {paste(dup_labels, collapse = ', ')}}"
    ))
  }

  unique_covariates <- unique(terms_df$variable)
  terms_df$covariate_id <- match(terms_df$variable, unique_covariates)
  terms_df$coef_name <- paste0("nl_", terms_df$component, "_", make.names(terms_df$variable))

  list(
    formula = nonlocal_formula,
    terms = terms_df,
    covariates = unique_covariates,
    needs_time = any(terms_df$component == "time_lag")
  )
}

.validate_nonlocal_terms <- function(nonlocal_formula, data, time, multi_family) {
  if (is.null(nonlocal_formula)) {
    return(NULL)
  }

  if (isTRUE(multi_family)) {
    cli_abort("`nonlocal_formula` is currently unsupported for multi-family models.")
  }

  if (isTRUE(nonlocal_formula$needs_time) && is.null(time)) {
    cli_abort(
      "Temporal `nonlocal_formula` terms require a `time` argument."
    )
  }

  missing_covariates <- setdiff(nonlocal_formula$covariates, names(data))
  if (length(missing_covariates)) {
    cli_abort(c(
      "Missing nonlocal covariate(s) in `data`.",
      "x" = "Missing: {.code {paste(missing_covariates, collapse = ', ')}}"
    ))
  }

  non_numeric <- nonlocal_formula$covariates[!vapply(nonlocal_formula$covariates, function(v) {
    is.numeric(data[[v]])
  }, logical(1L))]

  if (length(non_numeric)) {
    cli_abort(c(
      "Covariate diffusion covariates must be numeric.",
      "x" = "Non-numeric covariate(s): {.code {paste(non_numeric, collapse = ', ')}}"
    ))
  }

  nonlocal_formula
}

.build_vertex_time_covariates <- function(covariate_data,
                                          covariates,
                                          A_st,
                                          year_i,
                                          A_spatial_index = NULL,
                                          n_t = NULL) {
  # This is the single vertex-aggregation entry point for covariate diffusion.
  # A future raster-like covariate input can reuse it by supplying its own
  # projection matrix and index vectors.
  if (!inherits(covariate_data, "data.frame")) {
    cli_abort("`covariate_data` must be a data frame.")
  }
  A_st <- .as_dgC(A_st)

  n_obs <- nrow(covariate_data)
  if (length(year_i) != n_obs) {
    cli_abort("`year_i` length must equal `nrow(covariate_data)`.")
  }
  if (is.null(A_spatial_index)) {
    A_spatial_index <- seq_len(n_obs) - 1L
  }
  if (length(A_spatial_index) != n_obs) {
    cli_abort("`A_spatial_index` length must equal `nrow(covariate_data)`.")
  }
  if (!is.numeric(year_i) || anyNA(year_i) || any(!is.finite(year_i))) {
    cli_abort("`year_i` must contain finite numeric indices.")
  }
  if (!is.numeric(A_spatial_index) || anyNA(A_spatial_index) || any(!is.finite(A_spatial_index))) {
    cli_abort("`A_spatial_index` must contain finite numeric indices.")
  }
  if (any(year_i != round(year_i)) || any(A_spatial_index != round(A_spatial_index))) {
    cli_abort("`year_i` and `A_spatial_index` must contain whole-number indices.")
  }

  missing_covariates <- setdiff(covariates, names(covariate_data))
  if (length(missing_covariates)) {
    cli_abort(c(
      "Missing covariate(s) in `covariate_data`.",
      "x" = "Missing: {.code {paste(missing_covariates, collapse = ', ')}}"
    ))
  }

  non_numeric <- covariates[!vapply(covariates, function(v) is.numeric(covariate_data[[v]]), logical(1L))]
  if (length(non_numeric)) {
    cli_abort(c(
      "Covariates supplied to `.build_vertex_time_covariates()` must be numeric.",
      "x" = "Non-numeric covariate(s): {.code {paste(non_numeric, collapse = ', ')}}"
    ))
  }

  year_i <- as.integer(year_i)
  if (!length(year_i)) {
    cli_abort("`year_i` cannot be empty.")
  }
  if (any(year_i < 0L)) {
    cli_abort("`year_i` must use 0-based non-negative indexing.")
  }
  if (is.null(n_t)) {
    n_t <- max(year_i) + 1L
  } else {
    if (!is.numeric(n_t) || length(n_t) != 1L || !is.finite(n_t) || n_t != round(n_t) || n_t <= 0) {
      cli_abort("`n_t` must be a single positive integer.")
    }
    n_t <- as.integer(n_t)
    if (max(year_i) >= n_t) {
      cli_abort("`year_i` contains a time index >= `n_t`.")
    }
  }

  A_spatial_index <- as.integer(A_spatial_index) + 1L
  if (any(A_spatial_index < 1L | A_spatial_index > nrow(A_st))) {
    cli_abort("`A_spatial_index` contains indices outside valid range.")
  }

  n_vertices <- ncol(A_st)
  n_covariates <- length(covariates)
  out <- array(
    0,
    dim = c(n_vertices, n_t, n_covariates),
    dimnames = list(NULL, NULL, covariates)
  )

  for (cov_idx in seq_along(covariates)) {
    cov_name <- covariates[[cov_idx]]
    x <- covariate_data[[cov_name]]
    if (any(is.infinite(x), na.rm = TRUE)) {
      cli_abort("Covariate `{cov_name}` contains Inf/-Inf values.")
    }
    for (t_i in seq_len(n_t)) {
      obs_this_time <- which(year_i == (t_i - 1L))
      if (!length(obs_this_time)) {
        next
      }
      x_t <- x[obs_this_time]
      keep <- !is.na(x_t)
      if (!any(keep)) {
        next
      }
      A_t <- A_st[A_spatial_index[obs_this_time[keep]], , drop = FALSE]
      x_t <- x_t[keep]
      numerator <- as.vector(Matrix::crossprod(A_t, x_t))
      denominator <- as.vector(Matrix::crossprod(A_t, rep(1, length(x_t))))
      good <- denominator > 0
      if (any(good)) {
        out[good, t_i, cov_idx] <- numerator[good] / denominator[good]
      }
    }
  }

  list(
    covariate_vertex_time = out,
    covariates = covariates,
    n_t = n_t,
    n_vertices = n_vertices
  )
}

.prepare_nonlocal_grid_inputs <- function(grid,
                                                     nonlocal_formula,
                                                     mesh,
                                                     xy_cols,
                                                     time,
                                                     time_df,
                                                     full_time_vec) {
  if (!inherits(grid, "data.frame")) {
    cli_abort("The nonlocal grid data must be `NULL` or a data frame.")
  }
  if (is.null(xy_cols) || length(xy_cols) != 2L) {
    cli_abort("The nonlocal grid data requires a mesh built with known `xy_cols` (e.g., from `make_mesh()`).")
  }

  missing_xy <- setdiff(xy_cols, names(grid))
  if (length(missing_xy)) {
    cli_abort(c(
      "The nonlocal grid data is missing required coordinate column(s).",
      "x" = "Missing: {.code {paste(missing_xy, collapse = ', ')}}"
    ))
  }
  non_numeric_xy <- xy_cols[!vapply(xy_cols, function(col) is.numeric(grid[[col]]), logical(1L))]
  if (length(non_numeric_xy)) {
    cli_abort(c(
      "The nonlocal grid data coordinates must be numeric.",
      "x" = "Non-numeric coordinate column(s): {.code {paste(non_numeric_xy, collapse = ', ')}}"
    ))
  }
  invalid_xy <- xy_cols[!vapply(xy_cols, function(col) all(is.finite(grid[[col]])), logical(1L))]
  if (length(invalid_xy)) {
    cli_abort(c(
      "The nonlocal grid data coordinates must be finite and cannot contain `NA` values.",
      "x" = "Invalid coordinate column(s): {.code {paste(invalid_xy, collapse = ', ')}}"
    ))
  }

  missing_covariates <- setdiff(nonlocal_formula$covariates, names(grid))
  if (length(missing_covariates)) {
    cli_abort(c(
      "The nonlocal grid data is missing required covariate column(s).",
      "x" = "Missing: {.code {paste(missing_covariates, collapse = ', ')}}"
    ))
  }
  non_numeric <- nonlocal_formula$covariates[!vapply(nonlocal_formula$covariates, function(v) {
    is.numeric(grid[[v]])
  }, logical(1L))]
  if (length(non_numeric)) {
    cli_abort(c(
      "The nonlocal grid data covariates must be numeric.",
      "x" = "Non-numeric covariate(s): {.code {paste(non_numeric, collapse = ', ')}}"
    ))
  }

  if (isTRUE(nonlocal_formula$needs_time)) {
    if (!time %in% names(grid)) {
      cli_abort("The nonlocal grid data is missing the time column {.code {time}}.")
    }
    missing_slices <- setdiff(full_time_vec, grid[[time]])
    if (length(missing_slices)) {
      cli_abort(c(
        "The nonlocal grid data does not cover all fitted (+ `extra_time`) time slices.",
        "x" = "Missing time slice(s): {.code {paste(missing_slices, collapse = ', ')}}"
      ))
    }
    year_i <- time_df$year_i[match(grid[[time]], time_df$time_from_data)]
  } else {
    year_i <- rep(0L, nrow(grid))
  }

  A_st <- fmesher::fm_basis(mesh, loc = as.matrix(grid[, xy_cols, drop = FALSE]))
  A_spatial_index <- seq_len(nrow(grid)) - 1L

  list(
    data = grid,
    A_st = A_st,
    A_spatial_index = A_spatial_index,
    year_i = year_i,
    n_t = nrow(time_df)
  )
}

.build_nonlocal_tmb_data <- function(nonlocal_formula,
                                                data,
                                                A_st,
                                                A_spatial_index,
                                                year_i,
                                                n_t) {
  if (is.null(nonlocal_formula)) {
    return(NULL)
  }

  vertex_cov <- .build_vertex_time_covariates(
    covariate_data = data,
    covariates = nonlocal_formula$covariates,
    A_st = A_st,
    year_i = year_i,
    A_spatial_index = A_spatial_index,
    n_t = n_t
  )
  covariate_vertex_time <- vertex_cov$covariate_vertex_time

  component_levels <- c("diffusion", "time_lag")
  component_id <- match(nonlocal_formula$terms$component, component_levels)
  terms_df <- nonlocal_formula$terms
  covariates <- nonlocal_formula$covariates
  covariate_has_spatial <- integer(length(covariates))
  covariate_has_temporal <- integer(length(covariates))
  for (i in seq_along(covariates)) {
    components <- terms_df$component[terms_df$variable == covariates[i]]
    covariate_has_spatial[i] <- any(components == "diffusion")
    covariate_has_temporal[i] <- any(components == "time_lag")
  }

  list(
    covariate_vertex_time = covariate_vertex_time,
    covariates = covariates,
    covariate_has_spatial = covariate_has_spatial,
    covariate_has_temporal = covariate_has_temporal,
    term_component = nonlocal_formula$terms$component,
    term_component_id = as.integer(component_id),
    term_covariate_index = as.integer(nonlocal_formula$terms$covariate_id),
    term_covariate_index0 = as.integer(nonlocal_formula$terms$covariate_id - 1L),
    term_coef_name = nonlocal_formula$terms$coef_name,
    n_vertices = vertex_cov$n_vertices,
    n_t = vertex_cov$n_t,
    n_covariates = length(covariates),
    n_terms = nrow(nonlocal_formula$terms)
  )
}

.solve_nonlocal_vertex_time <- function(component, vertex_time_input, M0, M1,
                                                  kappaS, kappaT,
                                                  has_space = NULL,
                                                  has_time = NULL) {
  n_vertices <- nrow(vertex_time_input)
  n_t <- ncol(vertex_time_input)
  out <- matrix(0, nrow = n_vertices, ncol = n_t)

  solve_sparse <- function(A, rhs, label) {
    tryCatch(
      as.numeric(Matrix::solve(A, rhs)),
      error = function(e) {
        cli_abort(c(
          "Covariate diffusion sparse solve failed in diagnostic plotting.",
          "x" = paste0(label, ": ", conditionMessage(e))
        ))
      }
    )
  }

  if (component == "combined") {
    has_space <- isTRUE(has_space)
    has_time <- isTRUE(has_time)
    if (has_space && !has_time) {
      component <- "diffusion"
    } else if (has_time && !has_space) {
      component <- "time_lag"
    } else if (!has_space && !has_time) {
      cli_abort("`component = \"combined\"` requires at least one fitted `diffusion()` or `time_lag()` term.")
    } else {
      kappaS_scale <- 1 / (kappaS^2)
      kappaT_scale <- kappaT
      system_mat <- M0 + kappaS_scale * M1
      for (tt in seq_len(n_t)) {
        rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE])
        if (tt > 1L && kappaT_scale != 0) {
          rhs <- rhs + as.numeric(kappaT_scale * (M0 %*% out[, tt - 1L, drop = TRUE]))
        }
        out[, tt] <- solve_sparse(system_mat, rhs, "combined system (space + time)")
      }
      return(out)
    }
  }

  if (component == "diffusion") {
    kappaS_scale <- 1 / (kappaS^2)
    system_mat <- M0 + kappaS_scale * M1
    for (tt in seq_len(n_t)) {
      rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE])
      out[, tt] <- solve_sparse(system_mat, rhs, "spatial system (M0 + kappaS^-2 * M1)")
    }
    return(out)
  }

  if (component == "time_lag") {
    denom <- 1 + kappaT
    out[, 1L] <- vertex_time_input[, 1L] / denom
    if (n_t > 1L) {
      for (tt in 2:n_t) {
        out[, tt] <- (vertex_time_input[, tt] + kappaT * out[, tt - 1L]) / denom
      }
    }
    return(out)
  }

  cli_abort("Unknown covariate-diffusion component in solver.")
}

.project_nonlocal_vertex_time <- function(transformed_vertex_time,
                                                     A_st,
                                                     A_spatial_index,
                                                     year_i,
                                                     n_t) {
  A_st <- .as_dgC(A_st)
  n_i <- length(A_spatial_index)
  A_spatial_index <- as.integer(A_spatial_index) + 1L
  year_i <- as.integer(year_i)
  projected_by_t <- lapply(seq_len(n_t), function(tt) {
    as.numeric(A_st %*% transformed_vertex_time[, tt, drop = TRUE])
  })
  out <- numeric(n_i)
  for (i in seq_len(n_i)) {
    out[i] <- projected_by_t[[year_i[i] + 1L]][A_spatial_index[i]]
  }
  out
}

.nonlocal_predict_colnames <- function(term_coef_name) {
  term_coef_name
}

.compute_nonlocal_term_values <- function(nonlocal_parsed,
                                                 covariate_vertex_time,
                                                 A_st,
                                                 A_spatial_index,
                                                 year_i,
                                                 n_t,
                                                 M0,
                                                 M1,
                                                 log_kappaS_nl,
                                                 kappaT_nl_raw) {
  if (is.null(nonlocal_parsed)) {
    return(NULL)
  }
  n_terms <- nonlocal_parsed$n_terms
  if (is.null(n_terms) || n_terms < 1L) {
    return(NULL)
  }
  n_covariates <- nonlocal_parsed$n_covariates
  if (length(log_kappaS_nl) != n_covariates ||
      length(kappaT_nl_raw) != n_covariates) {
    cli_abort("Covariate diffusion parameter vectors did not match the expected number of lag covariates.")
  }

  if (length(dim(covariate_vertex_time)) != 3L) {
    cli_abort("`covariate_vertex_time` must be a 3D array: [vertices, time, covariates].")
  }
  if (dim(covariate_vertex_time)[3] != n_covariates) {
    cli_abort("`covariate_vertex_time` covariate dimension does not match covariate-diffusion metadata.")
  }

  term_out <- matrix(0, nrow = length(A_spatial_index), ncol = n_terms)
  colnames(term_out) <- .nonlocal_predict_colnames(nonlocal_parsed$term_coef_name)

  for (term_i in seq_len(n_terms)) {
    component <- nonlocal_parsed$term_component[[term_i]]
    cov_i <- nonlocal_parsed$term_covariate_index[[term_i]]
    cov_slice <- matrix(
      covariate_vertex_time[, , cov_i],
      nrow = dim(covariate_vertex_time)[1],
      ncol = dim(covariate_vertex_time)[2]
    )
    kappaS <- exp(log_kappaS_nl[[cov_i]])
    kappaT <- kappaT_nl_raw[[cov_i]]
    transformed_vertex_time <- .solve_nonlocal_vertex_time(
      component = component,
      vertex_time_input = cov_slice,
      M0 = M0,
      M1 = M1,
      kappaS = kappaS,
      kappaT = kappaT
    )
    term_out[, term_i] <- .project_nonlocal_vertex_time(
      transformed_vertex_time = transformed_vertex_time,
      A_st = A_st,
      A_spatial_index = A_spatial_index,
      year_i = year_i,
      n_t = n_t
    )
  }

  term_out
}

.nl_plot_extract_mesh <- function(mesh) {
  if (is.null(mesh$loc) || is.null(mesh$graph) || is.null(mesh$graph$tv)) {
    cli_abort("Could not find mesh vertices/triangles in `object$spde$mesh`.")
  }
  loc <- as.matrix(mesh$loc[, 1:2, drop = FALSE])
  tv <- as.matrix(mesh$graph$tv)
  if (!nrow(tv) || ncol(tv) != 3L) {
    cli_abort("Mesh triangle connectivity (`mesh$graph$tv`) must be an n x 3 matrix.")
  }
  if (min(tv) == 0L) tv <- tv + 1L
  if (any(tv < 1L) || any(tv > nrow(loc))) {
    cli_abort("Mesh triangle indices were out of bounds for mesh vertices.")
  }
  center <- colMeans(loc)
  vertex_i <- which.min(rowSums((loc - matrix(center, nrow = nrow(loc), ncol = 2L, byrow = TRUE))^2))
  list(loc = loc, tv = tv, vertex_i = vertex_i)
}

.nl_plot_resolve_time <- function(object, component, time_value, n_steps) {
  n_t <- object$tmb_data$n_t
  if (is.null(n_t) || !length(n_t) || n_t < 1L) {
    cli_abort("Could not determine the number of time slices from `object$tmb_data$n_t`.")
  }
  if (!is.null(object$time_lu) &&
      "time_from_data" %in% names(object$time_lu) &&
      nrow(object$time_lu) == n_t) {
    time_values <- object$time_lu$time_from_data
  } else {
    time_values <- seq_len(n_t)
  }
  if (is.null(time_value)) {
    time_i <- 1L
  } else {
    time_i <- match(time_value, time_values)
    if (is.na(time_i) && is.numeric(time_value) && length(time_value) == 1L) {
      candidate <- as.integer(round(time_value))
      if (is.finite(candidate) && candidate >= 1L && candidate <= n_t) {
        time_i <- candidate
      }
    }
    if (is.na(time_i)) {
      preview <- paste(utils::head(time_values, 12L), collapse = ", ")
      cli_abort(c(
        "Could not match `time_value` to modeled time slices.",
        "x" = "Available times include: {.code {preview}}."
      ))
    }
  }
  if (component == "diffusion") {
    time_idx <- time_i
  } else {
    time_idx <- seq.int(time_i, min(n_t, time_i + n_steps - 1L))
    if (length(time_idx) < n_steps) {
      cli_inform("Requested `n_steps` exceeded modeled time range; using available trailing slices.")
    }
  }
  list(n_t = n_t, time_values = time_values, time_i = time_i, time_idx = time_idx)
}

.nl_plot_extract_kappas <- function(object, cov_i) {
  if (!is.null(object$model) &&
      !is.null(object$model$par) &&
      !is.null(object$tmb_obj) &&
      !is.null(object$tmb_obj$env)) {
    params <- object$tmb_obj$env$parList(object$model$par)
  } else if (!is.null(object$tmb_params)) {
    params <- object$tmb_params
  } else {
    cli_abort("Could not extract covariate-diffusion parameters from `object`.")
  }
  list(
    kappaS = exp(params$log_kappaS_nl[cov_i]),
    kappaT = params$kappaT_nl_raw[cov_i]
  )
}

.nl_plot_context <- function(object, covariate, component, component_missing,
                             time_value, n_steps, function_name) {
  stopifnot(inherits(object, "sdmTMB"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli_abort("`ggplot2` must be installed to use `{function_name}()`.")
  }
  if (isTRUE(component_missing)) {
    cli_abort("`component` is required and must be one of `diffusion`, `time_lag`, or `combined`.")
  }
  if (!component %in% c("diffusion", "time_lag", "combined")) {
    cli_abort("`component` must be exactly one of `diffusion`, `time_lag`, or `combined`.")
  }
  if (is.null(object$nonlocal_parsed)) {
    cli_abort("`object` does not contain `nonlocal_parsed`.")
  }
  if (!is.numeric(n_steps) || length(n_steps) != 1L || !is.finite(n_steps) || n_steps < 1L) {
    cli_abort("`n_steps` must be a single positive integer.")
  }
  n_steps <- as.integer(round(n_steps))
  if (is.null(object$nonlocal_formula_parsed) ||
      is.null(object$nonlocal_formula_parsed$terms) ||
      nrow(object$nonlocal_formula_parsed$terms) == 0L) {
    cli_abort("`object` does not contain covariate-diffusion terms.")
  }

  terms_df <- object$nonlocal_formula_parsed$terms
  covariates <- unique(terms_df$variable)
  if (is.null(covariate)) {
    if (length(covariates) != 1L) {
      cli_abort(c(
        "Multiple covariate-diffusion covariates are present.",
        "x" = "Set `covariate` to one of: {.code {paste(covariates, collapse = ', ')}}."
      ))
    }
    covariate <- covariates[[1L]]
  }
  covariate <- as.character(covariate[[1L]])
  if (!covariate %in% covariates) {
    cli_abort(c(
      "Unknown covariate-diffusion covariate.",
      "x" = "Could not find `{covariate}` in the fitted covariate-diffusion terms."
    ))
  }
  components_for_covariate <- unique(terms_df$component[terms_df$variable == covariate])
  if (component != "combined" && !component %in% components_for_covariate) {
    cli_abort(c(
      "Requested component/covariate term was not fitted.",
      "x" = "No term `{component}({covariate})` in `object$nonlocal_formula`."
    ))
  }
  component_for_time <- if (component == "combined" && !"time_lag" %in% components_for_covariate) "diffusion" else component

  mesh_info <- .nl_plot_extract_mesh(object$spde$mesh)
  time_info <- .nl_plot_resolve_time(
    object, component_for_time, time_value, n_steps
  )

  cov_i <- match(covariate, object$nonlocal_parsed$covariates)
  if (is.na(cov_i)) {
    cli_abort("Internal mismatch: selected covariate not found in `nonlocal_parsed$covariates`.")
  }

  xy_cols <- object$spde$xy_cols
  list(
    covariate = covariate,
    component = component,
    mesh_info = mesh_info,
    time_info = time_info,
    cov_i = cov_i,
    kappas = .nl_plot_extract_kappas(object, cov_i),
    has_space = as.logical(object$nonlocal_parsed$covariate_has_spatial[cov_i]),
    has_time = as.logical(object$nonlocal_parsed$covariate_has_temporal[cov_i]),
    xlab = if (!is.null(xy_cols) && length(xy_cols) >= 2L) xy_cols[1] else "x",
    ylab = if (!is.null(xy_cols) && length(xy_cols) >= 2L) xy_cols[2] else "y"
  )
}

.nl_plot_time_panels <- function(first_field, first_title, transformed_vertex_time,
                                 time_i, time_idx, time_values) {
  panel_fields <- vector("list", length(time_idx) + 1L)
  panel_titles <- character(length(panel_fields))
  panel_fields[[1L]] <- first_field
  panel_titles[[1L]] <- first_title
  for (j in seq_along(time_idx)) {
    tt <- time_idx[j]
    lag <- tt - time_i
    panel_fields[[j + 1L]] <- transformed_vertex_time[, tt]
    panel_titles[[j + 1L]] <- if (lag == 0L) {
      paste0("diffused (t=", time_values[tt], ")")
    } else {
      paste0("lag+", lag, " (t=", time_values[tt], ")")
    }
  }
  list(fields = panel_fields, titles = panel_titles)
}

.nl_plot_locations <- function(object, newdata, type) {
  mesh <- object$spde$mesh
  xy_cols <- object$spde$xy_cols
  mesh_info <- .nl_plot_extract_mesh(mesh)
  if (is.null(newdata)) {
    if (identical(type, "raster")) {
      cli_abort('`type = "raster"` requires `newdata`.')
    }
    loc <- mesh_info$loc
    tv <- mesh_info$tv
    edge_i <- rbind(tv[, c(1L, 2L)], tv[, c(2L, 3L)], tv[, c(3L, 1L)])
    edge_i <- t(apply(edge_i, 1L, sort))
    edge_i <- unique(edge_i)
    edge_df <- data.frame(
      x = loc[edge_i[, 1L], 1],
      y = loc[edge_i[, 1L], 2],
      xend = loc[edge_i[, 2L], 1],
      yend = loc[edge_i[, 2L], 2]
    )
    return(list(
      loc = loc,
      A = Matrix::Diagonal(nrow(loc)),
      edge_df = edge_df
    ))
  }

  if (is.null(xy_cols) || length(xy_cols) != 2L) {
    cli_abort("`newdata` requires a mesh built with known `xy_cols` (e.g., from `make_mesh()`).")
  }
  if (!inherits(newdata, "data.frame") || any(!xy_cols %in% names(newdata))) {
    cli_abort(c(
      "`newdata` must be a data frame with coordinate columns matching the fitted mesh.",
      "x" = "Required columns: {.code {paste(xy_cols, collapse = ', ')}}."
    ))
  }
  loc <- unique(as.data.frame(newdata[, xy_cols, drop = FALSE]))
  loc <- as.matrix(loc)
  if (!is.numeric(loc) || any(!is.finite(loc))) {
    cli_abort("`newdata` coordinate columns must be finite numeric values.")
  }
  list(
    loc = loc,
    A = fmesher::fm_basis(mesh, loc = loc),
    edge_df = NULL
  )
}

.nl_plot_df <- function(loc, A, panel_fields, panel_titles, common_scale) {
  plot_values <- vapply(panel_fields, function(v) {
    as.numeric(A %*% v)
  }, numeric(nrow(loc)))
  if (!is.matrix(plot_values)) {
    plot_values <- matrix(plot_values, ncol = 1L)
  }
  colnames(plot_values) <- panel_titles

  plot_df <- do.call(rbind, lapply(seq_len(ncol(plot_values)), function(j) {
    data.frame(
      panel = panel_titles[j],
      x = loc[, 1],
      y = loc[, 2],
      value = plot_values[, j],
      stringsAsFactors = FALSE
    )
  }))
  plot_df$panel <- factor(plot_df$panel, levels = panel_titles)
  plot_df$value_plot <- plot_df$value

  fill_name <- "Value"
  if (!isTRUE(common_scale)) {
    fill_name <- "Relative value"
    for (p in levels(plot_df$panel)) {
      i <- which(plot_df$panel == p)
      rng <- range(plot_df$value_plot[i], finite = TRUE)
      if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        plot_df$value_plot[i] <- 0
      } else {
        plot_df$value_plot[i] <- (plot_df$value_plot[i] - rng[1]) / (rng[2] - rng[1])
      }
    }
  }
  fill_limits <- if (isTRUE(common_scale)) NULL else c(0, 1)

  list(
    plot_df = plot_df,
    fill_name = fill_name,
    fill_limits = fill_limits
  )
}

.nl_plot_ggplot <- function(plot_df, edge_df, type, fill_limits, fill_name, xlim, ylim, xlab, ylab,
                            scale = c("distiller", "viridis")) {
  scale <- match.arg(scale)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$y))
  if (!is.null(edge_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = edge_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        inherit.aes = FALSE,
        colour = "grey85",
        linewidth = 0.2
      )
  }
  if (identical(type, "raster")) {
    p <- p + ggplot2::geom_raster(ggplot2::aes(fill = .data$value_plot))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$value_plot))
  }
  p <- p +
    ggplot2::facet_wrap(stats::as.formula("~ panel"), nrow = 1L) +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = identical(type, "point"))
  p <- p +
    if (scale == "viridis") {
      ggplot2::scale_colour_viridis_c(
        limits = fill_limits,
        name = fill_name,
        aesthetics = c("colour", "fill")
      )
    } else {
      ggplot2::scale_colour_distiller(
        limits = fill_limits, palette = "Blues",
        name = fill_name, direction = 1,
        aesthetics = c("colour", "fill")
      )
    }
  p +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab)
}

#' Plot Covariate-Diffusion Diagnostics
#'
#' Visualize fitted covariate-diffusion transforms or impulse-response kernels
#' for one selected covariate-diffusion term.
#' By default, values are plotted at mesh vertices with the mesh edges shown in
#' light grey. Values can also be evaluated at supplied `newdata` coordinates
#' and plotted as points or a raster.
#'
#' @param object A fitted [sdmTMB()] model with `nonlocal_formula`.
#' @param covariate Optional covariate name from `nonlocal_formula`.
#'   Required when multiple lag covariates were fitted.
#' @param component Covariate-diffusion component name. Must be one of
#'   `"diffusion"`, `"time_lag"`, or `"combined"`. `"combined"`
#'   plots the joint response of all covariate-diffusion components fitted for
#'   `covariate`.
#' @param time_value Optional time slice to plot or use for the impulse. Supply
#'   either a modeled time value or a 1-based time index. Defaults to 1.
#' @param n_steps Number of transformed slices to plot starting at
#'   `time_value`.
#' @param common_scale Should plotted panels share a common color scale?
#'   Defaults to `TRUE` for `plot_nonlocal_covariate()` and `FALSE` for
#'   `plot_nonlocal_kernel()`. `component = "time_lag"` alone likely needs
#'   `common_scale = TRUE` to make sense.
#' @param newdata Optional data frame with x/y coordinate columns matching the
#'   fitted mesh. If supplied, values are evaluated at the unique `newdata`
#'   coordinates. If `NULL`, values are evaluated at mesh vertices.
#' @param type Plot type: `"point"` or `"raster"`. `"raster"` requires
#'   `newdata`.
#'
#' @details
#' `plot_nonlocal_covariate()` visualizes the original mesh-vertex covariate
#' field and its fitted covariate-diffusion transform of one selected covariate
#' time slice over one or more lagged output time slices.
#'
#' `plot_nonlocal_kernel()` visualizes an impulse covariate diffusing through
#' one covariate-diffusion component.
#'
#' @return A `ggplot` object.
#'
#' @examplesIf ggplot2_installed()
#'
#' # Simulate some data for fitting:
#' set.seed(1)
#' n_t <- 6
#' n_sites <- 80
#' sites <- data.frame(X = runif(n_sites), Y = runif(n_sites))
#' dat <- data.frame(
#'   X = rep(sites$X, times = n_t),
#'   Y = rep(sites$Y, times = n_t),
#'   year = rep(seq_len(n_t), each = n_sites)
#' )
#' dat$x1 <- as.numeric(scale(
#'   sin(2 * pi * (dat$X + dat$year / 6)) +
#'     cos(2 * pi * (dat$Y - dat$year / 8)) +
#'     0.4 * sin(4 * pi * dat$X) * cos(dat$year / 2) +
#'     rnorm(nrow(dat), sd = 0.15)
#' ))
#' mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.12)
#' sim <- simulate_new(
#'   formula = ~ 1,
#'   data = dat,
#'   mesh = mesh,
#'   time = "year",
#'   family = gaussian(),
#'   spatial = "off",
#'   spatiotemporal = "off",
#'   range = 0.3,
#'   sigma_O = 0,
#'   sigma_E = 0,
#'   phi = 0.1,
#'   B = c(0, 0.7, 0.6),
#'   nonlocal_formula = ~ diffusion(x1) + time_lag(x1),
#'   lags_kappaS = 4.4,
#'   lags_rhoT = 0.3,
#'   seed = 123
#' )
#' dat$observed <- sim$observed
#'
#' # Fit the model:
#' fit <- sdmTMB(
#'   observed ~ 1,
#'   data = dat,
#'   mesh = mesh,
#'   time = "year",
#'   spatial = "off", # keeping example simple
#'   spatiotemporal = "off", # keeping example simple
#'   family = gaussian(),
#'   nonlocal_formula = ~ diffusion(x1) + time_lag(x1) #<
#' )
#'
#' plot_nonlocal_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "diffusion"
#' )
#' plot_nonlocal_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "time_lag",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' plot_nonlocal_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "combined",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' plot_nonlocal_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "diffusion"
#' )
#' plot_nonlocal_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "time_lag",
#'   time_value = 1,
#'   n_steps = 2,
#'   common_scale = TRUE #<
#' )
#' plot_nonlocal_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "combined",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' @rdname nonlocal_formula_plots
#' @export
plot_nonlocal_covariate <- function(object,
                                    component,
                                    newdata = NULL,
                                    type = c("point", "raster"),
                                    covariate = NULL,
                                    time_value = 1,
                                    n_steps = 1L,
                                    common_scale = TRUE) {
  type <- match.arg(type)
  ctx <- .nl_plot_context(
    object = object,
    covariate = covariate,
    component = if (missing(component)) NULL else component,
    component_missing = missing(component),
    time_value = time_value,
    n_steps = n_steps,
    function_name = "plot_nonlocal_covariate"
  )
  plot_locations <- .nl_plot_locations(object, newdata, type)

  original_vertex_time <- matrix(
    object$nonlocal_parsed$covariate_vertex_time[, , ctx$cov_i],
    nrow = object$nonlocal_parsed$n_vertices,
    ncol = ctx$time_info$n_t
  )

  source_vertex_time <- matrix(0, nrow = nrow(original_vertex_time), ncol = ctx$time_info$n_t)
  source_vertex_time[, ctx$time_info$time_i] <- original_vertex_time[, ctx$time_info$time_i]

  transformed_vertex_time <- .solve_nonlocal_vertex_time(
    component = ctx$component,
    vertex_time_input = source_vertex_time,
    M0 = object$tmb_data$spde$M0,
    M1 = object$tmb_data$spde$M1,
    kappaS = ctx$kappas$kappaS,
    kappaT = ctx$kappas$kappaT,
    has_space = ctx$has_space,
    has_time = ctx$has_time
  )

  time_values <- ctx$time_info$time_values
  time_i <- ctx$time_info$time_i
  time_idx <- ctx$time_info$time_idx
  panels <- .nl_plot_time_panels(
    first_field = original_vertex_time[, time_i],
    first_title = paste0("original (t=", time_values[time_i], ")"),
    transformed_vertex_time = transformed_vertex_time,
    time_i = time_i,
    time_idx = time_idx,
    time_values = time_values
  )
  panel <- .nl_plot_df(
    loc = plot_locations$loc,
    A = plot_locations$A,
    panel_fields = panels$fields,
    panel_titles = panels$titles,
    common_scale = common_scale
  )

  .nl_plot_ggplot(
    plot_df = panel$plot_df,
    edge_df = plot_locations$edge_df,
    type = type,
    fill_limits = panel$fill_limits,
    fill_name = panel$fill_name,
    xlim = range(plot_locations$loc[, 1]),
    ylim = range(plot_locations$loc[, 2]),
    xlab = ctx$xlab,
    ylab = ctx$ylab,
    scale = "viridis"
  )
}

#' @rdname nonlocal_formula_plots
#' @export
plot_nonlocal_kernel <- function(object,
                                  component,
                                  newdata = NULL,
                                  type = c("point", "raster"),
                                  covariate = NULL,
                                  time_value = 1,
                                  n_steps = 3L,
                                  common_scale = FALSE) {
  type <- match.arg(type)
  ctx <- .nl_plot_context(
    object = object,
    covariate = covariate,
    component = if (missing(component)) NULL else component,
    component_missing = missing(component),
    time_value = time_value,
    n_steps = n_steps,
    function_name = "plot_nonlocal_kernel"
  )
  plot_locations <- .nl_plot_locations(object, newdata, type)

  time_values <- ctx$time_info$time_values
  n_t <- ctx$time_info$n_t
  time_i <- ctx$time_info$time_i
  time_idx <- ctx$time_info$time_idx
  vertex_i <- ctx$mesh_info$vertex_i
  n_vertices <- nrow(ctx$mesh_info$loc)
  impulse_vertex_time <- matrix(0, nrow = n_vertices, ncol = n_t)
  impulse_vertex_time[vertex_i, time_i] <- 1

  transformed_vertex_time <- .solve_nonlocal_vertex_time(
    component = ctx$component,
    vertex_time_input = impulse_vertex_time,
    M0 = object$tmb_data$spde$M0,
    M1 = object$tmb_data$spde$M1,
    kappaS = ctx$kappas$kappaS,
    kappaT = ctx$kappas$kappaT,
    has_space = ctx$has_space,
    has_time = ctx$has_time
  )

  panels <- .nl_plot_time_panels(
    first_field = impulse_vertex_time[, time_i],
    first_title = paste0("original (t=", time_values[time_i], ")"),
    transformed_vertex_time = transformed_vertex_time,
    time_i = time_i,
    time_idx = time_idx,
    time_values = time_values
  )
  panel <- .nl_plot_df(
    loc = plot_locations$loc,
    A = plot_locations$A,
    panel_fields = panels$fields,
    panel_titles = panels$titles,
    common_scale = common_scale
  )

  .nl_plot_ggplot(
    plot_df = panel$plot_df,
    edge_df = plot_locations$edge_df,
    type = type,
    fill_limits = panel$fill_limits,
    fill_name = panel$fill_name,
    xlim = range(plot_locations$loc[, 1]),
    ylim = range(plot_locations$loc[, 2]),
    xlab = ctx$xlab,
    ylab = ctx$ylab
  )
}
