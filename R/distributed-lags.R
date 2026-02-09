.extract_distributed_lag_term_exprs <- function(expr) {
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    if (identical(fn, "+")) {
      if (length(expr) == 2L) {
        return(.extract_distributed_lag_term_exprs(expr[[2]]))
      }
      if (length(expr) == 3L) {
        return(c(
          .extract_distributed_lag_term_exprs(expr[[2]]),
          .extract_distributed_lag_term_exprs(expr[[3]])
        ))
      }
    }
    if (identical(fn, "(") && length(expr) == 2L) {
      return(.extract_distributed_lag_term_exprs(expr[[2]]))
    }
  }
  list(expr)
}

.append_distributed_lag_coef_columns <- function(X, coef_names) {
  if (!length(coef_names)) {
    return(X)
  }
  existing <- intersect(colnames(X), coef_names)
  if (length(existing)) {
    cli_abort(c(
      "Distributed lag coefficient names collide with existing fixed-effect columns.",
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

.coerce_integerish <- function(x, name) {
  if (!is.numeric(x)) {
    cli_abort("`{name}` must be numeric and integer-valued.")
  }
  if (anyNA(x)) {
    cli_abort("`{name}` cannot contain `NA` values.")
  }
  if (!length(x)) {
    return(integer(0))
  }
  if (any(!is.finite(x))) {
    cli_abort("`{name}` cannot contain non-finite values.")
  }
  tol <- sqrt(.Machine$double.eps)
  rounded <- round(x)
  if (any(abs(x - rounded) > tol)) {
    cli_abort("`{name}` must contain whole-number indices.")
  }
  as.integer(rounded)
}

.parse_distributed_lags_formula <- function(distributed_lags) {
  if (is.null(distributed_lags)) {
    return(NULL)
  }

  if (!inherits(distributed_lags, "formula")) {
    cli_abort("`distributed_lags` must be `NULL` or a one-sided formula.")
  }

  if (length(distributed_lags) != 2L) {
    cli_abort("`distributed_lags` must be a one-sided formula such as `~ space(x) + time(x)`.")
  }

  term_exprs <- .extract_distributed_lag_term_exprs(distributed_lags[[2]])
  if (!length(term_exprs)) {
    cli_abort("`distributed_lags` must include at least one lag term.")
  }

  allowed_wrappers <- c("space", "time", "spacetime")

  parsed_terms <- lapply(term_exprs, function(expr) {
    term_label <- paste(deparse(expr), collapse = "")

    if (!is.call(expr)) {
      cli_abort(c(
        "Unsupported term in `distributed_lags`.",
        "i" = "Terms must be wrapped in `space()`, `time()`, or `spacetime()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    wrapper <- as.character(expr[[1]])
    if (!wrapper %in% allowed_wrappers) {
      cli_abort(c(
        "Unsupported wrapper in `distributed_lags`.",
        "i" = "Allowed wrappers are `space()`, `time()`, and `spacetime()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    if (length(expr) != 2L || !is.symbol(expr[[2]])) {
      cli_abort(c(
        "Unsupported `distributed_lags` term structure.",
        "i" = "Use a bare variable name inside each wrapper, e.g. `space(depth)`.",
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
      "Duplicate `distributed_lags` terms are not supported.",
      "x" = "Duplicated term(s): {.code {paste(dup_labels, collapse = ', ')}}"
    ))
  }

  unique_covariates <- unique(terms_df$variable)
  terms_df$covariate_id <- match(terms_df$variable, unique_covariates)
  terms_df$coef_name <- paste0("dl_", terms_df$component, "_", make.names(terms_df$variable))

  list(
    formula = distributed_lags,
    terms = terms_df,
    covariates = unique_covariates,
    needs_time = any(terms_df$component %in% c("time", "spacetime"))
  )
}

.validate_distributed_lag_terms <- function(distributed_lags, data, time, delta, multi_family) {
  if (is.null(distributed_lags)) {
    return(NULL)
  }

  if (isTRUE(multi_family)) {
    cli_abort("`distributed_lags` is currently unsupported for multi-family models.")
  }

  if (isTRUE(distributed_lags$needs_time) && is.null(time)) {
    cli_abort(
      "`distributed_lags` terms wrapped in `time()` or `spacetime()` require a `time` argument."
    )
  }

  missing_covariates <- setdiff(distributed_lags$covariates, names(data))
  if (length(missing_covariates)) {
    cli_abort(c(
      "Missing distributed lag covariate(s) in `data`.",
      "x" = "Missing: {.code {paste(missing_covariates, collapse = ', ')}}"
    ))
  }

  non_numeric <- distributed_lags$covariates[!vapply(distributed_lags$covariates, function(v) {
    is.numeric(data[[v]])
  }, logical(1L))]

  if (length(non_numeric)) {
    cli_abort(c(
      "Distributed lag covariates must be numeric.",
      "x" = "Non-numeric covariate(s): {.code {paste(non_numeric, collapse = ', ')}}"
    ))
  }

  distributed_lags
}

# Coerce integer-valued vector to 0-based indexing.
# Accepts either 0-based (min == 0) or 1-based (min == 1) input.
.to_zero_based <- function(x, name) {
  x <- .coerce_integerish(x, name = name)
  if (!length(x)) return(x)
  min_x <- min(x)
  if (min_x == 1L) {
    x <- x - 1L
  } else if (min_x != 0L) {
    cli_abort("`{name}` must use 0-based or 1-based indexing.")
  }
  x
}

# Normalize A_spatial_index to 1-based and validate range.
.normalize_dl_index <- function(x, n_max, name) {
  x <- .to_zero_based(x, name) + 1L
  if (!length(x)) return(x)
  if (any(x < 1L | x > n_max)) {
    cli_abort("`{name}` contains indices outside valid range.")
  }
  x
}

# Normalize year_i to 0-based and compute/validate n_t.
.normalize_dl_year_i <- function(year_i, n_t = NULL) {
  year_i <- .to_zero_based(year_i, name = "year_i")
  if (!length(year_i)) {
    if (is.null(n_t)) {
      cli_abort("`year_i` cannot be empty.")
    }
    if (!is.numeric(n_t) || length(n_t) != 1L || !is.finite(n_t) ||
      abs(n_t - round(n_t)) > sqrt(.Machine$double.eps) || n_t <= 0) {
      cli_abort("`n_t` must be a single positive integer.")
    }
    return(list(year_i = integer(0), n_t = as.integer(round(n_t))))
  }
  if (is.null(n_t)) {
    n_t <- max(year_i) + 1L
  } else {
    if (!is.numeric(n_t) || length(n_t) != 1L || !is.finite(n_t) ||
      abs(n_t - round(n_t)) > sqrt(.Machine$double.eps)) {
      cli_abort("`n_t` must be a single positive integer.")
    }
    n_t <- as.integer(round(n_t))
    if (n_t <= 0L) cli_abort("`n_t` must be > 0.")
    if (max(year_i) >= n_t) {
      cli_abort("`year_i` contains a time index >= `n_t`.")
    }
  }
  list(year_i = year_i, n_t = n_t)
}

.build_vertex_time_covariates <- function(covariate_data,
                                          covariates,
                                          A_st,
                                          year_i,
                                          A_spatial_index = NULL,
                                          n_t = NULL) {
  if (!inherits(covariate_data, "data.frame")) {
    cli_abort("`covariate_data` must be a data frame.")
  }
  if (!inherits(A_st, "sparseMatrix")) {
    A_st <- Matrix::Matrix(A_st, sparse = TRUE)
  }
  A_st <- methods::as(A_st, "dgCMatrix")

  n_obs <- nrow(covariate_data)
  if (length(year_i) != n_obs) {
    cli_abort("`year_i` length must equal `nrow(covariate_data)`.")
  }
  if (is.null(A_spatial_index)) {
    A_spatial_index <- seq_len(n_obs)
  }
  if (length(A_spatial_index) != n_obs) {
    cli_abort("`A_spatial_index` length must equal `nrow(covariate_data)`.")
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

  year_lu <- .normalize_dl_year_i(year_i, n_t = n_t)
  year_i <- year_lu$year_i
  n_t <- year_lu$n_t

  A_spatial_index <- .normalize_dl_index(
    A_spatial_index,
    n_max = nrow(A_st),
    name = "A_spatial_index"
  )

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

.build_distributed_lag_tmb_data <- function(distributed_lags,
                                            data,
                                            A_st,
                                            A_spatial_index,
                                            year_i,
                                            n_t,
                                            covariate_vertex_time = NULL) {
  if (is.null(distributed_lags)) {
    return(NULL)
  }

  n_vertices <- ncol(A_st)

  if (is.null(covariate_vertex_time)) {
    vertex_cov <- .build_vertex_time_covariates(
      covariate_data = data,
      covariates = distributed_lags$covariates,
      A_st = A_st,
      year_i = year_i,
      A_spatial_index = A_spatial_index,
      n_t = n_t
    )
    covariate_vertex_time <- vertex_cov$covariate_vertex_time
    n_vertices <- vertex_cov$n_vertices
    n_t <- vertex_cov$n_t
  } else {
    if (length(distributed_lags$covariates) != 1L) {
      cli_abort("`experimental$distributed_lag_covariate_vertex` currently supports exactly one distributed-lag covariate.")
    }
    if (!is.numeric(covariate_vertex_time) || any(!is.finite(covariate_vertex_time))) {
      cli_abort("`experimental$distributed_lag_covariate_vertex` must be a finite numeric vector.")
    }
    if (length(covariate_vertex_time) != n_vertices) {
      cli_abort("`experimental$distributed_lag_covariate_vertex` must have length equal to the number of mesh vertices.")
    }
    covariate_vertex_time <- array(
      rep(covariate_vertex_time, n_t),
      dim = c(n_vertices, n_t, 1L),
      dimnames = list(NULL, NULL, distributed_lags$covariates)
    )
  }

  component_levels <- c("space", "time", "spacetime")
  component_id <- match(distributed_lags$terms$component, component_levels)
  terms_df <- distributed_lags$terms
  covariates <- distributed_lags$covariates
  covariate_has_spatial <- as.integer(vapply(covariates, function(cov_name) {
    any(terms_df$variable == cov_name & terms_df$component %in% c("space", "spacetime"))
  }, logical(1L)))
  covariate_has_temporal <- as.integer(vapply(covariates, function(cov_name) {
    any(terms_df$variable == cov_name & terms_df$component == "time")
  }, logical(1L)))
  covariate_has_spacetime <- as.integer(vapply(covariates, function(cov_name) {
    any(terms_df$variable == cov_name & terms_df$component == "spacetime")
  }, logical(1L)))

  list(
    covariate_vertex_time = covariate_vertex_time,
    covariates = covariates,
    covariate_has_spatial = covariate_has_spatial,
    covariate_has_temporal = covariate_has_temporal,
    covariate_has_spacetime = covariate_has_spacetime,
    term_component = distributed_lags$terms$component,
    term_component_id = as.integer(component_id),
    term_covariate_index = as.integer(distributed_lags$terms$covariate_id),
    term_covariate_index0 = as.integer(distributed_lags$terms$covariate_id - 1L),
    term_coef_name = distributed_lags$terms$coef_name,
    n_vertices = n_vertices,
    n_t = n_t,
    n_covariates = length(covariates),
    n_terms = nrow(distributed_lags$terms)
  )
}

.solve_distributed_lag_vertex_time <- function(component, vertex_time_input, M0, M1, kappaS, kappaT, kappaST) {
  n_vertices <- nrow(vertex_time_input)
  n_t <- ncol(vertex_time_input)
  out <- matrix(0, nrow = n_vertices, ncol = n_t)

  solve_sparse <- function(A, rhs, label) {
    tryCatch(
      as.numeric(Matrix::solve(A, rhs)),
      error = function(e) {
        cli_abort(c(
          "Distributed-lag sparse solve failed in diagnostic plotting.",
          "x" = paste0(label, ": ", conditionMessage(e))
        ))
      }
    )
  }

  if (component == "space") {
    kappaS_scale <- 1 / (kappaS^2)
    system_mat <- M0 + kappaS_scale * M1
    for (tt in seq_len(n_t)) {
      rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE])
      out[, tt] <- solve_sparse(system_mat, rhs, "spatial system (M0 + kappaS^-2 * M1)")
    }
    return(out)
  }

  if (component == "time") {
    out[, 1L] <- vertex_time_input[, 1L]
    if (n_t > 1L) {
      for (tt in 2:n_t) {
        out[, tt] <- vertex_time_input[, tt] + kappaT * out[, tt - 1L]
      }
    }
    return(out)
  }

  if (component == "spacetime") {
    kappaS_scale <- 1 / (kappaS^2)
    kappaST_scale <- kappaST * kappaS_scale
    out[, 1L] <- vertex_time_input[, 1L]
    if (n_t > 1L) {
      for (tt in 2:n_t) {
        rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE] -
          kappaST_scale * (M1 %*% out[, tt - 1L, drop = TRUE]))
        out[, tt] <- solve_sparse(M0, rhs, "spatiotemporal system (M0)")
      }
    }
    return(out)
  }

  cli_abort("Unknown distributed-lag component in solver.")
}

.project_distributed_lag_vertex_time <- function(transformed_vertex_time,
                                                 A_st,
                                                 A_spatial_index,
                                                 year_i,
                                                 n_t) {
  if (!inherits(A_st, "sparseMatrix")) {
    A_st <- Matrix::Matrix(A_st, sparse = TRUE)
  }
  A_st <- methods::as(A_st, "dgCMatrix")
  n_i <- length(A_spatial_index)
  A_spatial_index <- .normalize_dl_index(
    A_spatial_index,
    n_max = nrow(A_st),
    name = "A_spatial_index"
  )
  year_lu <- .normalize_dl_year_i(year_i, n_t = n_t)
  year_i <- year_lu$year_i
  projected_by_t <- lapply(seq_len(year_lu$n_t), function(tt) {
    as.numeric(A_st %*% transformed_vertex_time[, tt, drop = TRUE])
  })
  out <- numeric(n_i)
  for (i in seq_len(n_i)) {
    out[i] <- projected_by_t[[year_i[i] + 1L]][A_spatial_index[i]]
  }
  out
}

.distributed_lag_predict_colnames <- function(term_coef_name) {
  paste0("dl_cov_", sub("^dl_", "", term_coef_name))
}

.compute_distributed_lag_term_values <- function(distributed_lags_data,
                                                 covariate_vertex_time,
                                                 A_st,
                                                 A_spatial_index,
                                                 year_i,
                                                 n_t,
                                                 M0,
                                                 M1,
                                                 log_kappaS_dl,
                                                 log_kappaT_dl,
                                                 kappaST_dl_unscaled) {
  if (is.null(distributed_lags_data)) {
    return(NULL)
  }
  n_terms <- distributed_lags_data$n_terms
  if (is.null(n_terms) || n_terms < 1L) {
    return(NULL)
  }
  n_covariates <- distributed_lags_data$n_covariates
  if (length(log_kappaS_dl) != n_covariates ||
      length(log_kappaT_dl) != n_covariates ||
      length(kappaST_dl_unscaled) != n_covariates) {
    cli_abort("Distributed-lag parameter vectors did not match the expected number of lag covariates.")
  }

  if (length(dim(covariate_vertex_time)) != 3L) {
    cli_abort("`covariate_vertex_time` must be a 3D array: [vertices, time, covariates].")
  }
  if (dim(covariate_vertex_time)[3] != n_covariates) {
    cli_abort("`covariate_vertex_time` covariate dimension does not match distributed-lag metadata.")
  }

  term_out <- matrix(0, nrow = length(A_spatial_index), ncol = n_terms)
  colnames(term_out) <- .distributed_lag_predict_colnames(distributed_lags_data$term_coef_name)

  for (term_i in seq_len(n_terms)) {
    component <- distributed_lags_data$term_component[[term_i]]
    cov_i <- distributed_lags_data$term_covariate_index[[term_i]]
    cov_slice <- matrix(
      covariate_vertex_time[, , cov_i],
      nrow = dim(covariate_vertex_time)[1],
      ncol = dim(covariate_vertex_time)[2]
    )
    kappaS <- exp(log_kappaS_dl[[cov_i]])
    kappaT <- exp(log_kappaT_dl[[cov_i]])
    kappaST <- -stats::plogis(kappaST_dl_unscaled[[cov_i]])
    transformed_vertex_time <- .solve_distributed_lag_vertex_time(
      component = component,
      vertex_time_input = cov_slice,
      M0 = M0,
      M1 = M1,
      kappaS = kappaS,
      kappaT = kappaT,
      kappaST = kappaST
    )
    term_out[, term_i] <- .project_distributed_lag_vertex_time(
      transformed_vertex_time = transformed_vertex_time,
      A_st = A_st,
      A_spatial_index = A_spatial_index,
      year_i = year_i,
      n_t = n_t
    )
  }

  term_out
}

#' Plot Distributed-Lag Diffusion from a Point on the Mesh
#'
#' Diagnostic visualization of an impulse covariate diffusing through one
#' distributed-lag component (`space()`, `time()`, or `spacetime()`).
#' Values are plotted as colored mesh triangles, so no prediction grid is
#' required.
#'
#' @param object A fitted [sdmTMB()] model with `distributed_lags`.
#' @param covariate Optional covariate name from `distributed_lags`.
#'   Required when multiple lag covariates were fitted.
#' @param component Distributed-lag component name. Must be exactly one of
#'   `"space"`, `"time"`, or `"spacetime"`.
#' @param time_value Optional time slice for the impulse. Supply either a
#'   modeled time value or a 1-based time index. Required only when the fitted
#'   distributed-lag specification includes `time()` or `spacetime()` terms with
#'   a model time column.
#' @param n_steps Number of transformed slices to plot starting at
#'   `time_value`. Defaults to 3.
#' @param common_scale Should transformed panels share a common color scale?
#'   Defaults to `FALSE`.
#' @param plot Logical. If `TRUE` (default), draw the plot.
#'
#' @return Invisibly returns a list with impulse/transformed fields on vertices,
#'   triangle summaries used for plotting, selected indices, and a `ggplot`
#'   object (when \pkg{ggplot2} is available).
#' @export
plot_distributed_lag_diffusion <- function(object,
                                           covariate = NULL,
                                           component,
                                           time_value = NULL,
                                           n_steps = 3L,
                                           common_scale = FALSE,
                                           plot = TRUE) {
  stopifnot(inherits(object, "sdmTMB"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli_abort("`ggplot2` must be installed to use `plot_distributed_lag_diffusion()`.")
  }
  if (missing(component)) {
    cli_abort("`component` is required and must be one of `space`, `time`, or `spacetime`.")
  }
  if (!component %in% c("space", "time", "spacetime")) {
    cli_abort("`component` must be exactly one of `space`, `time`, or `spacetime`.")
  }
  if (is.null(object$distributed_lags_data)) {
    cli_abort("`object` does not contain `distributed_lags_data`.")
  }
  if (!is.numeric(n_steps) || length(n_steps) != 1L || !is.finite(n_steps) || n_steps < 1L) {
    cli_abort("`n_steps` must be a single positive integer.")
  }
  n_steps <- as.integer(round(n_steps))
  if (!is.logical(plot) || length(plot) != 1L) {
    cli_abort("`plot` must be `TRUE` or `FALSE`.")
  }
  if (is.null(object$distributed_lags_parsed) ||
      is.null(object$distributed_lags_parsed$terms) ||
      nrow(object$distributed_lags_parsed$terms) == 0L) {
    cli_abort("`object` does not contain distributed-lag terms.")
  }
  terms_df <- object$distributed_lags_parsed$terms
  has_temporal_distributed_lag <- any(terms_df$component %in% c("time", "spacetime"))
  model_has_time <- !is.null(object$time) && !identical(object$time, "_sdmTMB_time")
  if (is.null(time_value) && model_has_time && has_temporal_distributed_lag) {
    cli_abort("`time_value` is required for models with temporal distributed-lag terms.")
  }
  covariates <- unique(terms_df$variable)
  if (is.null(covariate)) {
    if (length(covariates) != 1L) {
      cli_abort(c(
        "Multiple distributed-lag covariates are present.",
        "x" = "Set `covariate` to one of: {.code {paste(covariates, collapse = ', ')}}."
      ))
    }
    covariate <- covariates[[1]]
  }
  covariate <- as.character(covariate[[1]])
  if (!covariate %in% covariates) {
    cli_abort(c(
      "Unknown distributed-lag covariate.",
      "x" = "Could not find `{covariate}` in the fitted distributed-lag terms."
    ))
  }
  components_for_covariate <- unique(terms_df$component[terms_df$variable == covariate])
  if (!component %in% components_for_covariate) {
    cli_abort(c(
      "Requested component/covariate term was not fitted.",
      "x" = "No term `{component}({covariate})` in `object$distributed_lags`."
    ))
  }

  mesh <- object$spde$mesh
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
  if (component == "space") {
    time_idx <- time_i
  } else {
    time_idx <- seq.int(time_i, min(n_t, time_i + n_steps - 1L))
    if (length(time_idx) < n_steps) {
      cli_inform("Requested `n_steps` exceeded modeled time range; using available trailing slices.")
    }
  }

  covariates <- object$distributed_lags_data$covariates
  cov_i <- match(covariate, covariates)
  if (is.na(cov_i)) {
    cli_abort("Internal mismatch: selected covariate not found in `distributed_lags_data$covariates`.")
  }

  if (!is.null(object$model) &&
      !is.null(object$model$par) &&
      !is.null(object$tmb_obj) &&
      !is.null(object$tmb_obj$env)) {
    params <- object$tmb_obj$env$parList(object$model$par)
  } else if (!is.null(object$tmb_params)) {
    params <- object$tmb_params
  } else {
    cli_abort("Could not extract distributed-lag parameters from `object`.")
  }
  kappaS <- exp(params$log_kappaS_dl[cov_i])
  kappaT <- exp(params$log_kappaT_dl[cov_i])
  kappaST <- -stats::plogis(params$kappaST_dl_unscaled[cov_i])

  has_space <- as.logical(object$distributed_lags_data$covariate_has_spatial[cov_i])
  has_time <- as.logical(object$distributed_lags_data$covariate_has_temporal[cov_i])
  has_spacetime <- as.logical(object$distributed_lags_data$covariate_has_spacetime[cov_i])
  if (component == "space" && !has_space) {
    cli_abort("Selected covariate was not fitted with a spatial distributed-lag component.")
  }
  if (component == "time" && !has_time) {
    cli_abort("Selected covariate was not fitted with a temporal distributed-lag component.")
  }
  if (component == "spacetime" && !has_spacetime) {
    cli_abort("Selected covariate was not fitted with a spatiotemporal distributed-lag component.")
  }

  n_vertices <- nrow(loc)
  impulse_vertex_time <- matrix(0, nrow = n_vertices, ncol = n_t)
  impulse_vertex_time[vertex_i, time_i] <- 1

  M0 <- object$tmb_data$spde$M0
  M1 <- object$tmb_data$spde$M1
  transformed_vertex_time <- .solve_distributed_lag_vertex_time(
    component = component,
    vertex_time_input = impulse_vertex_time,
    M0 = M0,
    M1 = M1,
    kappaS = kappaS,
    kappaT = kappaT,
    kappaST = kappaST
  )

  panel_fields <- vector("list", length(time_idx) + 1L)
  panel_titles <- character(length(panel_fields))
  panel_fields[[1L]] <- impulse_vertex_time[, time_i]
  panel_titles[[1L]] <- paste0("original (t=", time_values[time_i], ")")
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

  n_tri <- nrow(tv)
  triangle_values <- vapply(panel_fields, function(v) {
    rowMeans(matrix(v[tv], nrow = n_tri, ncol = 3L))
  }, numeric(n_tri))
  if (!is.matrix(triangle_values)) {
    triangle_values <- matrix(triangle_values, ncol = 1L)
  }
  colnames(triangle_values) <- panel_titles

  tri_x <- matrix(loc[tv, 1], nrow = n_tri, ncol = 3L)
  tri_y <- matrix(loc[tv, 2], nrow = n_tri, ncol = 3L)
  xlim <- range(loc[, 1])
  ylim <- range(loc[, 2])

  triangle_df <- do.call(rbind, lapply(seq_len(ncol(triangle_values)), function(j) {
    data.frame(
      panel = panel_titles[j],
      tri = rep(seq_len(n_tri), each = 3L),
      x = as.vector(t(tri_x)),
      y = as.vector(t(tri_y)),
      value = rep(triangle_values[, j], each = 3L),
      stringsAsFactors = FALSE
    )
  }))
  triangle_df$panel <- factor(triangle_df$panel, levels = panel_titles)
  point_df <- data.frame(
    panel = factor(panel_titles, levels = panel_titles),
    x = rep(loc[vertex_i, 1], length(panel_titles)),
    y = rep(loc[vertex_i, 2], length(panel_titles))
  )

  triangle_df$value_plot <- triangle_df$value
  fill_name <- "Value"
  if (!isTRUE(common_scale)) {
    fill_name <- "Relative value"
    for (p in levels(triangle_df$panel)) {
      i <- which(triangle_df$panel == p)
      rng <- range(triangle_df$value[i], finite = TRUE)
      if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        triangle_df$value_plot[i] <- 0
      } else {
        triangle_df$value_plot[i] <- (triangle_df$value[i] - rng[1]) / (rng[2] - rng[1])
      }
    }
  }

  fill_limits <- if (isTRUE(common_scale)) NULL else c(0, 1)
  xlab <- if (!is.null(object$spde$xy_cols) && length(object$spde$xy_cols) >= 2L) object$spde$xy_cols[1] else "x"
  ylab <- if (!is.null(object$spde$xy_cols) && length(object$spde$xy_cols) >= 2L) object$spde$xy_cols[2] else "y"
  plot_obj <- ggplot2::ggplot(
    triangle_df,
    ggplot2::aes(x = .data$x, y = .data$y, group = interaction(.data$panel, .data$tri))
  ) +
    ggplot2::geom_polygon(
      ggplot2::aes(fill = .data$value_plot),
      color = "#FFFFFF10", linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      shape = 21,
      fill = "black",
      color = "black",
      size = 1.4
    ) +
    ggplot2::facet_wrap(stats::as.formula("~ panel"), nrow = 1L) +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::scale_fill_viridis_c(
      limits = fill_limits,
      name = fill_name, option = "C"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab)

  if (plot) {
    print(plot_obj)
  }

  invisible(list(
    covariate = covariate,
    component = component,
    vertex = vertex_i,
    vertex_xy = loc[vertex_i, , drop = TRUE],
    impulse_time_index = time_i,
    impulse_time_value = time_values[time_i],
    transformed_time_index = time_idx,
    transformed_time_values = time_values[time_idx],
    impulse_vertex_time = impulse_vertex_time,
    transformed_vertex_time = transformed_vertex_time,
    triangle_values = triangle_values,
    triangle_df = triangle_df,
    plot = plot_obj,
    mesh_loc = loc,
    mesh_triangles = tv
  ))
}
