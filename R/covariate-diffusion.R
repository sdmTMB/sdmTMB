.as_dgC <- function(x) {
  if (!inherits(x, "sparseMatrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  methods::as(x, "dgCMatrix")
}

.extract_covariate_diffusion_term_exprs <- function(expr) {
  if (is.call(expr)) {
    fn <- as.character(expr[[1]])
    if (identical(fn, "+")) {
      if (length(expr) == 2L) {
        return(.extract_covariate_diffusion_term_exprs(expr[[2]]))
      }
      if (length(expr) == 3L) {
        return(c(
          .extract_covariate_diffusion_term_exprs(expr[[2]]),
          .extract_covariate_diffusion_term_exprs(expr[[3]])
        ))
      }
    }
    if (identical(fn, "(") && length(expr) == 2L) {
      return(.extract_covariate_diffusion_term_exprs(expr[[2]]))
    }
  }
  list(expr)
}

.append_covariate_diffusion_coef_columns <- function(X, coef_names) {
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

.parse_covariate_diffusion_formula <- function(covariate_diffusion) {
  if (is.null(covariate_diffusion)) {
    return(NULL)
  }

  if (!inherits(covariate_diffusion, "formula")) {
    cli_abort("`covariate_diffusion` must be `NULL` or a one-sided formula.")
  }

  if (length(covariate_diffusion) != 2L) {
    cli_abort("`covariate_diffusion` must be a one-sided formula such as `~ space(x) + time(x)`.")
  }

  term_exprs <- .extract_covariate_diffusion_term_exprs(covariate_diffusion[[2]])
  if (!length(term_exprs)) {
    cli_abort("`covariate_diffusion` must include at least one lag term.")
  }

  allowed_wrappers <- c("space", "time", "spacetime")

  parsed_terms <- lapply(term_exprs, function(expr) {
    term_label <- paste(deparse(expr), collapse = "")

    if (!is.call(expr)) {
      cli_abort(c(
        "Unsupported term in `covariate_diffusion`.",
        "i" = "Terms must be wrapped in `space()` or `time()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    wrapper <- as.character(expr[[1]])
    if (!wrapper %in% allowed_wrappers) {
      cli_abort(c(
        "Unsupported wrapper in `covariate_diffusion`.",
        "i" = "Allowed wrappers are `space()` and `time()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    if (length(expr) != 2L || !is.symbol(expr[[2]])) {
      cli_abort(c(
        "Unsupported `covariate_diffusion` term structure.",
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
      "Duplicate `covariate_diffusion` terms are not supported.",
      "x" = "Duplicated term(s): {.code {paste(dup_labels, collapse = ', ')}}"
    ))
  }

  unique_covariates <- unique(terms_df$variable)
  terms_df$covariate_id <- match(terms_df$variable, unique_covariates)
  terms_df$coef_name <- paste0("dl_", terms_df$component, "_", make.names(terms_df$variable))

  list(
    formula = covariate_diffusion,
    terms = terms_df,
    covariates = unique_covariates,
    needs_time = any(terms_df$component %in% c("time", "spacetime"))
  )
}

.validate_covariate_diffusion_terms <- function(covariate_diffusion, data, time, multi_family) {
  if (is.null(covariate_diffusion)) {
    return(NULL)
  }

  if (isTRUE(multi_family)) {
    cli_abort("`covariate_diffusion` is currently unsupported for multi-family models.")
  }

  if (isTRUE(covariate_diffusion$needs_time) && is.null(time)) {
    cli_abort(
      "Temporal `covariate_diffusion` terms require a `time` argument."
    )
  }

  missing_covariates <- setdiff(covariate_diffusion$covariates, names(data))
  if (length(missing_covariates)) {
    cli_abort(c(
      "Missing covariate diffusion covariate(s) in `data`.",
      "x" = "Missing: {.code {paste(missing_covariates, collapse = ', ')}}"
    ))
  }

  non_numeric <- covariate_diffusion$covariates[!vapply(covariate_diffusion$covariates, function(v) {
    is.numeric(data[[v]])
  }, logical(1L))]

  if (length(non_numeric)) {
    cli_abort(c(
      "Distributed lag covariates must be numeric.",
      "x" = "Non-numeric covariate(s): {.code {paste(non_numeric, collapse = ', ')}}"
    ))
  }

  covariate_diffusion
}

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

.normalize_dl_index <- function(x, n_max, name) {
  x <- .to_zero_based(x, name) + 1L
  if (!length(x)) return(x)
  if (any(x < 1L | x > n_max)) {
    cli_abort("`{name}` contains indices outside valid range.")
  }
  x
}

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
  A_st <- .as_dgC(A_st)

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

.build_covariate_diffusion_tmb_data <- function(covariate_diffusion,
                                            data,
                                            A_st,
                                            A_spatial_index,
                                            year_i,
                                            n_t,
                                            covariate_vertex_time = NULL) {
  if (is.null(covariate_diffusion)) {
    return(NULL)
  }

  n_vertices <- ncol(A_st)

  if (is.null(covariate_vertex_time)) {
    vertex_cov <- .build_vertex_time_covariates(
      covariate_data = data,
      covariates = covariate_diffusion$covariates,
      A_st = A_st,
      year_i = year_i,
      A_spatial_index = A_spatial_index,
      n_t = n_t
    )
    covariate_vertex_time <- vertex_cov$covariate_vertex_time
    n_vertices <- vertex_cov$n_vertices
    n_t <- vertex_cov$n_t
  } else {
    if (!is.numeric(covariate_vertex_time) || any(!is.finite(covariate_vertex_time))) {
      cli_abort("`experimental$covariate_diffusion_covariate_vertex` must be finite numeric values.")
    }
    if (is.null(dim(covariate_vertex_time))) {
      if (length(covariate_diffusion$covariates) != 1L) {
        cli_abort("A vector `experimental$covariate_diffusion_covariate_vertex` supports exactly one covariate-diffusion covariate.")
      }
      if (length(covariate_vertex_time) != n_vertices) {
        cli_abort("A vector `experimental$covariate_diffusion_covariate_vertex` must have length equal to the number of mesh vertices.")
      }
      covariate_vertex_time <- array(
        rep(covariate_vertex_time, n_t),
        dim = c(n_vertices, n_t, 1L),
        dimnames = list(NULL, NULL, covariate_diffusion$covariates)
      )
    } else if (length(dim(covariate_vertex_time)) == 2L) {
      if (length(covariate_diffusion$covariates) != 1L) {
        cli_abort("A matrix `experimental$covariate_diffusion_covariate_vertex` supports exactly one covariate-diffusion covariate.")
      }
      if (!all(dim(covariate_vertex_time) == c(n_vertices, n_t))) {
        cli_abort("A matrix `experimental$covariate_diffusion_covariate_vertex` must have dimensions `n_vertices` by `n_t`.")
      }
      covariate_vertex_time <- array(
        covariate_vertex_time,
        dim = c(n_vertices, n_t, 1L),
        dimnames = list(NULL, NULL, covariate_diffusion$covariates)
      )
    } else if (length(dim(covariate_vertex_time)) == 3L) {
      expected_dim <- c(n_vertices, n_t, length(covariate_diffusion$covariates))
      if (!all(dim(covariate_vertex_time) == expected_dim)) {
        cli_abort("An array `experimental$covariate_diffusion_covariate_vertex` must have dimensions `n_vertices` by `n_t` by `n_covariates`.")
      }
      dimnames(covariate_vertex_time) <- list(NULL, NULL, covariate_diffusion$covariates)
    } else {
      cli_abort("`experimental$covariate_diffusion_covariate_vertex` must be a vector, matrix, or 3D array.")
    }
  }

  component_levels <- c("space", "time", "spacetime")
  component_id <- match(covariate_diffusion$terms$component, component_levels)
  terms_df <- covariate_diffusion$terms
  covariates <- covariate_diffusion$covariates
  covariate_has_spatial <- integer(length(covariates))
  covariate_has_temporal <- integer(length(covariates))
  covariate_has_spacetime <- integer(length(covariates))
  for (i in seq_along(covariates)) {
    components <- terms_df$component[terms_df$variable == covariates[i]]
    covariate_has_spatial[i]   <- any(components %in% c("space", "spacetime"))
    covariate_has_temporal[i]  <- any(components == "time")
    covariate_has_spacetime[i] <- any(components == "spacetime")
  }

  list(
    covariate_vertex_time = covariate_vertex_time,
    covariates = covariates,
    covariate_has_spatial = covariate_has_spatial,
    covariate_has_temporal = covariate_has_temporal,
    covariate_has_spacetime = covariate_has_spacetime,
    term_component = covariate_diffusion$terms$component,
    term_component_id = as.integer(component_id),
    term_covariate_index = as.integer(covariate_diffusion$terms$covariate_id),
    term_covariate_index0 = as.integer(covariate_diffusion$terms$covariate_id - 1L),
    term_coef_name = covariate_diffusion$terms$coef_name,
    n_vertices = n_vertices,
    n_t = n_t,
    n_covariates = length(covariates),
    n_terms = nrow(covariate_diffusion$terms)
  )
}

.solve_covariate_diffusion_vertex_time <- function(component, vertex_time_input, M0, M1,
                                                  kappaS, kappaT, kappaST,
                                                  has_space = NULL,
                                                  has_time = NULL,
                                                  has_spacetime = NULL) {
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
    has_spacetime <- isTRUE(has_spacetime)
    kappaS_scale <- if (has_space || has_spacetime) 1 / (kappaS^2) else 0
    kappaT_scale <- if (has_time) kappaT else 0
    kappaST_scale <- if (has_spacetime) kappaST * kappaS_scale else 0
    system_mat <- M0 + kappaS_scale * M1
    for (tt in seq_len(n_t)) {
      rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE])
      if (tt > 1L) {
        if (kappaT_scale != 0) {
          rhs <- rhs + as.numeric(kappaT_scale * (M0 %*% out[, tt - 1L, drop = TRUE]))
        }
        if (kappaST_scale != 0) {
          rhs <- rhs - as.numeric(kappaST_scale * (M1 %*% out[, tt - 1L, drop = TRUE]))
        }
      }
      out[, tt] <- solve_sparse(system_mat, rhs, "combined system (space + time + spacetime)")
    }
    return(out)
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
    denom <- 1 + kappaT
    out[, 1L] <- vertex_time_input[, 1L] / denom
    if (n_t > 1L) {
      for (tt in 2:n_t) {
        out[, tt] <- (vertex_time_input[, tt] + kappaT * out[, tt - 1L]) / denom
      }
    }
    return(out)
  }

  if (component == "spacetime") {
    kappaS_scale <- 1 / (kappaS^2)
    kappaST_scale <- kappaST * kappaS_scale
    system_mat <- M0 + (kappaS_scale - kappaST_scale) * M1
    for (tt in seq_len(n_t)) {
      rhs <- as.numeric(M0 %*% vertex_time_input[, tt, drop = TRUE])
      if (tt > 1L) {
        rhs <- rhs - as.numeric(kappaST_scale * (M1 %*% out[, tt - 1L, drop = TRUE]))
      }
      out[, tt] <- solve_sparse(system_mat, rhs, "spatiotemporal system (M0 + (kappaS^-2 - kappaST * kappaS^-2) * M1)")
    }
    return(out)
  }

  cli_abort("Unknown covariate-diffusion component in solver.")
}

.project_covariate_diffusion_vertex_time <- function(transformed_vertex_time,
                                                 A_st,
                                                 A_spatial_index,
                                                 year_i,
                                                 n_t) {
  A_st <- .as_dgC(A_st)
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

.covariate_diffusion_predict_colnames <- function(term_coef_name) {
  paste0("diffusion_cov_", sub("^dl_", "", term_coef_name))
}

.append_covariate_diffusion_term_values <- function(nd, object, tmb_obj, lp,
                                                covariate_vertex_time, A_st,
                                                A_spatial_index, year_i, n_t) {
  par_list <- tmb_obj$env$parList(lp)
  dl_term_values <- .compute_covariate_diffusion_term_values(
    covariate_diffusion_data = object$covariate_diffusion_data,
    covariate_vertex_time = covariate_vertex_time,
    A_st = A_st,
    A_spatial_index = A_spatial_index,
    year_i = year_i,
    n_t = n_t,
    M0 = object$tmb_data$spde$M0,
    M1 = object$tmb_data$spde$M1,
    log_kappaS_dl = par_list$log_kappaS_dl,
    kappaT_dl_raw = par_list$kappaT_dl_raw,
    kappaST_dl_raw = par_list$kappaST_dl_raw
  )
  cbind(nd, as.data.frame(dl_term_values))
}

.compute_covariate_diffusion_term_values <- function(covariate_diffusion_data,
                                                 covariate_vertex_time,
                                                 A_st,
                                                 A_spatial_index,
                                                 year_i,
                                                 n_t,
                                                 M0,
                                                 M1,
                                                 log_kappaS_dl,
                                                 kappaT_dl_raw,
                                                 kappaST_dl_raw) {
  if (is.null(covariate_diffusion_data)) {
    return(NULL)
  }
  n_terms <- covariate_diffusion_data$n_terms
  if (is.null(n_terms) || n_terms < 1L) {
    return(NULL)
  }
  n_covariates <- covariate_diffusion_data$n_covariates
  if (length(log_kappaS_dl) != n_covariates ||
      length(kappaT_dl_raw) != n_covariates ||
      length(kappaST_dl_raw) != n_covariates) {
    cli_abort("Covariate diffusion parameter vectors did not match the expected number of lag covariates.")
  }

  if (length(dim(covariate_vertex_time)) != 3L) {
    cli_abort("`covariate_vertex_time` must be a 3D array: [vertices, time, covariates].")
  }
  if (dim(covariate_vertex_time)[3] != n_covariates) {
    cli_abort("`covariate_vertex_time` covariate dimension does not match covariate-diffusion metadata.")
  }

  term_out <- matrix(0, nrow = length(A_spatial_index), ncol = n_terms)
  colnames(term_out) <- .covariate_diffusion_predict_colnames(covariate_diffusion_data$term_coef_name)

  for (term_i in seq_len(n_terms)) {
    component <- covariate_diffusion_data$term_component[[term_i]]
    cov_i <- covariate_diffusion_data$term_covariate_index[[term_i]]
    cov_slice <- matrix(
      covariate_vertex_time[, , cov_i],
      nrow = dim(covariate_vertex_time)[1],
      ncol = dim(covariate_vertex_time)[2]
    )
    kappaS <- exp(log_kappaS_dl[[cov_i]])
    kappaT <- kappaT_dl_raw[[cov_i]]
    kappaST <- kappaST_dl_raw[[cov_i]]
    transformed_vertex_time <- .solve_covariate_diffusion_vertex_time(
      component = component,
      vertex_time_input = cov_slice,
      M0 = M0,
      M1 = M1,
      kappaS = kappaS,
      kappaT = kappaT,
      kappaST = kappaST
    )
    term_out[, term_i] <- .project_covariate_diffusion_vertex_time(
      transformed_vertex_time = transformed_vertex_time,
      A_st = A_st,
      A_spatial_index = A_spatial_index,
      year_i = year_i,
      n_t = n_t
    )
  }

  term_out
}

.dl_plot_extract_mesh <- function(mesh) {
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

.dl_plot_resolve_time <- function(object, component, time_value, n_steps,
                                  allow_spatial_steps = FALSE) {
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
  if (component == "space" && !isTRUE(allow_spatial_steps)) {
    time_idx <- time_i
  } else {
    time_idx <- seq.int(time_i, min(n_t, time_i + n_steps - 1L))
    if (length(time_idx) < n_steps) {
      cli_inform("Requested `n_steps` exceeded modeled time range; using available trailing slices.")
    }
  }
  list(n_t = n_t, time_values = time_values, time_i = time_i, time_idx = time_idx)
}

.dl_plot_extract_kappas <- function(object, cov_i) {
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
    kappaS = exp(params$log_kappaS_dl[cov_i]),
    kappaT = params$kappaT_dl_raw[cov_i],
    kappaST = params$kappaST_dl_raw[cov_i]
  )
}

.dl_plot_context <- function(object, covariate, component, component_missing,
                             time_value, n_steps, allow_spatial_steps,
                             function_name) {
  stopifnot(inherits(object, "sdmTMB"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli_abort("`ggplot2` must be installed to use `{function_name}()`.")
  }
  if (isTRUE(component_missing)) {
    cli_abort("`component` is required and must be one of `space`, `time`, `spacetime`, or `combined`.")
  }
  if (!component %in% c("space", "time", "spacetime", "combined")) {
    cli_abort("`component` must be exactly one of `space`, `time`, `spacetime`, or `combined`.")
  }
  if (is.null(object$covariate_diffusion_data)) {
    cli_abort("`object` does not contain `covariate_diffusion_data`.")
  }
  if (!is.numeric(n_steps) || length(n_steps) != 1L || !is.finite(n_steps) || n_steps < 1L) {
    cli_abort("`n_steps` must be a single positive integer.")
  }
  n_steps <- as.integer(round(n_steps))
  if (is.null(object$covariate_diffusion_parsed) ||
      is.null(object$covariate_diffusion_parsed$terms) ||
      nrow(object$covariate_diffusion_parsed$terms) == 0L) {
    cli_abort("`object` does not contain covariate-diffusion terms.")
  }

  terms_df <- object$covariate_diffusion_parsed$terms
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
      "x" = "No term `{component}({covariate})` in `object$covariate_diffusion`."
    ))
  }
  component_for_time <- if (component == "combined" && !any(components_for_covariate %in% c("time", "spacetime"))) "space" else component

  mesh_info <- .dl_plot_extract_mesh(object$spde$mesh)
  time_info <- .dl_plot_resolve_time(
    object, component_for_time, time_value, n_steps,
    allow_spatial_steps = allow_spatial_steps
  )

  cov_i <- match(covariate, object$covariate_diffusion_data$covariates)
  if (is.na(cov_i)) {
    cli_abort("Internal mismatch: selected covariate not found in `covariate_diffusion_data$covariates`.")
  }

  xy_cols <- object$spde$xy_cols
  list(
    covariate = covariate,
    component = component,
    mesh_info = mesh_info,
    time_info = time_info,
    cov_i = cov_i,
    kappas = .dl_plot_extract_kappas(object, cov_i),
    has_space = as.logical(object$covariate_diffusion_data$covariate_has_spatial[cov_i]),
    has_time = as.logical(object$covariate_diffusion_data$covariate_has_temporal[cov_i]),
    has_spacetime = as.logical(object$covariate_diffusion_data$covariate_has_spacetime[cov_i]),
    xlab = if (!is.null(xy_cols) && length(xy_cols) >= 2L) xy_cols[1] else "x",
    ylab = if (!is.null(xy_cols) && length(xy_cols) >= 2L) xy_cols[2] else "y"
  )
}

.dl_plot_transform_values <- function(x, value_transform) {
  switch(value_transform,
    identity = x,
    sqrt = {
      if (any(x < 0, na.rm = TRUE)) {
        cli_abort("`value_transform = \"sqrt\"` requires non-negative values. Use `\"signed_sqrt\"` for signed values.")
      }
      sqrt(x)
    },
    signed_sqrt = sign(x) * sqrt(abs(x))
  )
}

.dl_plot_time_panels <- function(first_field, first_title, transformed_vertex_time,
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

.dl_plot_panel_dfs <- function(loc, tv, panel_fields, panel_titles, vertex_i,
                               common_scale, value_transform) {
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

  triangle_df$value_plot <- .dl_plot_transform_values(triangle_df$value, value_transform)
  fill_name <- if (value_transform == "identity") "Value" else paste0("Value (", value_transform, ")")
  if (!isTRUE(common_scale)) {
    fill_name <- if (value_transform == "identity") "Relative value" else paste0("Relative value (", value_transform, ")")
    for (p in levels(triangle_df$panel)) {
      i <- which(triangle_df$panel == p)
      rng <- range(triangle_df$value_plot[i], finite = TRUE)
      if (!all(is.finite(rng)) || rng[1] == rng[2]) {
        triangle_df$value_plot[i] <- 0
      } else {
        triangle_df$value_plot[i] <- (triangle_df$value_plot[i] - rng[1]) / (rng[2] - rng[1])
      }
    }
  }
  fill_limits <- if (isTRUE(common_scale)) NULL else c(0, 1)

  list(
    triangle_values = triangle_values,
    triangle_df = triangle_df,
    point_df = point_df,
    fill_name = fill_name,
    fill_limits = fill_limits
  )
}

.dl_plot_ggplot <- function(triangle_df, point_df, fill_limits, fill_name, xlim, ylim, xlab, ylab) {
  ggplot2::ggplot(
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
}

#' Plot Covariate-Diffusion Diagnostics on the Mesh
#'
#' Visualize fitted covariate-diffusion transforms or impulse-response kernels
#' for one selected covariate-diffusion term.
#' Values are plotted as colored mesh triangles, so no prediction grid is
#' required.
#'
#' @param object A fitted [sdmTMB()] model with `covariate_diffusion`.
#' @param covariate Optional covariate name from `covariate_diffusion`.
#'   Required when multiple lag covariates were fitted.
#' @param component Covariate-diffusion component name. Must be one of
#'   `"space"`, `"time"`, or `"combined"`. `"combined"`
#'   plots the joint response of all covariate-diffusion components fitted for
#'   `covariate`.
#' @param time_value Optional time slice to plot or use for the impulse. Supply
#'   either a modeled time value or a 1-based time index. Defaults to 1.
#' @param n_steps Number of transformed slices to plot starting at
#'   `time_value`. Defaults to 1 for `plot_diffused_covariate()` and 3 for
#'   `plot_diffusion_kernel()`.
#' @param common_scale Should plotted panels share a common color scale?
#'   Defaults to `TRUE` for `plot_diffused_covariate()` and `FALSE` for
#'   `plot_diffusion_kernel()`.
#' @param plot Should the plot be printed? Defaults to `TRUE`.
#'
#' @details
#' `plot_diffused_covariate()` visualizes the original mesh-vertex covariate
#' field and its fitted covariate-diffusion transform of one selected covariate
#' time slice over one or more lagged output time slices.
#'
#' `plot_diffusion_kernel()` visualizes an impulse covariate diffusing through
#' one covariate-diffusion component.
#'
#' @return Invisibly returns a list with fields on vertices, triangle summaries
#'   used for plotting, selected indices, and a `ggplot` object.
#'
#' @examplesIf ggplot2_installed()
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
#'
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
#'   covariate_diffusion = ~ space(x1) + time(x1),
#'   diffusion_kappaS = 4.4,
#'   diffusion_rhoT = 0.3,
#'   seed = 123
#' )
#' dat$observed <- sim$observed
#'
#' fit <- sdmTMB(
#'   observed ~ 1,
#'   data = dat,
#'   mesh = mesh,
#'   time = "year",
#'   spatial = "off",
#'   spatiotemporal = "off",
#'   family = gaussian(),
#'   covariate_diffusion = ~ space(x1) + time(x1)
#' )
#'
#' plot_diffused_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "space"
#' )
#' plot_diffused_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "time",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' plot_diffused_covariate(
#'   fit,
#'   covariate = "x1",
#'   component = "combined",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' plot_diffusion_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "space"
#' )
#' plot_diffusion_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "time",
#'   time_value = 1,
#'   n_steps = 2,
#'   common_scale = TRUE #<
#' )
#' plot_diffusion_kernel(
#'   fit,
#'   covariate = "x1",
#'   component = "combined",
#'   time_value = 1,
#'   n_steps = 2
#' )
#' @rdname covariate_diffusion_plots
#' @export
plot_diffused_covariate <- function(object,
                                    covariate = NULL,
                                    component,
                                    time_value = 1,
                                    n_steps = 1L,
                                    common_scale = TRUE,
                                    plot = TRUE) {
  ctx <- .dl_plot_context(
    object = object,
    covariate = covariate,
    component = if (missing(component)) NULL else component,
    component_missing = missing(component),
    time_value = time_value,
    n_steps = n_steps,
    function_name = "plot_diffused_covariate",
    allow_spatial_steps = FALSE
  )

  original_vertex_time <- matrix(
    object$covariate_diffusion_data$covariate_vertex_time[, , ctx$cov_i],
    nrow = object$covariate_diffusion_data$n_vertices,
    ncol = ctx$time_info$n_t
  )

  source_vertex_time <- matrix(0, nrow = nrow(original_vertex_time), ncol = ctx$time_info$n_t)
  source_vertex_time[, ctx$time_info$time_i] <- original_vertex_time[, ctx$time_info$time_i]

  transformed_vertex_time <- .solve_covariate_diffusion_vertex_time(
    component = ctx$component,
    vertex_time_input = source_vertex_time,
    M0 = object$tmb_data$spde$M0,
    M1 = object$tmb_data$spde$M1,
    kappaS = ctx$kappas$kappaS,
    kappaT = ctx$kappas$kappaT,
    kappaST = ctx$kappas$kappaST,
    has_space = ctx$has_space,
    has_time = ctx$has_time,
    has_spacetime = ctx$has_spacetime
  )

  time_values <- ctx$time_info$time_values
  time_i <- ctx$time_info$time_i
  time_idx <- ctx$time_info$time_idx
  panels <- .dl_plot_time_panels(
    first_field = original_vertex_time[, time_i],
    first_title = paste0("original (t=", time_values[time_i], ")"),
    transformed_vertex_time = transformed_vertex_time,
    time_i = time_i,
    time_idx = time_idx,
    time_values = time_values
  )
  panel <- .dl_plot_panel_dfs(
    loc = ctx$mesh_info$loc,
    tv = ctx$mesh_info$tv,
    panel_fields = panels$fields,
    panel_titles = panels$titles,
    vertex_i = ctx$mesh_info$vertex_i,
    common_scale = common_scale,
    value_transform = "identity"
  )
  point_df <- data.frame(
    panel = factor(character(0L), levels = panels$titles),
    x = numeric(0L),
    y = numeric(0L)
  )
  plot_obj <- .dl_plot_ggplot(
    triangle_df = panel$triangle_df,
    point_df = point_df,
    fill_limits = panel$fill_limits,
    fill_name = panel$fill_name,
    xlim = range(ctx$mesh_info$loc[, 1]),
    ylim = range(ctx$mesh_info$loc[, 2]),
    xlab = ctx$xlab,
    ylab = ctx$ylab
  )

  if (isTRUE(plot)) {
    print(plot_obj)
  }

  invisible(list(
    covariate = ctx$covariate,
    component = ctx$component,
    time_index = time_i,
    time_value = time_values[time_i],
    transformed_time_index = time_idx,
    transformed_time_values = time_values[time_idx],
    original_vertex_time = original_vertex_time,
    source_vertex_time = source_vertex_time,
    transformed_vertex_time = transformed_vertex_time,
    triangle_values = panel$triangle_values,
    triangle_df = panel$triangle_df,
    plot = plot_obj,
    mesh_loc = ctx$mesh_info$loc,
    mesh_triangles = ctx$mesh_info$tv
  ))
}

#' @rdname covariate_diffusion_plots
#' @export
plot_diffusion_kernel <- function(object,
                                           covariate = NULL,
                                           component,
                                           time_value = 1,
                                           n_steps = 3L,
                                           common_scale = FALSE,
                                           plot = TRUE) {
  ctx <- .dl_plot_context(
    object = object,
    covariate = covariate,
    component = if (missing(component)) NULL else component,
    component_missing = missing(component),
    time_value = time_value,
    n_steps = n_steps,
    function_name = "plot_diffusion_kernel",
    allow_spatial_steps = FALSE
  )

  time_values <- ctx$time_info$time_values
  n_t <- ctx$time_info$n_t
  time_i <- ctx$time_info$time_i
  time_idx <- ctx$time_info$time_idx
  vertex_i <- ctx$mesh_info$vertex_i
  n_vertices <- nrow(ctx$mesh_info$loc)
  impulse_vertex_time <- matrix(0, nrow = n_vertices, ncol = n_t)
  impulse_vertex_time[vertex_i, time_i] <- 1

  transformed_vertex_time <- .solve_covariate_diffusion_vertex_time(
    component = ctx$component,
    vertex_time_input = impulse_vertex_time,
    M0 = object$tmb_data$spde$M0,
    M1 = object$tmb_data$spde$M1,
    kappaS = ctx$kappas$kappaS,
    kappaT = ctx$kappas$kappaT,
    kappaST = ctx$kappas$kappaST,
    has_space = ctx$has_space,
    has_time = ctx$has_time,
    has_spacetime = ctx$has_spacetime
  )

  panels <- .dl_plot_time_panels(
    first_field = impulse_vertex_time[, time_i],
    first_title = paste0("original (t=", time_values[time_i], ")"),
    transformed_vertex_time = transformed_vertex_time,
    time_i = time_i,
    time_idx = time_idx,
    time_values = time_values
  )
  panel <- .dl_plot_panel_dfs(
    ctx$mesh_info$loc, ctx$mesh_info$tv, panels$fields, panels$titles, vertex_i,
    common_scale = common_scale,
    value_transform = "identity"
  )

  plot_obj <- .dl_plot_ggplot(
    triangle_df = panel$triangle_df,
    point_df = panel$point_df,
    fill_limits = panel$fill_limits,
    fill_name = panel$fill_name,
    xlim = range(ctx$mesh_info$loc[, 1]),
    ylim = range(ctx$mesh_info$loc[, 2]),
    xlab = ctx$xlab,
    ylab = ctx$ylab
  )

  if (isTRUE(plot)) {
    print(plot_obj)
  }

  invisible(list(
    covariate = ctx$covariate,
    component = ctx$component,
    vertex = vertex_i,
    vertex_xy = ctx$mesh_info$loc[vertex_i, , drop = TRUE],
    impulse_time_index = time_i,
    impulse_time_value = time_values[time_i],
    transformed_time_index = time_idx,
    transformed_time_values = time_values[time_idx],
    impulse_vertex_time = impulse_vertex_time,
    transformed_vertex_time = transformed_vertex_time,
    triangle_values = panel$triangle_values,
    triangle_df = panel$triangle_df,
    plot = plot_obj,
    mesh_loc = ctx$mesh_info$loc,
    mesh_triangles = ctx$mesh_info$tv
  ))
}
