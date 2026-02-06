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

  if (isTRUE(delta)) {
    cli_abort("`distributed_lags` is currently unsupported for delta/hurdle models.")
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
                                            n_t) {
  if (is.null(distributed_lags)) {
    return(NULL)
  }

  vertex_cov <- .build_vertex_time_covariates(
    covariate_data = data,
    covariates = distributed_lags$covariates,
    A_st = A_st,
    year_i = year_i,
    A_spatial_index = A_spatial_index,
    n_t = n_t
  )

  component_levels <- c("space", "time", "spacetime")
  component_id <- match(distributed_lags$terms$component, component_levels)

  list(
    covariate_vertex_time = vertex_cov$covariate_vertex_time,
    covariates = distributed_lags$covariates,
    term_component = distributed_lags$terms$component,
    term_component_id = as.integer(component_id),
    term_covariate_index = as.integer(distributed_lags$terms$covariate_id),
    term_covariate_index0 = as.integer(distributed_lags$terms$covariate_id - 1L),
    term_coef_name = distributed_lags$terms$coef_name,
    n_vertices = vertex_cov$n_vertices,
    n_t = vertex_cov$n_t,
    n_covariates = length(distributed_lags$covariates),
    n_terms = nrow(distributed_lags$terms)
  )
}
