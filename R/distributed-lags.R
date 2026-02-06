.parse_distributed_lags_formula <- function(distributed_lags) {
  if (is.null(distributed_lags)) {
    return(NULL)
  }

  if (!inherits(distributed_lags, "formula")) {
    cli_abort("`distributed_lags` must be `NULL` or a one-sided formula.")
  }

  if (length(distributed_lags) != 2L) {
    cli_abort("`distributed_lags` must be a one-sided formula such as `~ spatial(x) + temporal(x)`.")
  }

  term_labels <- attr(stats::terms(distributed_lags), "term.labels")
  if (!length(term_labels)) {
    cli_abort("`distributed_lags` must include at least one lag term.")
  }

  allowed_wrappers <- c("spatial", "temporal", "spatiotemporal")

  parsed_terms <- lapply(term_labels, function(term_label) {
    expr <- str2lang(term_label)

    if (!is.call(expr)) {
      cli_abort(c(
        "Unsupported term in `distributed_lags`.",
        "i" = "Terms must be wrapped in `spatial()`, `temporal()`, or `spatiotemporal()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    wrapper <- as.character(expr[[1]])
    if (!wrapper %in% allowed_wrappers) {
      cli_abort(c(
        "Unsupported wrapper in `distributed_lags`.",
        "i" = "Allowed wrappers are `spatial()`, `temporal()`, and `spatiotemporal()`.",
        "x" = "Problematic term: {.code {term_label}}"
      ))
    }

    if (length(expr) != 2L || !is.symbol(expr[[2]])) {
      cli_abort(c(
        "Unsupported `distributed_lags` term structure.",
        "i" = "Use a bare variable name inside each wrapper, e.g. `spatial(depth)`.",
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
    needs_time = any(terms_df$component %in% c("temporal", "spatiotemporal"))
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
      "`distributed_lags` terms wrapped in `temporal()` or `spatiotemporal()` require a `time` argument."
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

.normalize_dl_index <- function(x, n_max, name) {
  if (anyNA(x)) {
    cli_abort("`{name}` cannot contain `NA` values.")
  }
  x <- as.integer(x)
  if (!length(x)) return(x)
  min_x <- min(x)
  if (min_x == 0L) {
    x <- x + 1L
  } else if (min_x != 1L) {
    cli_abort("`{name}` must be 0-based or 1-based indexing.")
  }
  if (any(x < 1L | x > n_max)) {
    cli_abort("`{name}` contains indices outside valid range.")
  }
  x
}

.normalize_dl_year_i <- function(year_i, n_t = NULL) {
  if (anyNA(year_i)) {
    cli_abort("`year_i` cannot contain `NA` values.")
  }
  year_i <- as.integer(year_i)
  if (!length(year_i)) {
    if (is.null(n_t)) {
      cli_abort("`year_i` cannot be empty.")
    }
    return(list(year_i = integer(0), n_t = as.integer(n_t)))
  }
  min_year <- min(year_i)
  if (min_year == 1L) {
    year_i <- year_i - 1L
  } else if (min_year != 0L) {
    cli_abort("`year_i` must be 0-based or 1-based indexing.")
  }
  if (is.null(n_t)) {
    n_t <- max(year_i) + 1L
  } else {
    n_t <- as.integer(n_t)
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

  component_levels <- c("spatial", "temporal", "spatiotemporal")
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
