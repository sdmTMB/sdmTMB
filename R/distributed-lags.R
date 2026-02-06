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

.build_vertex_time_covariates <- function(...) {
  cli_abort("Internal error: `.build_vertex_time_covariates()` is not implemented yet.")
}

.build_distributed_lag_tmb_data <- function(...) {
  cli_abort("Internal error: `.build_distributed_lag_tmb_data()` is not implemented yet.")
}
