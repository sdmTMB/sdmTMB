.validate_multi_family_list <- function(family, data = NULL, distribution_column = NULL) {
  if (!is.list(family)) {
    cli_abort("`family` must be a named list of family objects for multi-likelihood models.")
  }
  if (length(family) == 0L) {
    cli_abort("`family` must contain at least one family object.")
  }
  if (is.null(names(family)) || anyNA(names(family)) || any(!nzchar(names(family)))) {
    cli_abort("`family` must be a named list with non-empty names.")
  }
  if (anyDuplicated(names(family))) {
    cli_abort("`family` list names must be unique.")
  }
  is_family <- vapply(family, inherits, logical(1), "family")
  if (!all(is_family)) {
    cli_abort("All elements of `family` must inherit from class 'family'.")
  }

  has_mix <- vapply(
    family,
    function(x) any(grepl("_mix$", x$family)),
    logical(1)
  )
  if (any(has_mix)) {
    cli_abort(paste0(
      "Families ending in `_mix` are not supported in multi-likelihood mode: ",
      paste(names(family)[has_mix], collapse = ", ")
    ))
  }

  family_enum <- lapply(family, function(x) {
    codes <- .enum_family(x$family)
    if (any(is.na(codes))) {
      cli_abort("Unsupported family supplied in `family` list.")
    }
    codes
  })
  link_enum <- lapply(family, function(x) {
    codes <- .enum_link(x$link)
    if (any(is.na(codes))) {
      cli_abort("Unsupported link supplied in `family` list.")
    }
    codes
  })

  out <- list(
    family_names = names(family),
    family_enum = family_enum,
    link_enum = link_enum
  )

  if (!is.null(distribution_column)) {
    if (is.null(data)) {
      cli_abort("`data` must be supplied when `distribution_column` is provided.")
    }
    if (!distribution_column %in% names(data)) {
      cli_abort("`distribution_column` must be a column in `data`.")
    }
    dist_values <- data[[distribution_column]]
    if (is.factor(dist_values)) {
      dist_values <- as.character(dist_values)
    }
    if (!is.character(dist_values)) {
      cli_abort("`distribution_column` must be a character or factor column.")
    }
    if (anyNA(dist_values)) {
      cli_abort("`distribution_column` must not contain missing values.")
    }
    missing_names <- setdiff(unique(dist_values), names(family))
    if (length(missing_names) > 0L) {
      cli_abort(paste0(
        "Unknown family names in `distribution_column`: ",
        paste(missing_names, collapse = ", ")
      ))
    }
    out$e_i <- match(dist_values, names(family))
  }

  out
}
