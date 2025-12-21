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

.parse_multi_family <- function(family, data, distribution_column) {
  multi_family <- is.list(family) && !inherits(family, "family")
  family_list <- NULL
  multi_info <- NULL
  family_names <- NULL
  link_names <- NULL
  e_i <- NULL

  if (multi_family) {
    if (is.null(distribution_column)) {
      cli_abort("`distribution_column` must be supplied when `family` is a named list.")
    }
    family_list <- family
    multi_info <- .validate_multi_family_list(
      family_list,
      data = data,
      distribution_column = distribution_column
    )
    family_names <- vapply(family_list, function(x) x$family[1], character(1))
    link_names <- vapply(family_list, function(x) x$link[1], character(1))
    if (any(vapply(family_list, function(x) length(x$family) > 1L, logical(1)))) {
      cli_abort("Multi-likelihood models do not support multi-component families yet.")
    }
    if (any(vapply(family_list, function(x) isTRUE(x$delta), logical(1)))) {
      cli_abort("Delta families are not supported in multi-likelihood models yet.")
    }
    allowed_families <- c(
      "gaussian", "poisson", "binomial",
      "nbinom1", "nbinom2", "Gamma",
      "lognormal", "tweedie", "student",
      "gengamma"
    )
    if (any(!family_names %in% allowed_families)) {
      bad <- family_names[!family_names %in% allowed_families]
      cli_abort(paste0(
        "Only gaussian, poisson, and binomial are supported in multi-likelihood models for now. ",
        "Unsupported families: ", paste(unique(bad), collapse = ", ")
      ))
    }
    e_i <- as.integer(multi_info$e_i) - 1L
    family <- family_list[[1]]
  } else {
    if (!is.null(distribution_column)) {
      cli_abort("`distribution_column` is reserved for multi-likelihood models and is not yet supported.")
    }
  }

  list(
    multi_family = multi_family,
    family_list = family_list,
    multi_info = multi_info,
    family_names = family_names,
    link_names = link_names,
    e_i = e_i,
    family = family,
    family_input = if (multi_family) family_list else family
  )
}

.multi_family_param_offsets <- function(family_list) {
  n_fam <- length(family_list)
  family_names <- vapply(family_list, function(x) x$family[1], character(1))

  student_fixed <- vapply(
    family_list,
    function(x) identical(x$family[1], "student") && !is.null(x$df),
    logical(1)
  )
  if (any(student_fixed)) {
    cli_abort("Fixed student df is not supported in multi-likelihood models yet.")
  }

  uses_phi <- family_names %in% c(
    "gaussian", "Gamma", "lognormal", "nbinom1", "nbinom2",
    "tweedie", "student", "gengamma"
  )
  uses_thetaf <- family_names %in% "tweedie"
  uses_student_df <- family_names %in% "student"
  uses_gengamma_Q <- family_names %in% "gengamma"

  make_offsets <- function(uses) {
    len <- ifelse(uses, 1L, 0L)
    start <- rep(-1L, n_fam)
    if (any(uses)) {
      start[uses] <- seq_len(sum(uses)) - 1L
    }
    list(start = start, len = len, total = sum(len))
  }

  list(
    ln_phi = make_offsets(uses_phi),
    thetaf = make_offsets(uses_thetaf),
    ln_student_df = make_offsets(uses_student_df),
    gengamma_Q = make_offsets(uses_gengamma_Q)
  )
}
