.is_named_family_list <- function(family) {
  is.list(family) && !inherits(family, "family")
}

.validate_family_spec_list <- function(family_list) {
  if (!is.list(family_list) || inherits(family_list, "family")) {
    cli_abort("`family` must be a family object or a named list of family objects.")
  }
  if (length(family_list) == 0L) {
    cli_abort("`family` must contain at least one family object.")
  }
  if (is.null(names(family_list)) || anyNA(names(family_list)) || any(!nzchar(names(family_list)))) {
    cli_abort("`family` must be a named list with non-empty names.")
  }
  if (anyDuplicated(names(family_list))) {
    cli_abort("`family` list names must be unique.")
  }
  is_family <- vapply(family_list, inherits, logical(1), what = "family")
  if (!all(is_family)) {
    cli_abort("All elements of `family` must inherit from class 'family'.")
  }
  invisible(family_list)
}

.map_family_codes <- function(values, valid, what) {
  codes <- unname(valid[values])
  if (anyNA(codes)) {
    bad <- unique(values[is.na(codes)])
    cli_abort(paste0(
      "Unsupported ", what, " supplied in `family`: ",
      paste(bad, collapse = ", ")
    ))
  }
  as.integer(codes)
}

.make_family_param_slots <- function(family_names) {
  slot_index <- function(uses) {
    out <- rep(NA_integer_, length(uses))
    if (any(uses)) {
      out[uses] <- seq_len(sum(uses))
    }
    out
  }

  uses_phi <- !family_names %in% c("binomial", "poisson", "censored_poisson")

  list(
    ln_phi = slot_index(uses_phi),
    thetaf = slot_index(family_names == "tweedie"),
    ln_student_df = slot_index(family_names == "student"),
    gengamma_Q = slot_index(family_names == "gengamma")
  )
}

.family_spec_slot_length <- function(slot) {
  used <- slot[!is.na(slot)]
  if (!length(used)) {
    return(0L)
  }
  as.integer(max(used))
}

.family_spec_tmb_data <- function(family_spec) {
  zero_based_slot <- function(slot) {
    out <- rep(-1L, length(slot))
    used <- !is.na(slot)
    out[used] <- slot[used] - 1L
    out
  }

  family_code <- matrix(0L, nrow = family_spec$n_f, ncol = family_spec$n_m)
  link_code <- matrix(0L, nrow = family_spec$n_f, ncol = family_spec$n_m)
  family_code[family_spec$active] <- family_spec$family_code[family_spec$active]
  link_code[family_spec$active] <- family_spec$link_code[family_spec$active]

  combine_kind_codes <- c(
    single = 0L,
    delta = 1L,
    poisson_link_delta = 2L
  )

  list(
    obs_family_id = family_spec$family_id_i - 1L,
    component_active = matrix(
      as.integer(family_spec$active),
      nrow = family_spec$n_f,
      ncol = family_spec$n_m
    ),
    family_code = family_code,
    link_code = link_code,
    combine_kind = unname(as.integer(combine_kind_codes[family_spec$combine_kind])),
    ln_phi_slot = zero_based_slot(family_spec$param_slot$ln_phi),
    thetaf_slot = zero_based_slot(family_spec$param_slot$thetaf),
    ln_student_df_slot = zero_based_slot(family_spec$param_slot$ln_student_df),
    gengamma_Q_slot = zero_based_slot(family_spec$param_slot$gengamma_Q)
  )
}

.family_spec_validate_response <- function(y_i, family_spec, upr = NULL) {
  if (length(family_spec$family_id_i) != length(y_i)) {
    cli_abort("Internal family spec error: `family_id_i` must match the response length.")
  }

  row_family <- family_spec$family_id_i
  single_rows <- family_spec$combine_kind[row_family] == "single"
  family_name <- family_spec$family_name[cbind(row_family, 1L)]
  link_name <- family_spec$link_name[cbind(row_family, 1L)]

  positive_rows <- single_rows & family_name %in% c("Gamma", "lognormal")
  if (any(y_i[positive_rows] <= 0, na.rm = TRUE)) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }

  log_link_rows <- single_rows & link_name == "log"
  if (any(y_i[log_link_rows] < 0, na.rm = TRUE)) {
    cli_abort("`link = 'log'` but the reponse data include values < 0.")
  }

  if (!is.null(upr)) {
    censored_rows <- single_rows & family_name == "censored_poisson"
    if (any(y_i[censored_rows] > upr[censored_rows], na.rm = TRUE)) {
      cli_abort("Observed values must be <= `control$censored_upper` for censored Poisson rows.")
    }
  }

  invisible(NULL)
}

.validate_distribution_column <- function(data, distribution_column, family_labels) {
  if (is.null(data)) {
    cli_abort("`data` must be supplied when `distribution_column` is used.")
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
  unknown <- setdiff(unique(dist_values), family_labels)
  if (length(unknown) > 0L) {
    cli_abort(paste0(
      "Unknown family names in `distribution_column`: ",
      paste(unknown, collapse = ", ")
    ))
  }
  match(dist_values, family_labels)
}

.build_family_spec <- function(family, data = NULL, distribution_column = NULL) {
  if (inherits(family, "family")) {
    if (!is.null(distribution_column)) {
      cli_abort("`distribution_column` is only supported for named `family` lists.")
    }
    family_list <- list(family1 = family)
    user_family <- family
  } else if (.is_named_family_list(family)) {
    .validate_family_spec_list(family)
    family_list <- family
    user_family <- family
  } else {
    cli_abort("`family` must be a family object or a named list of family objects.")
  }

  n_f <- length(family_list)
  family_labels <- names(family_list)
  components_per_family <- vapply(family_list, function(x) length(x$family), integer(1))
  if (any(components_per_family > 2L)) {
    cli_abort("Families with more than 2 components are not supported.")
  }

  is_delta_family <- vapply(family_list, function(x) isTRUE(x$delta), logical(1))
  has_two_components <- components_per_family == 2L
  if (any(has_two_components & !is_delta_family)) {
    cli_abort("Only delta families can have 2 components.")
  }
  if (any(is_delta_family & components_per_family != 2L)) {
    cli_abort("Delta families must include exactly 2 components.")
  }

  if (n_f > 1L) {
    has_mix <- vapply(
      family_list,
      function(x) any(grepl("_mix$", x$family)),
      logical(1)
    )
    if (any(has_mix)) {
      cli_abort(paste0(
        "Families ending in `_mix` are not supported in multi-family mode: ",
        paste(family_labels[has_mix], collapse = ", ")
      ))
    }
  }

  n_m <- max(components_per_family)
  active <- matrix(FALSE, nrow = n_f, ncol = n_m)
  family_name <- matrix(NA_character_, nrow = n_f, ncol = n_m)
  link_name <- matrix(NA_character_, nrow = n_f, ncol = n_m)
  for (i in seq_len(n_f)) {
    active[i, seq_len(components_per_family[i])] <- TRUE
    family_name[i, seq_len(components_per_family[i])] <- family_list[[i]]$family
    link_name[i, seq_len(components_per_family[i])] <- family_list[[i]]$link
  }

  family_code <- matrix(NA_integer_, nrow = n_f, ncol = n_m)
  link_code <- matrix(NA_integer_, nrow = n_f, ncol = n_m)
  family_code[active] <- .map_family_codes(family_name[active], .valid_family, "family")
  link_code[active] <- .map_family_codes(link_name[active], .valid_link, "link")

  combine_kind <- ifelse(
    has_two_components,
    ifelse(
      vapply(family_list, function(x) identical(x$type, "poisson_link_delta"), logical(1)),
      "poisson_link_delta",
      "delta"
    ),
    "single"
  )

  target_family <- ifelse(has_two_components, family_name[, 2], family_name[, 1])
  fixed_student_df <- vapply(
    family_list,
    function(x) {
      if (length(x$family) == 2L) {
        identical(x$family[2], "student") && !is.null(x$df)
      } else {
        identical(x$family[1], "student") && !is.null(x$df)
      }
    },
    logical(1)
  )
  if (n_f > 1L && any(fixed_student_df)) {
    cli_abort("Fixed student df is not supported in multi-family models yet.")
  }

  param_slot <- .make_family_param_slots(target_family)

  if (n_f > 1L) {
    if (is.null(distribution_column)) {
      cli_abort("`distribution_column` must be supplied when `family` has more than one entry.")
    }
    family_id_i <- .validate_distribution_column(data, distribution_column, family_labels)
  } else {
    if (!is.null(distribution_column)) {
      cli_abort("`distribution_column` is only supported when `family` has more than one entry.")
    }
    family_id_i <- if (is.null(data)) integer(0) else rep.int(1L, nrow(data))
  }

  list(
    n_f = n_f,
    n_m = n_m,
    family_list = family_list,
    family_labels = family_labels,
    distribution_column = if (n_f > 1L) distribution_column else NULL,
    family_id_i = as.integer(family_id_i),
    active = active,
    family_name = family_name,
    link_name = link_name,
    family_code = family_code,
    link_code = link_code,
    combine_kind = combine_kind,
    param_slot = param_slot,
    family = family_list[[1]],
    family_input = user_family
  )
}

.family_spec_process_response <- function(y_i, size, weights, family_spec) {
  if (length(family_spec$family_id_i) != length(y_i)) {
    cli_abort("Internal family spec error: `family_id_i` must match the response length.")
  }
  row_family <- family_spec$family_id_i
  family_name <- family_spec$family_name[cbind(row_family, 1L)]
  non_delta_rows <- family_spec$combine_kind[row_family] == "single"
  binom_rows <- non_delta_rows & family_name == "binomial"
  betabinom_rows <- non_delta_rows & family_name == "betabinomial"

  process_binomial_like <- function(rows, family_label, allow_counts) {
    if (!any(rows)) {
      return(list(y_i = y_i, size = size, weights = weights))
    }
    y_vals <- y_i[rows]
    y_vals <- y_vals[!is.na(y_vals)]
    if (!is.numeric(y_vals)) {
      cli_abort(paste0(family_label, " rows must have numeric response values in multi-family models."))
    }
    if (any(y_vals < 0)) {
      cli_abort(paste0(family_label, " rows must have non-negative response values in multi-family models."))
    }
    if (!allow_counts && any(y_vals > 1)) {
      cli_abort("Binomial rows must have values between 0 and 1 in multi-family models.")
    }

    counts_rows <- if (allow_counts) rows & !is.na(y_i) & y_i > 1 else rep(FALSE, length(rows))
    prop_rows <- rows & !is.na(y_i) & y_i > 0 & y_i < 1
    if (any(counts_rows) || any(prop_rows)) {
      if (is.null(weights)) {
        suffix <- if (allow_counts) "proportions or counts" else "proportions"
        cli_abort(paste0(
          family_label, " rows with ", suffix,
          " require `weights` to supply the binomial size in multi-family models."
        ))
      }
      weights_vec <- weights
      if (anyNA(weights_vec[counts_rows | prop_rows])) {
        cli_abort(paste0(
          "`weights` must not contain missing values for ",
          tolower(family_label), " rows in multi-family models."
        ))
      }
      if (any(weights_vec[counts_rows | prop_rows] <= 0)) {
        cli_abort(paste0(
          "`weights` must be > 0 for ",
          tolower(family_label), " rows in multi-family models."
        ))
      }
      if (any(counts_rows)) {
        if (any(weights_vec[counts_rows] < y_i[counts_rows], na.rm = TRUE)) {
          cli_abort(paste0(
            family_label, " counts must be <= `weights` (binomial size) in multi-family models."
          ))
        }
        size[counts_rows] <- weights_vec[counts_rows]
        weights_vec[counts_rows] <- 1
      }
      if (any(prop_rows)) {
        size[prop_rows] <- weights_vec[prop_rows]
        y_i[prop_rows] <- y_i[prop_rows] * weights_vec[prop_rows]
        weights_vec[prop_rows] <- 1
      }
      weights <- weights_vec
    }

    list(y_i = y_i, size = size, weights = weights)
  }

  res <- process_binomial_like(binom_rows, family_label = "Binomial", allow_counts = FALSE)
  y_i <- res$y_i
  size <- res$size
  weights <- res$weights

  res <- process_binomial_like(betabinom_rows, family_label = "Betabinomial", allow_counts = TRUE)
  y_i <- res$y_i
  size <- res$size
  weights <- res$weights

  list(y_i = y_i, size = size, weights = weights)
}

.family_spec_build_response <- function(y_i, family_spec) {
  if (length(family_spec$family_id_i) != length(y_i)) {
    cli_abort("Internal family spec error: `family_id_i` must match the response length.")
  }
  y_out <- matrix(NA_real_, nrow = length(y_i), ncol = family_spec$n_m)
  row_family <- family_spec$family_id_i

  single_rows <- family_spec$combine_kind[row_family] == "single"
  if (any(single_rows)) {
    y_out[single_rows, 1L] <- y_i[single_rows]
  }

  two_component_rows <- !single_rows
  if (any(two_component_rows)) {
    y_out[two_component_rows, 1L] <- ifelse(y_i[two_component_rows] > 0, 1, 0)
    y_out[two_component_rows, 2L] <- ifelse(y_i[two_component_rows] > 0, y_i[two_component_rows], NA_real_)
  }

  y_out
}
