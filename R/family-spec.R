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
    cli_abort("Unsupported {what} supplied in `family`: {paste(bad, collapse = ', ')}")
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

.family_spec_response_family_id <- function(family_spec, y_i) {
  if (length(family_spec$family_id_i) != length(y_i)) {
    cli_abort("Internal family spec error: `family_id_i` must match the response length.")
  }
  family_spec$family_id_i
}

.family_spec_subset_rows <- function(family_spec, rows) {
  if (!length(family_spec$family_id_i)) {
    return(family_spec)
  }
  family_spec$family_id_i <- family_spec$family_id_i[rows]
  family_spec
}

.family_spec_validate_response <- function(y_i, family_spec, upr = NULL) {
  row_family <- .family_spec_response_family_id(family_spec, y_i)
  single_rows <- family_spec$combine_kind[row_family] == "single"
  family_name <- family_spec$family_name[cbind(row_family, 1L)]
  link_name <- family_spec$link_name[cbind(row_family, 1L)]

  positive_rows <- single_rows & family_name %in% c("Gamma", "lognormal")
  if (any(y_i[positive_rows] <= 0, na.rm = TRUE)) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }

  log_link_rows <- single_rows & link_name == "log"
  if (any(y_i[log_link_rows] < 0, na.rm = TRUE)) {
    cli_abort("`link = 'log'` but the response data include values < 0.")
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
    cli_abort("Unknown family names in `distribution_column`: {paste(unknown, collapse = ', ')}")
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
      cli_abort(
        "Families ending in `_mix` are not supported in multi-family mode: {paste(family_labels[has_mix], collapse = ', ')}"
      )
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
  row_family <- .family_spec_response_family_id(family_spec, y_i)
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
      cli_abort("{family_label} rows must have numeric response values in multi-family models.")
    }
    if (any(y_vals < 0)) {
      cli_abort("{family_label} rows must have non-negative response values in multi-family models.")
    }
    if (!allow_counts && any(y_vals > 1)) {
      cli_abort("Binomial rows must have values between 0 and 1 in multi-family models.")
    }

    counts_rows <- if (allow_counts) rows & !is.na(y_i) & y_i > 1 else rep(FALSE, length(rows))
    prop_rows <- rows & !is.na(y_i) & y_i > 0 & y_i < 1
    bernoulli_rows <- rows & !is.na(y_i) & y_i %in% c(0, 1)

    # Bernoulli rows do not require explicit trial sizes. If a shared
    # `weights` vector is present, default missing Bernoulli entries to 1.
    if (any(bernoulli_rows) && !is.null(weights)) {
      weights_vec <- weights
      missing_bernoulli_weights <- bernoulli_rows & is.na(weights_vec)
      if (any(missing_bernoulli_weights)) {
        weights_vec[missing_bernoulli_weights] <- 1
      }
      weights <- weights_vec
    }

    if (any(counts_rows) || any(prop_rows)) {
      if (is.null(weights)) {
        suffix <- if (allow_counts) "proportions or counts" else "proportions"
        cli_abort(
          "{family_label} rows with {suffix} require `weights` to supply the binomial size in multi-family models."
        )
      }
      weights_vec <- weights
      if (anyNA(weights_vec[counts_rows | prop_rows])) {
        cli_abort("`weights` must not contain missing values for {tolower(family_label)} rows in multi-family models.")
      }
      if (any(weights_vec[counts_rows | prop_rows] <= 0)) {
        cli_abort("`weights` must be > 0 for {tolower(family_label)} rows in multi-family models.")
      }
      if (any(counts_rows)) {
        if (any(weights_vec[counts_rows] < y_i[counts_rows], na.rm = TRUE)) {
          cli_abort("{family_label} counts must be <= `weights` (binomial size) in multi-family models.")
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
  row_family <- .family_spec_response_family_id(family_spec, y_i)
  y_out <- matrix(NA_real_, nrow = length(y_i), ncol = family_spec$n_m)

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

.object_family_spec <- function(object, caller = "This method") {
  if (!is.null(object$family_spec)) {
    return(object$family_spec)
  }

  if (.is_named_family_list(object$family) && length(object$family) > 1L) {
    cli_abort(
      "{caller} requires canonical `family_spec` metadata for multi-family objects. Refit this model with the current version of sdmTMB."
    )
  }

  .build_family_spec(object$family, data = object$data)
}

.family_spec_is_multi_family <- function(family_spec) {
  family_spec$n_f > 1L
}

.family_spec_has_two_components <- function(family_spec) {
  family_spec$n_m == 2L
}

.object_is_multi_family <- function(object, caller = "This method") {
  .family_spec_is_multi_family(.object_family_spec(object, caller = caller))
}

.object_has_two_components <- function(object, caller = "This method") {
  .family_spec_has_two_components(.object_family_spec(object, caller = caller))
}

.family_spec_row_family_id <- function(family_spec, data) {
  if (.family_spec_is_multi_family(family_spec)) {
    if (is.null(data)) {
      cli_abort("`newdata` is required to resolve row-wise families for this multi-family model.")
    }
    return(.validate_distribution_column(
      data = data,
      distribution_column = family_spec$distribution_column,
      family_labels = family_spec$family_labels
    ))
  }
  rep.int(1L, nrow(data))
}

.family_spec_inverse_link <- function(eta, link) {
  switch(link,
    identity = eta,
    log = exp(eta),
    logit = stats::plogis(eta),
    inverse = 1 / eta,
    cloglog = 1 - exp(-exp(eta)),
    cli_abort("Unsupported link in family spec: {.val {link}}")
  )
}

.family_spec_link <- function(mu, link) {
  switch(link,
    identity = mu,
    log = log(mu),
    logit = stats::qlogis(mu),
    inverse = 1 / mu,
    cloglog = log(-log1p(-mu)),
    cli_abort("Unsupported link in family spec: {.val {link}}")
  )
}

.family_spec_apply_link <- function(x, link, inverse = TRUE) {
  if (!length(x)) {
    return(x)
  }
  out <- x
  link_vals <- unique(link[!is.na(link)])
  for (this_link in link_vals) {
    ii <- link == this_link
    out[ii] <- if (inverse) {
      .family_spec_inverse_link(x[ii], this_link)
    } else {
      .family_spec_link(x[ii], this_link)
    }
  }
  out
}

.family_spec_prediction_output <- function(x, family_spec, row_family_id,
  type = c("link", "response"), model = NA_integer_, simulated = FALSE) {

  type <- match.arg(type)
  x <- as.matrix(x)
  n <- nrow(x)
  n_m <- family_spec$n_m
  if (ncol(x) < n_m) {
    cli_abort("Internal family spec error: prediction matrix has fewer components than expected.")
  }
  active <- family_spec$active[row_family_id, , drop = FALSE]
  combine_kind <- family_spec$combine_kind[row_family_id]
  link1 <- family_spec$link_name[cbind(row_family_id, 1L)]
  link2 <- if (n_m > 1L) family_spec$link_name[cbind(row_family_id, 2L)] else rep(NA_character_, n)
  raw1 <- x[, 1L]
  raw2 <- if (n_m > 1L) x[, 2L] else rep(NA_real_, n)

  if (simulated) {
    est1 <- raw1
    est2 <- if (n_m > 1L) raw2 else rep(NA_real_, n)
    if (n_m > 1L) est2[!active[, 2L]] <- NA_real_
    combined <- est1
    if (n_m > 1L) {
      two_component_rows <- combine_kind %in% c("delta", "poisson_link_delta")
      combined[two_component_rows] <- est1[two_component_rows] * est2[two_component_rows]
    }
  } else if (type == "response") {
    est1_raw <- .family_spec_apply_link(raw1, link1, inverse = TRUE)
    est2_raw <- rep(NA_real_, n)
    if (n_m > 1L && any(active[, 2L])) {
      est2_raw[active[, 2L]] <- .family_spec_apply_link(raw2[active[, 2L]], link2[active[, 2L]], inverse = TRUE)
    }
    est1 <- est1_raw
    est2 <- est2_raw
    poisson_link_rows <- if (n_m > 1L) combine_kind == "poisson_link_delta" else rep(FALSE, n)
    if (any(poisson_link_rows)) {
      n_groups <- est1_raw[poisson_link_rows]
      p_encounter <- 1 - exp(-n_groups)
      pos_mean <- est2_raw[poisson_link_rows]
      est1[poisson_link_rows] <- p_encounter
      est2[poisson_link_rows] <- (n_groups * pos_mean) / p_encounter
    }
    combined <- est1
    if (n_m > 1L) {
      delta_rows <- combine_kind == "delta"
      if (any(delta_rows)) {
        combined[delta_rows] <- est1[delta_rows] * est2[delta_rows]
      }
      if (any(poisson_link_rows)) {
        combined[poisson_link_rows] <- exp(raw1[poisson_link_rows] + raw2[poisson_link_rows])
      }
    }
  } else {
    est1 <- raw1
    est2 <- rep(NA_real_, n)
    if (n_m > 1L && any(active[, 2L])) {
      est2[active[, 2L]] <- raw2[active[, 2L]]
    }
    combined <- est1
    if (n_m > 1L) {
      delta_rows <- combine_kind == "delta"
      if (any(delta_rows)) {
        mu_prod <- .family_spec_apply_link(raw1[delta_rows], link1[delta_rows], inverse = TRUE) *
          .family_spec_apply_link(raw2[delta_rows], link2[delta_rows], inverse = TRUE)
        combined[delta_rows] <- .family_spec_apply_link(mu_prod, link2[delta_rows], inverse = FALSE)
      }
      poisson_link_rows <- combine_kind == "poisson_link_delta"
      if (any(poisson_link_rows)) {
        combined[poisson_link_rows] <- raw1[poisson_link_rows] + raw2[poisson_link_rows]
      }
    }
  }

  est <- if (is.na(model)) {
    combined
  } else if (isTRUE(model == 1L)) {
    est1
  } else if (isTRUE(model == 2L)) {
    est2
  } else {
    cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
  }

  list(est = est, est1 = est1, est2 = est2)
}

.family_spec_prediction_link_name <- function(family_spec, row_family_id, model = NA_integer_, simulated = FALSE) {
  if (simulated) {
    return("response")
  }
  link1 <- family_spec$link_name[cbind(row_family_id, 1L)]
  if (family_spec$n_m == 1L) {
    links <- link1
  } else if (is.na(model)) {
    link2 <- family_spec$link_name[cbind(row_family_id, 2L)]
    combine_kind <- family_spec$combine_kind[row_family_id]
    links <- ifelse(
      combine_kind == "single",
      link1,
      ifelse(combine_kind == "poisson_link_delta", "log", link2)
    )
  } else if (isTRUE(model == 1L)) {
    links <- link1
  } else if (isTRUE(model == 2L)) {
    active2 <- family_spec$active[row_family_id, 2L]
    link2 <- family_spec$link_name[cbind(row_family_id, 2L)]
    links <- link2[active2]
  } else {
    cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
  }
  links <- unique(stats::na.omit(links))
  if (length(links) == 1L) {
    links
  } else {
    "mixed"
  }
}
