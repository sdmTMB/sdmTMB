.process_binomial_like_rows <- function(y_i, size, weights, rows,
  family_label, allow_counts, weighted_binary_counts = FALSE) {
  if (!any(rows)) return(list(y_i = y_i, size = size, weights = weights))
  values <- y_i[rows & !is.na(y_i)]
  if (!is.numeric(values) || any(values < 0)) {
    cli_abort("{family_label} rows must have non-negative numeric response values.")
  }

  classified <- .classify_binomial_like_numeric_rows(
    y_i, weights, allow_counts, weighted_binary_counts
  )
  counts <- rows & classified$count_rows
  bernoulli <- rows & classified$bernoulli_rows
  proportions <- rows & classified$proportion_rows

  # A numeric binomial response that contains proportions uses `weights` as
  # the number of trials. Binary rows with an available weight are then
  # boundary proportions (e.g. 3/3), which keeps the fit identical to an
  # equivalent cbind(success, failure) representation. Binary rows with a
  # missing weight stay Bernoulli, and the all-binary case keeps `weights`
  # as likelihood weights.
  if (!allow_counts && any(proportions) && !is.null(weights)) {
    weighted_binary <- bernoulli & !is.na(weights)
    proportions <- proportions | weighted_binary
    bernoulli <- bernoulli & !weighted_binary
  }

  if (!allow_counts && any(rows & !is.na(y_i) & y_i > 1)) {
    cli_abort("Binomial rows must have values between 0 and 1.")
  }
  if (allow_counts && any(proportions & !(y_i > 0 & y_i < 1))) {
    cli_abort("{family_label} rows must be integer counts or proportions in (0, 1).")
  }
  if (any(bernoulli) && !is.null(weights)) weights[bernoulli & is.na(weights)] <- 1

  converted <- counts | proportions
  if (any(converted)) {
    if (is.null(weights)) {
      kind <- if (allow_counts) "proportions or counts" else "proportions"
      cli_abort("{family_label} rows with {kind} require `weights` to supply binomial size.")
    }
    missing_size <- converted & is.na(weights)
    if (any(missing_size)) {
      # Match model-frame NA handling: a missing trial size makes the
      # observation unavailable to the likelihood, just like a missing
      # response. Use harmless placeholders for TMB's auxiliary vectors.
      y_i[missing_size] <- NA_real_
      size[missing_size] <- 1
      weights[missing_size] <- 1
      converted[missing_size] <- FALSE
    }
    if (any(weights[converted] <= 0)) {
      cli_abort("`weights` must be non-missing and > 0 for {tolower(family_label)} rows.")
    }
    if (any(counts & weights < y_i, na.rm = TRUE)) {
      cli_abort("{family_label} counts must be <= `weights` (binomial size).")
    }
    size[converted] <- weights[converted]
    y_i[proportions] <- y_i[proportions] * weights[proportions]
    weights[converted] <- 1
  }
  list(y_i = y_i, size = size, weights = weights)
}

.family_spec_process_response <- function(y_i, size, weights, family_spec) {
  row_family <- .family_spec_response_family_id(family_spec, y_i)
  component1 <- family_spec$components[family_spec$components$component == 1L, ]
  family_name <- component1$family_name[match(row_family, component1$family_id)]
  single <- family_spec$families$combine_kind[row_family] == "single"

  out <- .process_binomial_like_rows(
    y_i, size, weights, single & family_name == "binomial",
    "Binomial", allow_counts = FALSE
  )
  .process_binomial_like_rows(
    out$y_i, out$size, out$weights, single & family_name == "betabinomial",
    "Betabinomial", allow_counts = TRUE, weighted_binary_counts = TRUE
  )
}

.family_spec_validate_response <- function(y_i, family_spec, upr = NULL) {
  row_family <- .family_spec_response_family_id(family_spec, y_i)
  component1 <- family_spec$components[family_spec$components$component == 1L, ]
  family_name <- component1$family_name[match(row_family, component1$family_id)]
  link_name <- component1$link_name[match(row_family, component1$family_id)]
  single <- family_spec$families$combine_kind[row_family] == "single"

  if (any(y_i[single & family_name %in% c("Gamma", "lognormal")] <= 0, na.rm = TRUE)) {
    cli_abort("Gamma and lognormal must have response values > 0.")
  }
  ordered_beta <- single & family_name == "ordbeta"
  if (any(y_i[ordered_beta] < 0 | y_i[ordered_beta] > 1, na.rm = TRUE)) {
    cli_abort("Ordered beta requires response values in [0, 1].")
  }
  if (any(y_i[single & link_name == "log"] < 0, na.rm = TRUE)) {
    cli_abort("`link = 'log'` but the response data include values < 0.")
  }
  censored <- single & family_name == "censored_poisson"
  if (!is.null(upr) && any(y_i[censored] > upr[censored], na.rm = TRUE)) {
    cli_abort("Observed values must be <= `control$censored_upper` for censored Poisson rows.")
  }
  invisible(NULL)
}

.family_spec_build_response <- function(y_i, family_spec) {
  row_family <- .family_spec_response_family_id(family_spec, y_i)
  response <- matrix(NA_real_, nrow = length(y_i), ncol = family_spec$n_m)
  single <- family_spec$families$combine_kind[row_family] == "single"
  response[single, 1L] <- y_i[single]
  response[!single, 1L] <- as.numeric(y_i[!single] > 0)
  if (family_spec$n_m == 2L) {
    positive <- !single & y_i > 0
    response[positive, 2L] <- y_i[positive]
  }
  response
}

.prepare_family_response <- function(y_i, weights, family_spec, upr = NULL) {
  component1 <- .family_spec_component_value(family_spec, 1L, 1L, "family_name")
  ordinary_binomial_like <- family_spec$n_f == 1L && family_spec$n_m == 1L &&
    component1 %in% c("binomial", "betabinomial")
  if (ordinary_binomial_like && (is.character(y_i) || is.factor(y_i))) {
    y_i <- factor(y_i)
    if (nlevels(y_i) > 2L) cli_abort("More than 2 levels detected for response")
    y_i <- pmin(as.numeric(y_i) - 1, 1)
  } else if (ordinary_binomial_like && is.matrix(y_i)) {
    if (ncol(y_i) != 2L) cli_abort("Binomial matrix responses must have two columns.")
    size <- rowSums(y_i)
    y_i <- y_i[, 1L]
    .family_spec_validate_response(y_i, family_spec, upr)
    return(list(
      y_i = y_i, size = size, weights = weights,
      response = .family_spec_build_response(y_i, family_spec)
    ))
  }
  size <- rep(1, NROW(y_i))
  processed <- .family_spec_process_response(y_i, size, weights, family_spec)
  .family_spec_validate_response(processed$y_i, family_spec, upr)
  processed$response <- .family_spec_build_response(processed$y_i, family_spec)
  processed
}
