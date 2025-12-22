.validate_multi_family_list <- function(family, data = NULL, distribution_column = NULL) {
  if (!is.list(family)) {
    cli_abort("`family` must be a named list of family objects for multi-family models.")
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
      "Families ending in `_mix` are not supported in multi-family mode: ",
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
  delta_family <- NULL
  poisson_link_delta <- NULL
  family_names2 <- NULL
  link_names2 <- NULL
  family_enum1 <- NULL
  family_enum2 <- NULL
  link_enum1 <- NULL
  link_enum2 <- NULL
  has_delta <- FALSE

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
    family_names2 <- vapply(
      family_list,
      function(x) if (length(x$family) > 1L) x$family[2] else NA_character_,
      character(1)
    )
    link_names2 <- vapply(
      family_list,
      function(x) if (length(x$link) > 1L) x$link[2] else NA_character_,
      character(1)
    )
    delta_family <- vapply(family_list, function(x) isTRUE(x$delta), logical(1))
    poisson_link_delta <- vapply(
      family_list,
      function(x) isTRUE(x$type == "poisson_link_delta"),
      logical(1)
    )
    has_delta <- any(delta_family)
    if (any(vapply(family_list, function(x) length(x$family) > 2L, logical(1)))) {
      cli_abort("Multi-family models with families that have more than 2 components are not supported.")
    }
    if (any(!delta_family & vapply(family_list, function(x) length(x$family) > 1L, logical(1)))) {
      cli_abort("Multi-component families are only supported for delta families.")
    }
    if (any(delta_family & vapply(family_list, function(x) length(x$family) != 2L, logical(1)))) {
      cli_abort("Delta families must include exactly 2 components.")
    }
    if (any(delta_family & family_names != "binomial")) {
      cli_abort("Delta families must use binomial for the first component.")
    }
    allowed_families <- c(
      "gaussian", "poisson", "binomial",
      "nbinom1", "nbinom2", "Gamma",
      "lognormal", "tweedie", "student",
      "gengamma", "Beta", "betabinomial"
    )
    family_names_all <- c(family_names, family_names2)
    family_names_all <- family_names_all[!is.na(family_names_all)]
    if (any(!family_names_all %in% allowed_families)) {
      bad <- family_names_all[!family_names_all %in% allowed_families]
      cli_abort(paste0(
        "Unsupported families in multi-family models: ",
        paste(unique(bad), collapse = ", ")
      ))
    }
    e_i <- as.integer(multi_info$e_i) - 1L
    family_enum1 <- as.integer(vapply(multi_info$family_enum, function(x) x[1], numeric(1)))
    link_enum1 <- as.integer(vapply(multi_info$link_enum, function(x) x[1], numeric(1)))
    family_enum2 <- vapply(
      seq_along(multi_info$family_enum),
      function(i) if (delta_family[i]) multi_info$family_enum[[i]][2] else -1,
      numeric(1)
    )
    link_enum2 <- vapply(
      seq_along(multi_info$link_enum),
      function(i) if (delta_family[i]) multi_info$link_enum[[i]][2] else -1,
      numeric(1)
    )
    family_enum2 <- as.integer(family_enum2)
    link_enum2 <- as.integer(link_enum2)
    family <- family_list[[1]]
  } else {
    if (!is.null(distribution_column)) {
      cli_abort("`distribution_column` is only supported for multi-family models (named `family` list).")
    }
  }

  list(
    multi_family = multi_family,
    family_list = family_list,
    multi_info = multi_info,
    family_names = family_names,
    link_names = link_names,
    family_names2 = family_names2,
    link_names2 = link_names2,
    delta_family = delta_family,
    poisson_link_delta = poisson_link_delta,
    family_enum1 = family_enum1,
    family_enum2 = family_enum2,
    link_enum1 = link_enum1,
    link_enum2 = link_enum2,
    has_delta = has_delta,
    e_i = e_i,
    family = family,
    family_input = if (multi_family) family_list else family
  )
}

.multi_family_param_offsets <- function(family_list) {
  n_fam <- length(family_list)
  family_names <- vapply(
    family_list,
    function(x) if (isTRUE(x$delta)) x$family[2] else x$family[1],
    character(1)
  )

  student_fixed <- vapply(
    family_list,
    function(x) {
      target_family <- if (isTRUE(x$delta)) x$family[2] else x$family[1]
      identical(target_family, "student") && !is.null(x$df)
    },
    logical(1)
  )
  if (any(student_fixed)) {
    cli_abort("Fixed student df is not supported in multi-family models yet.")
  }

  uses_phi <- family_names %in% c(
    "gaussian", "Gamma", "lognormal", "nbinom1", "nbinom2",
    "tweedie", "student", "gengamma", "Beta", "betabinomial"
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

.multi_family_build_response <- function(y_i, e_i, delta_family) {
  idx <- e_i + 1L
  delta_rows <- delta_family[idx]
  y1 <- ifelse(delta_rows, ifelse(y_i > 0, 1, 0), y_i)
  y2 <- ifelse(delta_rows & y_i > 0, y_i, NA_real_)
  cbind(y1, y2)
}

.multi_family_process_response <- function(
  y_i,
  size,
  weights,
  e_i,
  family_names,
  link_names,
  delta_family
) {
  e_i_idx <- e_i + 1L
  non_delta_rows <- !delta_family[e_i_idx]
  binom_rows <- non_delta_rows & family_names[e_i_idx] == "binomial"
  betabinom_rows <- non_delta_rows & family_names[e_i_idx] == "betabinomial"

  process_binomial_like <- function(rows, family_label, allow_counts) {
    if (!any(rows)) return(list(y_i = y_i, size = size, weights = weights))
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

  out <- process_binomial_like(binom_rows, "Binomial", allow_counts = FALSE)
  y_i <- out$y_i
  size <- out$size
  weights <- out$weights

  out <- process_binomial_like(betabinom_rows, "Betabinomial", allow_counts = TRUE)
  y_i <- out$y_i
  size <- out$size
  weights <- out$weights

  log_link_rows <- non_delta_rows & link_names[e_i_idx] == "log"
  if (any(log_link_rows) && min(y_i[log_link_rows], na.rm = TRUE) < 0) {
    cli_abort("`link = 'log'` but the response data include values < 0.")
  }
  positive_only_rows <- non_delta_rows & family_names[e_i_idx] %in% c("Gamma", "lognormal")
  if (any(positive_only_rows) && min(y_i[positive_only_rows], na.rm = TRUE) <= 0) {
    cli_abort("Gamma and lognormal rows must have response values > 0.")
  }

  list(y_i = y_i, size = size, weights = weights)
}

.multi_family_predict_e_g <- function(object, newdata, distribution_column) {
  if (is.null(distribution_column)) {
    cli_abort("`distribution_column` must be supplied for multi-family predictions.")
  }
  if (!is.list(object$family) || inherits(object$family, "family")) {
    cli_abort("Multi-family predictions require a named family list on the fitted object.")
  }
  info <- .validate_multi_family_list(
    object$family,
    data = newdata,
    distribution_column = distribution_column
  )
  as.integer(info$e_i) - 1L
}

.multi_family_predict_est <- function(
  eta1,
  eta2,
  family_list,
  fam_index,
  delta_family,
  poisson_link_delta,
  type = c("link", "response")
) {
  type <- match.arg(type)
  n <- length(eta1)
  est <- eta1
  est1 <- rep(NA_real_, n)
  est2 <- rep(NA_real_, n)
  delta_rows <- delta_family[fam_index]

  if (type == "link") {
    est1 <- eta1
    if (!is.null(eta2)) est2[delta_rows] <- eta2[delta_rows]
    if (any(delta_rows)) {
      plink <- delta_rows & poisson_link_delta[fam_index]
      if (any(plink)) {
        est[plink] <- eta1[plink] + eta2[plink]
      }
      standard <- delta_rows & !plink
      if (any(standard)) {
        for (k in seq_along(family_list)) {
          idx <- standard & fam_index == k
          if (!any(idx)) next
          fam <- family_list[[k]]
          p1 <- fam[[1]]$linkinv(eta1[idx])
          p2 <- fam[[2]]$linkinv(eta2[idx])
          est[idx] <- fam[[2]]$linkfun(p1 * p2)
        }
      }
    }
  } else {
    for (k in seq_along(family_list)) {
      idx <- fam_index == k
      if (!any(idx)) next
      fam <- family_list[[k]]
      if (isTRUE(fam$delta)) {
        est1[idx] <- fam[[1]]$linkinv(eta1[idx])
        est2[idx] <- fam[[2]]$linkinv(eta2[idx])
      } else {
        est1[idx] <- fam$linkinv(eta1[idx])
      }
    }
    est <- est1
    if (any(delta_rows)) {
      plink <- delta_rows & poisson_link_delta[fam_index]
      if (any(plink)) {
        n_val <- est1[plink]
        p <- 1 - exp(-n_val)
        w <- est2[plink]
        r <- (n_val * w) / p
        est1[plink] <- p
        est2[plink] <- r
        est[plink] <- n_val * w
      }
      standard <- delta_rows & !plink
      if (any(standard)) {
        est[standard] <- est1[standard] * est2[standard]
      }
    }
  }

  list(est = est, est1 = est1, est2 = est2)
}
