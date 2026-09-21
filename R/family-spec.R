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

# Family subsystem architecture:
# user input -> .compile_family_spec() -> .establish_analysis_rows() ->
# .prepare_family_response() -> .as_tmb_family_data() -> TMB resolver.
#
# Invariants: family IDs are one-based in R and zero-based only at the TMB
# boundary; a fit has one or two latent components; single families use LP1;
# delta families use LP1 and LP2; auxiliary parameters belong to families;
# and every analysis row has exactly one family ID. `families` and `components`
# are the sole declarative representation. Dense matrices and -1 sentinels are
# constructed only by `.as_tmb_family_data()`.

.family_registry <- local({
  registry <- data.frame(
    family_name = names(.valid_family),
    uses_phi = !names(.valid_family) %in% c("binomial", "poisson", "censored_poisson"),
    auxiliary = NA_character_,
    stringsAsFactors = FALSE
  )
  registry$auxiliary[registry$family_name == "tweedie"] <- "thetaf"
  registry$auxiliary[registry$family_name == "student"] <- "ln_student_df"
  registry$auxiliary[registry$family_name == "gengamma"] <- "gengamma_Q"
  function() registry
})

.family_metadata <- function(family_name) {
  registry <- .family_registry()
  out <- registry[match(family_name, registry$family_name), , drop = FALSE]
  if (anyNA(out$family_name)) {
    bad <- unique(family_name[is.na(out$family_name)])
    cli_abort("Unsupported family supplied in `family`: {paste(bad, collapse = ', ')}")
  }
  out
}

.make_family_param_slots <- function(family_names) {
  slot_index <- function(uses) {
    out <- rep(NA_integer_, length(uses))
    if (any(uses)) {
      out[uses] <- seq_len(sum(uses))
    }
    out
  }

  metadata <- .family_metadata(family_names)
  uses_phi <- metadata$uses_phi

  list(
    ln_phi = slot_index(uses_phi),
    thetaf = slot_index(metadata$auxiliary == "thetaf" & !is.na(metadata$auxiliary)),
    ln_student_df = slot_index(metadata$auxiliary == "ln_student_df" & !is.na(metadata$auxiliary)),
    gengamma_Q = slot_index(metadata$auxiliary == "gengamma_Q" & !is.na(metadata$auxiliary))
  )
}

.family_spec_slot_length <- function(slot) {
  used <- slot[!is.na(slot)]
  if (!length(used)) {
    return(0L)
  }
  as.integer(max(used))
}

.as_tmb_family_data <- function(family_spec) {
  zero_based_slot <- function(slot) {
    out <- rep(-1L, length(slot))
    used <- !is.na(slot)
    out[used] <- slot[used] - 1L
    out
  }

  component_active <- matrix(0L, nrow = family_spec$n_f, ncol = family_spec$n_m)
  family_code <- link_code <- component_active
  index <- cbind(family_spec$components$family_id, family_spec$components$component)
  component_active[index] <- 1L
  family_code[index] <- family_spec$components$family_code
  link_code[index] <- family_spec$components$link_code

  combine_kind_codes <- c(
    single = 0L,
    delta = 1L,
    poisson_link_delta = 2L
  )

  list(
    obs_family_id = family_spec$family_id_i - 1L,
    component_active = component_active,
    family_code = family_code,
    link_code = link_code,
    combine_kind = unname(as.integer(combine_kind_codes[family_spec$families$combine_kind])),
    ln_phi_slot = zero_based_slot(family_spec$param_slot$ln_phi),
    thetaf_slot = zero_based_slot(family_spec$param_slot$thetaf),
    ln_student_df_slot = zero_based_slot(family_spec$param_slot$ln_student_df),
    gengamma_Q_slot = zero_based_slot(family_spec$param_slot$gengamma_Q)
  )
}

.family_parameter_values <- function(family_spec, estimate_student_df,
  student_df_fixed) {
  n <- lapply(family_spec$param_slot, .family_spec_slot_length)
  has_ordbeta <- any(family_spec$components$family_name == "ordbeta")
  list(
    thetaf = rep(0, n$thetaf),
    ln_student_df = if (n$ln_student_df > 0L) {
      rep(if (estimate_student_df) log(2) else log(student_df_fixed - 1), n$ln_student_df)
    } else {
      numeric(0)
    },
    gengamma_Q = rep(0.5, n$gengamma_Q),
    psi = if (has_ordbeta) c(-1, 1) else numeric(0),
    ln_phi = rep(0, n$ln_phi)
  )
}

.map_family_parameters <- function(tmb_map, tmb_params, has_dispformula,
  estimate_student_df) {
  if (length(tmb_params$thetaf)) tmb_map$thetaf <- NULL
  if (length(tmb_params$ln_student_df) && estimate_student_df) {
    tmb_map$ln_student_df <- NULL
  }
  if (length(tmb_params$gengamma_Q)) tmb_map$gengamma_Q <- NULL
  if (length(tmb_params$psi)) tmb_map$psi <- NULL
  if (length(tmb_params$ln_phi)) {
    tmb_map$ln_phi <- if (has_dispformula) {
      factor(rep(NA, length(tmb_params$ln_phi)))
    } else {
      NULL
    }
  }
  tmb_map
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

.family_spec_component_active <- function(family_spec, row_family_id = family_spec$family_id_i) {
  if (!length(row_family_id)) {
    return(matrix(FALSE, nrow = 0L, ncol = family_spec$n_m))
  }
  if (anyNA(row_family_id) || any(row_family_id < 1L | row_family_id > family_spec$n_f)) {
    cli_abort("Internal family spec error: row-wise family ids are out of bounds.")
  }
  out <- matrix(FALSE, nrow = length(row_family_id), ncol = family_spec$n_m)
  for (component in seq_len(family_spec$n_m)) {
    active_families <- family_spec$components$family_id[
      family_spec$components$component == component
    ]
    out[, component] <- row_family_id %in% active_families
  }
  out
}

.family_spec_observed_response <- function(response, family_spec, model = NA_integer_) {
  response <- as.matrix(response)
  if (nrow(response) != length(family_spec$family_id_i)) {
    cli_abort("Internal family spec error: response rows do not match family metadata.")
  }
  if (ncol(response) < family_spec$n_m) {
    cli_abort("Internal family spec error: response columns do not match family metadata.")
  }

  active <- .family_spec_component_active(family_spec)
  if (is.na(model)) {
    out <- response[, 1L]
    if (family_spec$n_m > 1L) {
      use_model2 <- active[, 2L] & !is.na(response[, 2L])
      out[use_model2] <- response[use_model2, 2L]
    }
    return(as.numeric(out))
  }

  model <- as.integer(model)
  if (!model %in% seq_len(family_spec$n_m)) {
    cli_abort("`model` argument isn't valid for the fitted family structure.")
  }
  out <- response[, model]
  out[!active[, model]] <- NA_real_
  as.numeric(out)
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

.compile_family_spec <- function(family, data = NULL, distribution_column = NULL) {
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
    has_ordbeta <- vapply(
      family_list,
      function(x) any(x$family == "ordbeta"),
      logical(1)
    )
    if (any(has_ordbeta)) {
      cli_abort(
        "The `ordbeta` family is not supported in multi-family mode: {paste(family_labels[has_ordbeta], collapse = ', ')}"
      )
    }
  }

  n_m <- max(components_per_family)
  components <- do.call(rbind, lapply(seq_len(n_f), function(i) {
    data.frame(
      family_id = i,
      component = seq_len(components_per_family[[i]]),
      family_name = family_list[[i]]$family,
      link_name = family_list[[i]]$link,
      stringsAsFactors = FALSE
    )
  }))
  components$family_code <- .map_family_codes(components$family_name, .valid_family, "family")
  components$link_code <- .map_family_codes(components$link_name, .valid_link, "link")

  combine_kind <- ifelse(
    has_two_components,
    ifelse(
      vapply(family_list, function(x) identical(x$type, "poisson_link_delta"), logical(1)),
      "poisson_link_delta",
      "delta"
    ),
    "single"
  )

  target_family <- vapply(seq_len(n_f), function(i) {
    tail(components$family_name[components$family_id == i], 1L)
  }, character(1))
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

  families <- data.frame(
    family_id = seq_len(n_f),
    label = family_labels,
    combine_kind = unname(combine_kind),
    stringsAsFactors = FALSE
  )
  spec <- list(
    version = 1L,
    n_f = n_f,
    n_m = n_m,
    family_list = family_list,
    family_labels = family_labels,
    distribution_column = if (n_f > 1L) distribution_column else NULL,
    family_id_i = as.integer(family_id_i),
    families = families,
    components = components,
    param_slot = param_slot,
    family = family_list[[1]],
    family_input = user_family
  )
  .validate_family_spec(spec)
  spec
}

.validate_family_spec <- function(family_spec) {
  required <- c("version", "n_f", "n_m", "families", "components", "family_id_i")
  if (!all(required %in% names(family_spec))) cli_abort("Internal family spec error: incomplete specification.")
  if (!family_spec$n_m %in% 1:2) cli_abort("Internal family spec error: fits require one or two components.")
  if (!identical(family_spec$families$family_id, seq_len(family_spec$n_f))) {
    cli_abort("Internal family spec error: family IDs must be consecutive and one-based.")
  }
  if (anyNA(family_spec$family_id_i) || any(!family_spec$family_id_i %in% family_spec$families$family_id)) {
    cli_abort("Internal family spec error: every row must have one valid family ID.")
  }
  counts <- tabulate(family_spec$components$family_id, nbins = family_spec$n_f)
  if (any(!counts %in% 1:2) || any(family_spec$components$component < 1L | family_spec$components$component > 2L)) {
    cli_abort("Internal family spec error: every family must have one or two active components.")
  }
  single <- family_spec$families$combine_kind == "single"
  if (any(counts[single] != 1L) || any(counts[!single] != 2L)) {
    cli_abort("Internal family spec error: combine kind does not match active components.")
  }
  invisible(family_spec)
}

.family_spec_component <- function(family_spec, family_id, component) {
  key <- family_spec$components$family_id %in% family_id &
    family_spec$components$component == component
  family_spec$components[key, , drop = FALSE]
}

.family_spec_component_value <- function(family_spec, family_id, component, column) {
  rows <- .family_spec_component(family_spec, family_id, component)
  out <- rep(NA, length(family_id))
  matched <- match(family_id, rows$family_id)
  present <- !is.na(matched)
  out[present] <- rows[[column]][matched[present]]
  out
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

  .compile_family_spec(object$family, data = object$data)
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

.family_spec_component_prediction_output <- function(x, family_spec, row_family_id,
  type = c("link", "response"), model = NA_integer_, family_list = NULL) {

  type <- match.arg(type)
  x <- as.matrix(x)
  n <- nrow(x)
  n_m <- family_spec$n_m
  if (ncol(x) < n_m) {
    cli_abort("Internal family spec error: prediction matrix has fewer components than expected.")
  }
  active <- .family_spec_component_active(family_spec, row_family_id)
  combine_kind <- family_spec$families$combine_kind[row_family_id]
  link1 <- .family_spec_component_value(family_spec, row_family_id, 1L, "link_name")
  link2 <- if (n_m > 1L) .family_spec_component_value(family_spec, row_family_id, 2L, "link_name") else rep(NA_character_, n)
  raw1 <- x[, 1L]
  raw2 <- if (n_m > 1L) x[, 2L] else rep(NA_real_, n)

  if (type == "response") {
    if (!is.null(family_list)) {
      # Use each family's own linkinv (handles families like truncated_nbinom1/2
      # whose linkinv includes a truncation correction with phi in the closure)
      est1_raw <- rep(NA_real_, n)
      est2_raw <- rep(NA_real_, n)
      for (fid in seq_len(family_spec$n_f)) {
        rows <- row_family_id == fid
        if (!any(rows)) next
        fam <- family_list[[fid]]
        fam_has_two_components <- isTRUE(fam$delta) || length(fam$family) == 2L
        linkinv1 <- if (fam_has_two_components) fam[[1]]$linkinv else fam$linkinv
        est1_raw[rows] <- linkinv1(raw1[rows])
        if (fam_has_two_components) {
          active_rows <- rows & active[, 2L]
          if (any(active_rows)) {
            est2_raw[active_rows] <- fam[[2]]$linkinv(raw2[active_rows])
          }
        }
      }
    } else {
      est1_raw <- .family_spec_apply_link(raw1, link1, inverse = TRUE)
      est2_raw <- rep(NA_real_, n)
      if (n_m > 1L && any(active[, 2L])) {
        est2_raw[active[, 2L]] <- .family_spec_apply_link(raw2[active[, 2L]], link2[active[, 2L]], inverse = TRUE)
      }
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
  } else {
    est1 <- raw1
    est2 <- rep(NA_real_, n)
    if (n_m > 1L && any(active[, 2L])) {
      est2[active[, 2L]] <- raw2[active[, 2L]]
    }
  }

  est <- if (is.na(model) || isTRUE(model == 1L)) {
    est1
  } else if (isTRUE(model == 2L)) {
    est2
  } else {
    cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
  }

  list(est = est, est1 = est1, est2 = est2)
}

# TMB simulations return realized component responses rather than prediction
# reports. Combining those draws is intentionally separate from prediction
# formatting and is limited to multiplication of realized two-part outcomes.
.family_spec_combine_simulated <- function(x, family_spec, row_family_id,
  model = NA_integer_) {

  x <- as.matrix(x)
  est1 <- x[, 1L]
  if (family_spec$n_m == 1L) {
    return(est1)
  }
  active2 <- .family_spec_component_active(family_spec, row_family_id)[, 2L]
  est2 <- x[, 2L]
  est2[!active2] <- NA_real_
  if (is.na(model)) {
    combine_kind <- family_spec$families$combine_kind[row_family_id]
    est1[combine_kind != "single"] <- est1[combine_kind != "single"] *
      est2[combine_kind != "single"]
    return(est1)
  }
  if (isTRUE(model == 1L)) return(est1)
  if (isTRUE(model == 2L)) return(est2)
  cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
}

.family_spec_prediction_link_name <- function(family_spec, row_family_id, model = NA_integer_, simulated = FALSE) {
  if (simulated) {
    return("response")
  }
  link1 <- .family_spec_component_value(family_spec, row_family_id, 1L, "link_name")
  if (family_spec$n_m == 1L) {
    links <- link1
  } else if (is.na(model)) {
    link2 <- .family_spec_component_value(family_spec, row_family_id, 2L, "link_name")
    combine_kind <- family_spec$families$combine_kind[row_family_id]
    links <- ifelse(
      combine_kind == "single",
      link1,
      ifelse(combine_kind == "poisson_link_delta", "log", link2)
    )
  } else if (isTRUE(model == 1L)) {
    links <- link1
  } else if (isTRUE(model == 2L)) {
    active2 <- .family_spec_component_active(family_spec, row_family_id)[, 2L]
    link2 <- .family_spec_component_value(family_spec, row_family_id, 2L, "link_name")
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
