# Row-wise family mapping and prediction formatting. C++ reports are
# authoritative for combined predictions; these helpers only map rows, format
# component reports, and combine realized simulation draws.

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

.family_spec_component_prediction_output <- function(x, family_spec, row_family_id,
  type = c("link", "response"), model = NA_integer_, offset = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  n <- nrow(x)
  n_m <- family_spec$n_m
  if (ncol(x) < n_m) cli_abort("Internal family prediction error: prediction matrix has fewer components than expected.")
  active <- .family_spec_component_active(family_spec, row_family_id)
  combine_kind <- family_spec$families$combine_kind[row_family_id]
  raw1 <- x[, 1L]
  raw2 <- if (n_m > 1L) x[, 2L] else rep(NA_real_, n)
  if (type == "response") {
    est1_raw <- est2_raw <- rep(NA_real_, n)
    for (fid in seq_len(family_spec$n_f)) {
      rows <- row_family_id == fid
      if (!any(rows)) next
      fam <- family_spec$family_list[[fid]]
      has_two_components <- isTRUE(fam$delta) || length(fam$family) == 2L
      linkinv1 <- if (has_two_components) fam[[1L]]$linkinv else fam$linkinv
      est1_raw[rows] <- linkinv1(raw1[rows])
      if (has_two_components) {
        active_rows <- rows & active[, 2L]
        if (any(active_rows)) {
          est2_raw[active_rows] <- fam[[2L]]$linkinv(raw2[active_rows])
        }
      }
    }
    est1 <- est1_raw
    est2 <- est2_raw
    poisson_link_rows <- if (n_m > 1L) {
      combine_kind == "poisson_link_delta"
    } else {
      rep(FALSE, n)
    }
    if (any(poisson_link_rows)) {
      # The offset (log area swept) is carried in component 2, but encounter
      # probability depends on it too: p = 1 - exp(-a * n) and the positive
      # expectation is a * n * w / p, where est2_raw = a * w.
      if (is.null(offset)) offset <- rep(0, n)
      a <- exp(offset[poisson_link_rows])
      n_groups <- est1_raw[poisson_link_rows]
      p_encounter <- -expm1(-a * n_groups)
      est1[poisson_link_rows] <- p_encounter
      est2[poisson_link_rows] <- n_groups * est2_raw[poisson_link_rows] / p_encounter
    }
  } else {
    est1 <- raw1
    est2 <- rep(NA_real_, n)
    if (n_m > 1L && any(active[, 2L])) est2[active[, 2L]] <- raw2[active[, 2L]]
  }
  est <- if (is.na(model) || isTRUE(model == 1L)) est1 else if (isTRUE(model == 2L)) est2 else cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
  list(est = est, est1 = est1, est2 = est2)
}

# TMB simulations return realized component responses. Their combination is
# deliberately separate from prediction formatting.
.family_spec_combine_simulated <- function(x, family_spec, row_family_id, model = NA_integer_) {
  x <- as.matrix(x)
  est1 <- x[, 1L]
  if (family_spec$n_m == 1L) return(est1)
  est2 <- x[, 2L]
  est2[!.family_spec_component_active(family_spec, row_family_id)[, 2L]] <- NA_real_
  if (is.na(model)) {
    combine_kind <- family_spec$families$combine_kind[row_family_id]
    est1[combine_kind != "single"] <- est1[combine_kind != "single"] * est2[combine_kind != "single"]
    return(est1)
  }
  if (isTRUE(model == 1L)) return(est1)
  if (isTRUE(model == 2L)) return(est2)
  cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
}

.family_spec_prediction_link_name <- function(family_spec, row_family_id, model = NA_integer_, simulated = FALSE) {
  if (simulated) return("response")
  link1 <- .family_spec_component_value(family_spec, row_family_id, 1L, "link_name")
  if (family_spec$n_m == 1L || isTRUE(model == 1L)) {
    links <- link1
  } else if (is.na(model)) {
    link2 <- .family_spec_component_value(family_spec, row_family_id, 2L, "link_name")
    combine_kind <- family_spec$families$combine_kind[row_family_id]
    links <- ifelse(combine_kind == "single", link1, ifelse(combine_kind == "poisson_link_delta", "log", link2))
  } else if (isTRUE(model == 2L)) {
    active2 <- .family_spec_component_active(family_spec, row_family_id)[, 2L]
    links <- .family_spec_component_value(family_spec, row_family_id, 2L, "link_name")[active2]
  } else {
    cli_abort("`model` argument isn't valid; should be `NA`, `1`, or `2`.")
  }
  links <- unique(stats::na.omit(links))
  if (length(links) == 1L) links else "mixed"
}
