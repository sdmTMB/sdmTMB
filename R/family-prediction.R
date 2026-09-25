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

.family_spec_inverse_link <- function(eta, link) {
  switch(link,
    identity = eta, log = exp(eta), logit = stats::plogis(eta), inverse = 1 / eta,
    cloglog = 1 - exp(-exp(eta)),
    cli_abort("Unsupported link in family prediction: {.val {link}}")
  )
}

.family_spec_link <- function(mu, link) {
  switch(link,
    identity = mu, log = log(mu), logit = stats::qlogis(mu), inverse = 1 / mu,
    cloglog = log(-log1p(-mu)),
    cli_abort("Unsupported link in family prediction: {.val {link}}")
  )
}

.family_spec_apply_link <- function(x, link, inverse = TRUE) {
  if (!length(x)) return(x)
  out <- x
  for (this_link in unique(link[!is.na(link)])) {
    ii <- link == this_link
    out[ii] <- if (inverse) .family_spec_inverse_link(x[ii], this_link) else .family_spec_link(x[ii], this_link)
  }
  out
}

.family_spec_component_prediction_output <- function(x, family_spec, row_family_id,
  type = c("link", "response"), model = NA_integer_, family_list = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  n <- nrow(x)
  n_m <- family_spec$n_m
  if (ncol(x) < n_m) cli_abort("Internal family prediction error: prediction matrix has fewer components than expected.")
  active <- .family_spec_component_active(family_spec, row_family_id)
  combine_kind <- family_spec$families$combine_kind[row_family_id]
  link1 <- .family_spec_component_value(family_spec, row_family_id, 1L, "link_name")
  link2 <- if (n_m > 1L) .family_spec_component_value(family_spec, row_family_id, 2L, "link_name") else rep(NA_character_, n)
  raw1 <- x[, 1L]
  raw2 <- if (n_m > 1L) x[, 2L] else rep(NA_real_, n)
  if (type == "response") {
    if (!is.null(family_list)) {
      est1_raw <- est2_raw <- rep(NA_real_, n)
      for (fid in seq_len(family_spec$n_f)) {
        rows <- row_family_id == fid
        if (!any(rows)) next
        fam <- family_list[[fid]]
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
    } else {
      est1_raw <- .family_spec_apply_link(raw1, link1)
      est2_raw <- rep(NA_real_, n)
      if (n_m > 1L && any(active[, 2L])) est2_raw[active[, 2L]] <- .family_spec_apply_link(raw2[active[, 2L]], link2[active[, 2L]])
    }
    est1 <- est1_raw
    est2 <- est2_raw
    poisson_link_rows <- if (n_m > 1L) {
      combine_kind == "poisson_link_delta"
    } else {
      rep(FALSE, n)
    }
    if (any(poisson_link_rows)) {
      n_groups <- est1_raw[poisson_link_rows]
      p_encounter <- 1 - exp(-n_groups)
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
