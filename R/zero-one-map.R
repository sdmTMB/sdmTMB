#' Fix factor-level intercepts for groups with all zeros or ones
#'
#' A grouping factor level (e.g. factor(year), but also possibly a region,
#' stratum, region-by-year interaction, or any other factor level) with all-zero
#' observations can drive the corresponding factor-coded coefficient to `-Inf`,
#' causing convergence problems. In delta/hurdle and binomial or beta-binomial
#' models, an all-positive or all-success level can similarly drive the
#' encounter coefficient to `+Inf`. Tweedie will have the same issue with all zeros.This helper returns `map` and `start` lists
#' suitable for passing through [sdmTMBcontrol()] to fix those problem
#' coefficients at large positive/negative values. This mirrors an approach
#' taken in \pkg{VAST}.
#'
#' Inspect the returned `all_zero_levels` and `all_one_levels` elements before
#' fitting to make sure the detected problem levels match your expectation.
#'
#' Although the motivating use case is years in fisheries index standardization
#' (e.g., `density ~ 0 + factor(year)`), the same mechanism applies to any
#' factor whose levels each get a column in the design matrix. For an
#' interaction such as region by year, create a combined column (e.g.,
#' `data$region_year <- interaction(data$region, data$year, drop = TRUE)`) and
#' use `~ 0 + region_year`.
#'
#' For Poisson-link delta families (e.g., `delta_gamma(type = "poisson-link")`),
#' the first linear predictor controls both the encounter probability and the
#' positive rate, so fixing it shifts the meaning of the second-linear-predictor
#' coefficients in those levels. The helper still works for all-zero levels (the
#' large negative value makes encounter probability ~0, and the level
#' contributes ~0 to the index), but emits a message for all-positive levels.
#' Although this affects interpretation of the coefficients, the overall
#' combined prediction remains OK.
#'
#' @param formula The fixed-effect formula to be used in [sdmTMB()]. For delta
#'   models, either a single formula applied to both linear predictors, or a
#'   `list()` of two formulas. The formula should contain `factor(group)` (or
#'   `as.factor(group)`, or a pre-converted factor column) so each level has its
#'   own column, typically `response ~ 0 + factor(year)`.
#' @param data The data frame passed to [sdmTMB()].
#' @param group Character. Name of the grouping column in `data` whose levels
#'   should be inspected for all-zero/all-positive observations. Most commonly
#'   the time column.
#' @param family The `sdmTMB` family object (e.g., `tweedie()`,
#'   `delta_gamma()`).
#' @param weights Optional numeric vector of trial sizes for non-delta binomial
#'   or beta-binomial models fit with proportion data. This should match the
#'   `weights` argument passed to [sdmTMB()]. Ignored for other families.
#' @param value Positive numeric. Magnitude at which to fix the intercept;
#'   `plogis(20)`, `plogis(-20)`, and `exp(-20)` are all within ~2e-9 of 0 or 1,
#'   which is plenty.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`map`}{A list with `b_j` (and `b_j2` for delta families) factors
#'       suitable for `sdmTMBcontrol(map = ...)`.}
#'     \item{`start`}{A list with `b_j` (and `b_j2`) starting-value vectors
#'       suitable for `sdmTMBcontrol(start = ...)`.}
#'     \item{`all_zero_levels`}{Group levels with all-zero observations.}
#'     \item{`all_one_levels`}{Group levels with all-positive observations
#'       for delta families, or all-success observations for non-delta
#'       binomial/beta-binomial families. Empty for other families.}
#'   }
#'
#' @examples
#' d <- pcod
#' d$density[d$year == 2003] <- 0 # force an all-zero year
#'
#' m <- make_zero_one_map(
#'   formula = density ~ 0 + factor(year),
#'   data = d,
#'   group = "year",
#'   family = tweedie()
#' )
#' m$all_zero_levels
#' m$all_one_levels
#' m$start
#' m$map
#'
#' \donttest{
#' mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#' fit <- sdmTMB(
#'   density ~ 0 + factor(year),
#'   data = d,
#'   mesh = mesh,
#'   time = "year",
#'   family = tweedie(),
#'   spatial = "on",
#'   spatiotemporal = "off", # off for example speed
#'   control = sdmTMBcontrol(map = m$map, start = m$start) #<
#' )
#' fit
#' coef(fit)
#' }
#' @seealso [sdmTMBcontrol()], [sdmTMB()]
#' @export
make_zero_one_map <- function(formula, data, group, family,
                               weights = NULL, value = 20) {
  assert_that(inherits(data, "data.frame"))
  assert_that(is.character(group), length(group) == 1L)
  if (!group %in% names(data)) {
    cli_abort("`group` column {.val {group}} not found in `data`.")
  }
  assert_that(is.numeric(value), length(value) == 1L, value > 0)
  if (!is.null(weights)) {
    assert_that(is.numeric(weights), length(weights) == nrow(data))
  }

  if (!inherits(family, "family") && !inherits(family, "sdmTMBfamily")) {
    cli_abort("`family` must be an sdmTMB family object.")
  }
  delta <- isTRUE(family$delta)
  binomial_like <- identical(family$family[1], "binomial") ||
    identical(family$family[1], "betabinomial")

  # Normalize formula to list of length 1 (non-delta) or 2 (delta)
  formulas <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(formulas) || !all(vapply(formulas, inherits, logical(1), "formula"))) {
    cli_abort("`formula` must be a formula or a list of formulas.")
  }
  if (delta && length(formulas) == 1L) formulas <- list(formulas[[1]], formulas[[1]])
  if (!delta && length(formulas) != 1L) {
    cli_abort("Non-delta family requires a single formula.")
  }
  if (delta && length(formulas) != 2L) {
    cli_abort("Delta family requires a single formula or a list of two formulas.")
  }

  f1 <- formulas[[1]]
  if (length(f1) < 3L) {
    cli_abort("The first formula must be two-sided so the response can be inferred.")
  }
  response <- as.character(f1[[2]])
  if (length(response) != 1L) {
    cli_abort("Could not infer the response from the first formula.")
  }
  if (!response %in% names(data)) {
    cli_abort("`response` column {.val {response}} not found in `data`.")
  }

  to_native <- function(x, gvec) {
    if (is.integer(gvec)) suppressWarnings(as.integer(x))
    else if (is.numeric(gvec)) suppressWarnings(as.numeric(x))
    else x
  }

  build_model_frame <- function(form) {
    f_nobars <- reformulas::nobars(form)
    f_clean <- remove_s_and_t2(f_nobars)
    if (length(f_clean) < 3L) {
      f_clean <- stats::as.formula(
        paste(response, paste(deparse(f_clean), collapse = " "))
      )
    }
    mf <- stats::model.frame(f_clean, data = data)
    na_action <- attr(mf, "na.action")
    weights_sub <- .subset_by_na_action(weights, na_action)
    rows <- .get_model_frame_rows(mf, data)
    list(
      formula = f_clean,
      mf = mf,
      weights = weights_sub,
      group = data[[group]][rows]
    )
  }

  detect_group_columns <- function(X, form, gvec) {
    tt <- stats::terms(form)
    term_labels <- gsub(" ", "", attr(tt, "term.labels"))
    supported_terms <- gsub(" ", "", c(
      group,
      paste0("factor(", group, ")"),
      paste0("as.factor(", group, ")")
    ))
    term_match <- which(term_labels %in% supported_terms)
    if (length(term_match) != 1L) {
      cli_abort(c(
        "Could not identify a matching factor-coded term for {.val {group}} in the design matrix.",
        "i" = "Use something like {.code ~ 0 + factor(group)} or a pre-converted factor column with {.code ~ 0 + group}."
      ))
    }

    if (term_labels[term_match] == group && !is.factor(gvec)) {
      cli_abort(c(
        "The bare {.val {group}} term must already be a factor column.",
        "i" = "Use {.code factor(group)} in the formula or convert the column in `data` first."
      ))
    }

    gfac <- if (is.factor(gvec)) gvec else factor(gvec)
    G <- stats::model.matrix(~ 0 + gfac)
    cols <- which(attr(X, "assign") == term_match)
    if (length(cols) != ncol(G) ||
        !isTRUE(all.equal(
          as.matrix(X[, cols, drop = FALSE]),
          as.matrix(G),
          check.attributes = FALSE
        ))) {
      cli_abort(c(
        "Could not match {.val {group}} levels to one-hot columns in the design matrix.",
        "i" = "Make sure the formula includes {.code 0 + factor(group)} or {.code 0 + group} with a factor column."
      ))
    }

    list(columns = cols, levels = levels(gfac))
  }

  frame1 <- build_model_frame(formulas[[1]])
  gvec1 <- frame1$group
  if (binomial_like && !delta) {
    encounter <- .process_binomial_response(
      frame1$mf,
      weights = frame1$weights,
      weights_arg = "`weights`"
    )
    y <- encounter$y_i
    size <- encounter$size
  } else {
    y <- stats::model.response(frame1$mf, type = "numeric")
    size <- rep(1, length(y))
  }

  g_split <- split(seq_along(y), gvec1)
  all_zero <- names(g_split)[vapply(g_split, function(idx) {
    length(idx) > 0L && all(y[idx] == 0)
  }, logical(1))]

  if (delta) {
    all_one <- names(g_split)[vapply(g_split, function(idx) {
      length(idx) > 0L && all(y[idx] > 0)
    }, logical(1))]
  } else if (binomial_like) {
    all_one <- names(g_split)[vapply(g_split, function(idx) {
      length(idx) > 0L && all(y[idx] == size[idx])
    }, logical(1))]
  } else {
    all_one <- character(0)
  }

  all_zero_native <- to_native(all_zero, gvec1)
  all_one_native <- to_native(all_one, gvec1)

  # Poisson-link warning
  is_poisson_link <- isTRUE(family$type == "poisson_link_delta")
  if (is_poisson_link && length(all_one) > 0L) {
    cli::cli_inform(c(
      "Detected all-positive {group} level(s) with a Poisson-link delta family.",
      "i" = "Fixing the first linear predictor for those levels will shift the meaning of the second-linear-predictor coefficients. The combined prediction is still OK. See {.code ?make_zero_one_map}."
    ))
  }

  if (length(all_zero) == 0L && length(all_one) == 0L) {
    cli::cli_inform("No all-zero or all-positive {group} levels detected; returning empty `map`/`start`.")
  }

  # Build map/start per linear predictor by inspecting design-matrix columns
  # for the group term on the same NA-filtered rows used by model.matrix().
  build_for <- function(form, fix_zero, fix_one) {
    frame <- build_model_frame(form)
    gvec <- frame$group
    X <- stats::model.matrix(frame$formula, frame$mf)
    group_cols <- detect_group_columns(X, frame$formula, gvec)

    map_vec <- seq_len(ncol(X))
    start_vec <- rep(0, ncol(X))

    locate_col <- function(lvl) {
      idx <- match(to_native(lvl, gvec), to_native(group_cols$levels, gvec))
      if (is.na(idx)) {
        cli_abort(c(
          "Could not find a column in the design matrix corresponding to {group} = {.val {lvl}}.",
          "i" = "Make sure the formula uses the same grouped factor as the fitted model."
        ))
      }
      group_cols$columns[idx]
    }

    for (lvl in fix_zero) {
      j <- locate_col(lvl)
      map_vec[j] <- NA
      start_vec[j] <- -value
    }
    for (lvl in fix_one) {
      j <- locate_col(lvl)
      map_vec[j] <- NA
      start_vec[j] <- +value
    }
    list(map = factor(map_vec), start = start_vec)
  }

  m1 <- build_for(formulas[[1]], fix_zero = all_zero, fix_one = all_one)

  map_out <- list(b_j = m1$map)
  start_out <- list(b_j = m1$start)

  if (delta) {
    # Model 2 (positive component): all-zero levels have no positive data to
    # inform b_j2; mapping b_j2 off there prevents drift. All-one levels
    # contain only positive obs, so b_j2 is free there.
    m2 <- build_for(formulas[[2]], fix_zero = all_zero, fix_one = character(0))
    map_out$b_j2 <- m2$map
    start_out$b_j2 <- m2$start
  }

  list(
    map = map_out,
    start = start_out,
    all_zero_levels = all_zero_native,
    all_one_levels = all_one_native
  )
}
