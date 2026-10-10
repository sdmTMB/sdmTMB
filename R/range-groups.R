# Matern range groups --------------------------------------------------------
#
# `resolve_range_fields()` resolves which fields share a range into a table
# with one row per field and component. Everything else is derived from that
# table: the `ln_kappa` layout and map (`kappa_groups()`, `get_kappa_map()`),
# the PC prior flags (`matern_prior_flags()`), and post-fit output, which
# reads the table with `range_fields()`.

# Resolved range fields: one row per field (spatial, spatiotemporal, and each
# spatially varying coefficient, SVC) and component, with columns
# * `component`, `field` (field or SVC name), and `type` ("spatial",
#   "spatiotemporal", or "svc");
# * `active`: is the field estimated?
# * `label`: the user's `range_groups` label, or NA for a default;
# * `group`: integer range group (fields in a group share a range), or NA if
#   no active field has the field's range;
# * `range_used`: is the field's range used? As `active`, except that a
#   spatial range is also used by SVCs that share it (e.g., with
#   `spatial = "off"`);
# * `kappa_row`: the `ln_kappa` row holding the range of an active field.
# Rows are ordered spatial and spatiotemporal fields by component, then SVCs
# by component, which is the order of precedence for PC prior range parts.
# Unnamed spatiotemporal fields share their component's spatial range if
# `share_range`, and unnamed SVCs always do.
resolve_range_fields <- function(n_m, spatial, spatiotemporal, share_range,
                                 range_groups = NULL, svc = character(0),
                                 omit_spatial_intercept = FALSE) {
  user <- parse_range_groups(range_groups, n_m, svc)
  n_z <- length(svc)
  m <- c(rep(seq_len(n_m), each = 2L), rep(seq_len(n_m), each = n_z))
  type <- c(rep(c("spatial", "spatiotemporal"), n_m), rep("svc", n_z * n_m))
  f <- data.frame(component = m,
    field = c(rep(c("spatial", "spatiotemporal"), n_m), rep(svc, n_m)),
    type = type)
  f$active <- ifelse(type == "spatial",
    spatial[m] == "on" & !omit_spatial_intercept,
    type == "svc" | spatiotemporal[m] != "off")
  f$label <- vapply(seq_len(nrow(f)), function(i) {
    unname(user[[m[i]]][f$field[i]])
  }, character(1L))

  # Fields share a range if they have the same key. Prefixes keep user labels
  # from matching generated defaults.
  key <- ifelse(is.na(f$label), paste0("default:", f$field, m),
    paste0("user:", f$label))
  inherit <- is.na(f$label) &
    (type == "svc" | (type == "spatiotemporal" & share_range[m]))
  key[inherit] <- key[type == "spatial"][m[inherit]]
  f$group <- match(key, unique(key[f$active]))
  svc_on_spatial <- type == "svc" &
    (f$group == f$group[type == "spatial"][m]) %in% TRUE
  f$range_used <- f$active | (type == "spatial" & m %in% m[svc_on_spatial])

  kappa <- kappa_groups(f)
  f$kappa_row <- ifelse(type == "spatiotemporal", 2L, 1L)
  z <- which(type == "svc")
  f$kappa_row[z] <- vapply(z, function(i) {
    match(f$group[i], kappa[, m[i]])
  }, integer(1L))
  f
}

# `range_groups` as one named character vector of labels per component
parse_range_groups <- function(range_groups, n_m, svc) {
  if (is.null(range_groups)) return(rep(list(character(0)), n_m))
  fields <- c("spatial", "spatiotemporal")
  if (any(svc %in% fields)) {
    cli_abort("`range_groups` can't be used with a `spatial_varying` coefficient named `spatial` or `spatiotemporal`.")
  }
  if (!is.list(range_groups)) range_groups <- rep(list(range_groups), n_m)
  if (length(range_groups) != n_m) {
    cli_abort("`range_groups` must be a list with one element per model component ({n_m}).")
  }
  valid <- c(fields, svc)
  lapply(range_groups, function(g) {
    if (is.null(g)) return(character(0))
    if (!is.character(g) || anyNA(g) || is.null(names(g)) ||
        !all(names(g) %in% valid) || anyDuplicated(names(g))) {
      cli_abort(c("Each element of `range_groups` must be a character vector named by field.",
        "i" = "Valid names: {.val {valid}}."))
    }
    g
  })
}

# Range group held by each `ln_kappa` entry (a column per component), or NA
# if unused. This is the backend parameter layout:
# * Row 1 is the spatial range, used by the spatial field and by SVCs that
#   share it; row 2 is the spatiotemporal range.
# * If only one of rows 1 and 2 is used, the other takes its group, so the
#   `range` report is complete and `share_range` lets the backends reuse the
#   spatial precision matrix.
# * One row per SVC follows only if some SVC's range differs from row 1 of its
#   component.
kappa_groups <- function(fields) {
  groups <- field_matrix(fields$group, fields)
  used <- field_matrix(fields$range_used, fields)
  kappa <- ifelse(used[1:2, , drop = FALSE], groups[1:2, , drop = FALSE],
    NA_integer_)
  kappa[1L, ] <- ifelse(is.na(kappa[1L, ]), kappa[2L, ], kappa[1L, ])
  kappa[2L, ] <- ifelse(is.na(kappa[2L, ]), kappa[1L, ], kappa[2L, ])
  svc <- groups[-(1:2), , drop = FALSE]
  spatial <- kappa[rep(1L, nrow(svc)), , drop = FALSE]
  if (!all((svc == spatial) %in% TRUE)) kappa <- rbind(kappa, svc)
  kappa
}

# `ln_kappa` map factor from `kappa_groups()`
get_kappa_map <- function(kappa) {
  factor(match(kappa, unique(kappa[!is.na(kappa)])))
}

# Values by field as a matrix with rows spatial, spatiotemporal, then each
# SVC, and a column per component (the shape of the prior flags)
field_matrix <- function(x, fields) {
  n_m <- max(fields$component)
  svc <- fields$type == "svc"
  rbind(matrix(x[!svc], ncol = n_m), matrix(x[svc], ncol = n_m))
}

# Is a `pc_matern()` prior set? NULL for `sdmTMBpriors()` lists that predate it.
has_pc_prior <- function(prior) !is.null(prior) && !anyNA(prior[1:2])

# Which parts of the PC Matern priors apply, as matrices from `field_matrix()`:
# the sigma part for active fields, and the range part once per range group,
# from the first eligible field in `fields` order. Without `matern_svc`, the
# range part of `matern_s` also covers a spatial range used only by SVCs
# (e.g., with `spatial = "off"`).
matern_prior_flags <- function(fields, priors) {
  prior <- c(spatial = "matern_s", spatiotemporal = "matern_st",
    svc = "matern_svc")[fields$type]
  has_prior <- vapply(prior, function(p) has_pc_prior(priors[[p]]), logical(1L))
  on <- if (has_pc_prior(priors$matern_svc)) fields$active else fields$range_used
  eligible <- on & has_prior
  range_group <- ifelse(eligible, fields$group, NA)
  first <- eligible & !duplicated(range_group, incomparables = NA)
  check_range_prior_conflicts(fields, prior, range_group, priors)
  list(sigma_prior = field_matrix(fields$active & has_prior, fields),
    range_prior = field_matrix(first, fields))
}

# Warn if fields in a range group have PC Matern priors with different range
# parts (`range_gt`, `range_prob`), since only the first is applied. Different
# sigma parts aren't a conflict. `range_group` is the group of each field that
# contributes a range prior (NA otherwise).
check_range_prior_conflicts <- function(fields, prior, range_group, priors) {
  for (g in unique(range_group[!is.na(range_group)])) {
    i <- which(range_group == g)
    specs <- vapply(prior[i], function(p) {
      paste(priors[[p]][c(1L, 3L)], collapse = "/")
    }, character(1L))
    if (length(unique(specs)) == 1L) next
    who <- paste0(fields$field[i], " (`", prior[i], "`)")
    if (max(fields$component) > 1L) {
      who <- paste0("model ", fields$component[i], " ", who)
    }
    cli_warn(c(
      "Fields sharing a Mat\u00e9rn range have PC priors with different range parts (`range_gt`, `range_prob`).",
      "i" = "Fields: {who}.",
      "i" = "Only the range part from {who[1]} is applied.",
      "i" = "Use matching range settings or estimate separate ranges (`share_range` or `range_groups`)."
    ))
  }
  invisible()
}

# Resolved range fields of a fit
range_fields <- function(object) {
  if (is.null(object$range_fields)) {
    cli_abort("This model was fit with an sdmTMB version that is too old; please refit it.")
  }
  object$range_fields
}
