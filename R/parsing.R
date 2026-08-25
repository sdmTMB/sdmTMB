# barnames is an internal function to return the group (bar) names in a formula,
# e.g. (1|group)
# @param bars The character string representing the random effects, e.g. `1 |
#   group`
barnames <- function(bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

# termnames is an internal function to return the names of the intercept and/or
# slope terms before a bar in a formula, e.g. (1|group), (1+slope|group),
# or (0+slope|group)
# @param bars The character string representing the random effects, e.g. `1 |
#   group`
termnames <- function(bars) {
  vapply(bars, function(x) safe_deparse(x[[2]]), "")
}

# make_indices is an internal function to build lower triangular matrices for
# correlated random effects
# @param vec Vector of indices to generate row and col positions for
make_indices <- function(vec) {
  group_indices <- unlist(lapply(seq_along(vec), function(i) rep(i, sum(1:vec[i]))))

  # Function to generate lower triangular indices
  get_lower_tri_indices <- function(n) {
    cols <- unlist(lapply(seq_len(n), function(i) rep(i, n - i + 1)))
    rows <- unlist(lapply(seq_len(n), function(i) seq(i, n)))
    list(rows = rows, cols = cols)
  }

  # Getting row and column indices for each group
  col_indices <- unlist(lapply(vec, function(x) get_lower_tri_indices(x)$cols))
  row_indices <- unlist(lapply(vec, function(x) get_lower_tri_indices(x)$rows))
  # return list of indices, all starting at 0 for TMB
  list(group_indices = group_indices - 1L, rows = row_indices - 1L, cols = col_indices - 1L)
}

# parse_formula is an internal function to parse random effects in a formula and
# return objects for estimation
# @param f formula object
# @param data data frame used to build the random effects
parse_formula <- function(f, data) {
  b <- reformulas::findbars(f) # find expressions separated by |, NULL if no RE
  bn <- barnames(b) # names of groups
  tn <- termnames(b) # names of terms
  fe_form <- reformulas::nobars(f) # fixed effect formula, no bars
  re_cov_terms <- NULL

  re_cov_terms <- list(
    Zt = NULL, theta = NULL,
    Lind = NULL, Gp = NULL,
    lower = NULL, Lambdat = NULL,
    flist = NULL, cnms = NULL,
    Ztlist = NULL, nl = NULL
  )
  re_cov_terms$re_df <- data.frame(
    group_indices = integer(0),
    rows = integer(0), cols = integer(0),
    is_sd = integer(0), par_index = integer(0)
  )
  re_cov_terms$re_cov_term_map <- data.frame(
    group = integer(0),
    dim = integer(0),
    start = integer(0),
    end = integer(0)
  )
  re_cov_terms$re_b_df <- data.frame(
    level_ids = character(0),
    level = character(0),
    start = integer(0), # index of beta vec in TMB
    end = integer(0), # index of beta vec in TMB
    group_indices = integer(0), # which group are these levels associated with
    var_start = integer(0), # index of variances to use for these betas and groups
    var_end = integer(0)
  ) # index of variances to use for these betas and groups
  re_cov_terms$re_b_map <- data.frame(
    group = integer(0),
    dim = integer(0),
    start = integer(0),
    end = integer(0)
  )
  var_indx_vector <- 0

  if (length(bn) > 0) {
    mf <- model.frame(reformulas::subbars(f), data)
    re_cov_terms <- reformulas::mkReTrms(b, mf,
      drop.unused.levels = TRUE,
      reorder.terms = FALSE, # default is true, reorder based on dec levels
      reorder.vars = FALSE, # keep not alphabetical
      calc.lambdat = FALSE # Lambdat and Lind needed for lme4, but not glmmTMB (according to Bolker),
      # so probably not needed for sdmTMB
    )
    # re_cov_terms$theta gives the total number of params across v-cov matrices
    # see lme4 vignettes for construction details. These are indexes with
    # Lind which maps elements of theta to the VCov matrices.

    # this function creates replicated indices per element
    group_dims <- unlist(lapply(re_cov_terms$cnms, length)) # dimensions of RE for each group
    re_cov_terms$re_df <- as.data.frame(make_indices(group_dims))
    re_cov_terms$re_df$is_sd <- ifelse(
      re_cov_terms$re_df$rows == re_cov_terms$re_df$cols, 1L, 0L
    ) # used for TMB transform
    # index, e.g. 0 - 15 to be used for TMB indexing:
    re_cov_terms$re_df$par_index <- seq_len(nrow(re_cov_terms$re_df)) - 1L

    # also need to map these indices to the vector of estimated covariance parameters
    re_cov_terms$re_cov_term_map <- data.frame(
      group = unique(re_cov_terms$re_df$group_indices),
      dim = group_dims,
      start = NA, end = NA
    )
    for (i in seq_len(nrow(re_cov_terms$re_cov_term_map))) {
      indx <- which(re_cov_terms$re_df$group_indices == re_cov_terms$re_cov_term_map$group[i])
      re_cov_terms$re_cov_term_map$start[i] <- min(indx) - 1L # start at 0
      re_cov_terms$re_cov_term_map$end[i] <- max(indx) - 1L # start at 0
    }

    # index the level / group of the elements of Zt
    for (i in seq_len(length(re_cov_terms$Ztlist))) {
      # Add the term(s) and group to the name for each (separated by ":") -- otherwise
      # the level names might be the same, especially if the levels ids are integers for 2+ groups.
      levels <- rownames(re_cov_terms$Ztlist[[i]])
      level_ids <- paste0(tn[[i]], ":", bn[[i]], ":", levels)
      if (i == 1) {
        df <- data.frame(
          index = seq_len(length(level_ids)),
          level_ids = level_ids,
          level = levels,
          group_indices = i
        )
      } else {
        df <- rbind(df, data.frame(
          index = seq_len(length(level_ids)),
          level_ids = level_ids,
          level = levels,
          group_indices = i
        ))
      }
    }
    df$index <- seq(1, nrow(df))
    re_cov_terms$re_b_df <- data.frame(
      level_ids = unique(df$level_ids),
      level = df$level[!duplicated(df$level_ids)],
      start = NA, end = NA, group_indices = NA
    )
    for (i in seq_len(nrow(re_cov_terms$re_b_df))) {
      indx <- which(df$level_ids == re_cov_terms$re_b_df$level_ids[i])
      re_cov_terms$re_b_df$start[i] <- min(indx) - 1L # start at 0
      re_cov_terms$re_b_df$end[i] <- max(indx) - 1L # start at 0
      re_cov_terms$re_b_df$group_indices[i] <- df$group_indices[indx[1]]
    }
    # add the variance index -- largely for groups with 1 type of RE
    re_cov_terms$re_b_df$var_start <- 0
    re_cov_terms$re_b_df$var_end <- 0
    group_index <- 0

    for (i in seq_len(nrow(re_cov_terms$re_b_df))) {
      if (re_cov_terms$re_b_df$group_indices[i] > group_index) {
        if (i == 1) {
          start_index <- re_cov_terms$re_b_df$start[i]
          end_index <- re_cov_terms$re_b_df$end[i]
        } else {
          start_index <- end_index + 1L
          end_index <- start_index + (re_cov_terms$re_b_df$end[i] -
            re_cov_terms$re_b_df$start[i])
        }
        group_index <- group_index + 1L
      }
      re_cov_terms$re_b_df$var_start[i] <- start_index
      re_cov_terms$re_b_df$var_end[i] <- end_index
    }
    var_indx_vector <- unlist(
      mapply(seq,
        from = re_cov_terms$re_b_df$var_start,
        to = re_cov_terms$re_b_df$var_end, SIMPLIFY = FALSE
      )
    )
    re_cov_terms$re_b_df <- re_cov_terms$re_b_df[, c("level_ids", "level", "start", "end", "group_indices")]
    # index the groups
    # also need to map these indices to the vector of estimated covariance parameters
    re_cov_terms$re_b_map <- data.frame(
      group = unique(re_cov_terms$re_b_df$group_indices),
      start = NA, end = NA
    )
    for (i in seq_len(nrow(re_cov_terms$re_b_map))) {
      indx <- which(re_cov_terms$re_b_df$group_indices == re_cov_terms$re_b_map$group[i])
      re_cov_terms$re_b_map$start[i] <- min(indx) - 1L # start at 0
      re_cov_terms$re_b_map$end[i] <- max(indx) - 1L # start at 0
    }
  }
  list(
    bars = b, barnames = bn, form_no_bars = fe_form, n_bars = length(bn),
    re_cov_terms = re_cov_terms, var_indx_vector = var_indx_vector
  )
}

add_model_index <- function(split_formula, dataframe_name) {
  lapply(seq_along(split_formula), function(i) {
    # Access the specific data frame using the provided data frame name
    df <- split_formula[[i]]$re_cov_terms[[dataframe_name]]
    nrows <- nrow(df)
    if (nrows > 0) {
      df$model <- i # Add new column with the position in the list
    } else {
      df$model <- integer(0)
    }
    # Assign the modified data frame back to the original list structure
    split_formula[[i]]$re_cov_terms[[dataframe_name]] <- df

    df
  })
}

#' Mark a model-derived time index in a spatially varying formula
#'
#' `model_index()` marks one complete fixed-effect term whose model-derived,
#' centered time index will be used as a single spatially varying covariate.
#' It is only valid inside the `spatial_varying` argument to [sdmTMB()].
#' A term with multiple model-matrix columns, such as a factor, still creates
#' exactly one spatially varying coefficient field.
#'
#' @param x One complete term that also occurs in the main model formula.
#'
#' @return This function is a formula marker and is not evaluated directly.
#'
#' @details
#' For model component \eqn{m} and fitted time \eqn{t}, `model_index()` first
#' constructs
#' \deqn{a_{tm} = X^I_{tm} b_m + \bar{\epsilon}_{tm},}
#' where \eqn{X^I} contains the fixed-effect columns belonging to `x` plus the
#' fixed-effect intercept when present. \eqn{\bar{\epsilon}_{tm}} is the
#' arithmetic mean of the raw spatiotemporal field across mesh vertices. The
#' resulting expected response is placed on the log scale and centered across
#' fitted times. For a standard delta model, component expectations are
#' multiplied before taking the log and centering; Poisson-link delta models
#' use the sum of their two component predictors.
#'
#' The centered index is jointly estimated with the model and interacts with
#' one existing SVC field per model component. It does not replace the selected
#' fixed effect or the local spatiotemporal field. Because the spatial mean is
#' an equal-vertex mean, its value can depend on mesh construction.
#'
#' Exactly one `model_index()` term is supported and `time` is required. The
#' wrapped expression must exactly match a complete fixed-effect term in every
#' component, and its design must be constant within time. `extra_time` and
#' [project()] are unsupported; [predict()] and [simulate.sdmTMB()] support
#' fitted time values.
#'
#' The resulting SVC describes an endogenous association between the fitted
#' model-derived index and spatial redistribution. It is not, by itself,
#' evidence of a causal density-dependent mechanism.
#'
#' @references
#' Thorson, J. T. (2022). Development and simulation testing for a new
#' approach to density dependence in species distribution models. *ICES
#' Journal of Marine Science*, 79, 117--128. \doi{10.1093/icesjms/fsab247}.
#'
#' @examples
#' \dontrun{
#' pcod$year_factor <- factor(pcod$year)
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
#' fit <- sdmTMB(
#'   density ~ 0 + year_factor,
#'   spatial_varying = ~ 0 + model_index(year_factor),
#'   data = pcod,
#'   mesh = mesh,
#'   time = "year",
#'   family = tweedie(link = "log")
#' )
#' }
#' @export
model_index <- function(x) {
  cli_abort("`model_index()` can only be used inside the `spatial_varying` argument to `sdmTMB()`.")
}

.parse_model_index_svc <- function(spatial_varying) {
  inactive <- list(
    ordinary_formula = spatial_varying,
    active = FALSE,
    term = NULL,
    label = NULL
  )
  if (is.null(spatial_varying)) return(inactive)
  if (!inherits(spatial_varying, "formula") || length(spatial_varying) != 2L) {
    cli_abort("`spatial_varying` must be a one-sided formula when using `model_index()`.")
  }

  tt <- stats::terms(spatial_varying, specials = "model_index", keep.order = TRUE)
  labels <- attr(tt, "term.labels")
  n_specials <- length(attr(tt, "specials")$model_index)
  special_labels <- labels[grepl("(^|:)model_index\\(", labels)]
  if (!n_specials) return(inactive)
  if (n_specials > 1L || length(special_labels) > 1L) {
    cli_abort("At most one `model_index()` term may be used in `spatial_varying`.")
  }

  label <- special_labels[[1L]]
  call <- tryCatch(str2lang(label), error = function(e) NULL)
  direct <- is.call(call) && identical(call[[1L]], as.name("model_index"))
  if (!direct || length(call) != 2L || grepl("model_index", safe_deparse(call[[2L]]), fixed = TRUE)) {
    cli_abort("Use `model_index()` as one direct term with exactly one argument in `spatial_varying`.")
  }
  term <- safe_deparse(call[[2L]])
  other_labels <- labels[labels != label]
  ordinary_formula <- stats::reformulate(
    other_labels,
    intercept = attr(tt, "intercept") == 1L,
    env = environment(spatial_varying)
  )

  list(
    ordinary_formula = ordinary_formula,
    active = TRUE,
    term = term,
    label = label
  )
}

.model_index_coefficient_indices <- function(model_index, fixed_terms, X_ij) {
  if (!isTRUE(model_index$active)) return(NULL)
  lapply(seq_along(X_ij), function(m) {
    labels <- attr(fixed_terms[[m]], "term.labels")
    term_number <- match(model_index$term, labels)
    if (is.na(term_number)) {
      cli_abort("The term `{model_index$term}` inside `model_index()` must occur as a complete term in every model formula.")
    }
    assign <- attr(X_ij[[m]], "assign")
    which(assign == 0L | assign == term_number)
  })
}

.build_model_index_X_t <- function(X_ij, coefficient_indices, year_i, n_t, term) {
  lapply(seq_along(X_ij), function(m) {
    out <- matrix(0, nrow = n_t, ncol = ncol(X_ij[[m]]))
    cols <- coefficient_indices[[m]]
    for (t in seq_len(n_t) - 1L) {
      rows <- which(year_i == t)
      values <- X_ij[[m]][rows, cols, drop = FALSE]
      if (!length(rows) || any(apply(values, 2L, function(x) any(x != x[[1L]])))) {
        cli_abort("The fixed-effect design for `model_index({term})` must be constant within each fitted time value.")
      }
      out[t + 1L, cols] <- values[1L, , drop = TRUE]
    }
    colnames(out) <- colnames(X_ij[[m]])
    out
  })
}
