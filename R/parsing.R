# barnames is an internal function to return the group (bar) names in a formula,
# e.g. (1|group)
# @param bars The character string representing the random effects, e.g. `1 |
#   group`
barnames <- function(bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
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
  fe_form <- reformulas::nobars(f) # fixed effect formula, no bars
  re_cov_terms <- NULL

  re_cov_terms <- list(
    Zt = NULL, theta = NULL, Lind = NULL, Gp = NULL,
    lower = NULL, Lambdat = NULL, flist = NULL,
    cnms = NULL, Ztlist = NULL, nl = NULL
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
    level_ids = integer(0),
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
      reorder.vars = FALSE
    ) # keep not alphabetical
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
      levels <- levels(data[, bn[[i]]])
      # Add the group and ":" to the name for each -- otherwise the level names might be the
      # same, especially if the levels ids are integers for 2+ groups
      level_ids <- paste0(bn[[i]], ":", rownames(re_cov_terms$Ztlist[[i]]))
      if (i == 1) {
        df <- data.frame(
          index = seq_len(length(level_ids)),
          level_ids = level_ids,
          group_indices = i
        )
      } else {
        df <- rbind(df, data.frame(
          index = seq_len(length(level_ids)),
          level_ids = level_ids,
          group_indices = i
        ))
      }
    }
    df$index <- seq(1, nrow(df))
    re_cov_terms$re_b_df <- data.frame(
      level_ids = unique(df$level_ids),
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
    re_cov_terms$re_b_df <- re_cov_terms$re_b_df[, c("level_ids", "start", "end", "group_indices")]
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
