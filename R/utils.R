#' Optimization control options
#'
#' [sdmTMB()] and [stats::nlminb()] control options.
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param nlminb_loops How many times to run [stats::nlminb()] optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#' @param newton_loops How many Newton optimization steps to try after running
#'   [stats::nlminb()]. This sometimes aids convergence by further reducing the
#'   log-likelihood gradient with respect to the fixed effects. This calculates
#'   the Hessian at the current MLE with [stats::optimHess()] using a
#'   finite-difference approach and uses this to update the fixed effect
#'   estimates.
#' @param map A named list with factor `NA`s specifying parameter values that
#'   should be fixed at a constant value. See the documentation in
#'   [TMB::MakeADFun()]. This should usually be used with `start` to specify the
#'   fixed value.
#' @param start A named list specifying the starting values for parameters. You
#'   can see the necessary structure by fitting the model once and inspecting
#'   `your_model$tmb_obj$env$parList()`. Elements of `start` that are specified
#'   will replace the default starting values.
#' @param quadratic_roots Experimental feature for internal use right now; may
#'   be moved to a branch. Logical: should quadratic roots be calculated? Note:
#'   on the sdmTMB side, the first two coefficients are used to generate the
#'   quadratic parameters. This means that if you want to generate a quadratic
#'   profile for depth, and depth and depth^2 are part of your formula, you need
#'   to make sure these are listed first and that an intercept isn't included.
#'   For example, `formula = cpue ~ 0 + depth + depth2 + as.factor(year)`.
#' @param getsd Logical indicating whether to call [TMB::sdreport()].
#' @param normalize Logical: use [TMB::normalize()] to normalize the process
#'   likelihood using the Laplace approximation? Can result in a substantial
#'   speed boost in some cases. This used to default to `FALSE` prior to
#'   May 2021. Currently not working for models fit with REML or random intercepts.
#' @param multiphase Logical: estimate the fixed and random effects in phases?
#'   Phases are usually faster and more stable.
#' @param profile Logical: should population-level/fixed effects be profiled
#'   out of the likelihood? These are then appended to the random effects
#'   vector without the Laplace approximation. See [TMB::MakeADFun()]. *This
#'   can dramatically speed up model fit if there are many fixed effects but is
#'   experimental at this stage.*
#' @param lower An optional named list of lower bounds within the optimization.
#'   Parameter vectors with the same name (e.g., `b_j` or `ln_kappa` in some
#'   cases) can be specified as a numeric vector. E.g.
#'   `lower = list(b_j = c(-5, -5))`. Note that [stats::optimHess()] does not
#'   implement lower and upper bounds, so you must set `newton_loops = 0` if
#'   setting limits.
#' @param upper An optional named list of upper bounds within the optimization.
#' @param censored_upper An optional vector of upper bounds for
#'   [sdmTMBcontrol()]. Values of `NA` indicate an unbounded right-censored to
#'   distribution, values greater that the observation indicate and upper bound,
#'   and values equal to the observation indicate no censoring.
#' @param get_joint_precision Logical. Passed to `getJointPrecision` in
#'   [TMB::sdreport()]. Must be `TRUE` to use simulation-based methods in
#'   [predict.sdmTMB()] or `[get_index_sims()]`. If not needed, setting this
#'   `FALSE` will reduce object size.
#' @param parallel Argument currently ignored. For parallel processing with 3
#'   cores, as an example, use `TMB::openmp(n = 3, DLL = "sdmTMB")`. But be
#'   careful, because it's not always faster with more cores and there is
#'   definitely an upper limit.
#' @param suppress_nlminb_warnings Suppress uninformative warnings
#'   from [stats::nlminb()] arising when a function evaluation is `NA`, which
#'   are then replaced with `Inf` and avoided during estimation?
#' @param collapse_spatial_variance Logical: should spatial and/or spatiotemporal
#'   random fields be automatically dropped if their estimated standard deviation
#'   is effectively zero (i.e., below `collapse_threshold`)? This helps prevent
#'   overfitting and numerical instability when the data provide little evidence
#'   for spatial or spatiotemporal variation. I.e., when the variance parameter is
#'   estimated on or near the boundary of zero. When enabled, the model will be
#'   automatically refitted via [update.sdmTMB()] with the corresponding field(s)
#'   disabled. This adds a computational cost (a single model refit if
#'   collapsing occurs) but can yield a simpler, more stable model and more
#'   reliable inference. Default is `FALSE` for backwards compatibility.
#' @param collapse_threshold Numeric: the standard deviation threshold below which random
#'   fields are considered to be collapsing to zero. Only used when
#'   `collapse_spatial_variance = TRUE`. Values are on the standard deviation
#'   scale (i.e., square root of variance). Default is 0.01.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @return A list of control arguments
#' @export
#' @details
#' Usually used within [sdmTMB()]. For example:
#'
#' ```
#' sdmTMB(..., control = sdmTMBcontrol(newton_loops = 2))
#' ```
#' @examples
#' sdmTMBcontrol()
sdmTMBcontrol <- function(
  eval.max = 2e3L,
  iter.max = 1e3L,
  normalize = FALSE,
  nlminb_loops = 1L,
  newton_loops = 1L,
  getsd = TRUE,
  quadratic_roots = FALSE,
  start = NULL,
  map = NULL,
  lower = NULL,
  upper = NULL,
  censored_upper = NULL,
  multiphase = TRUE,
  profile = FALSE,
  get_joint_precision = TRUE,
  parallel = getOption("sdmTMB.cores", 1L),
  suppress_nlminb_warnings = TRUE,
  collapse_spatial_variance = FALSE,
  collapse_threshold = 0.01,
  ...) {

  assert_that(is.numeric(nlminb_loops), is.numeric(newton_loops))
  assert_that(nlminb_loops >= 1L)
  assert_that(newton_loops >= 0L)
  # if (newton_loops > 1L) {
  #   cli::cli_inform("There is rarely a benefit to making `newton_loops` > 1.")
  # }

  if (!is.null(parallel)) {
    assert_that(!is.na(parallel), parallel > 0)
    parallel <- as.integer(parallel)
  }

  assert_that(is.logical(profile) || is.character(profile))
  assert_that(is.logical(collapse_spatial_variance))
  assert_that(is.numeric(collapse_threshold), collapse_threshold > 0)

  out <- named_list(
    eval.max,
    iter.max,
    normalize,
    nlminb_loops,
    newton_loops,
    getsd,
    profile,
    quadratic_roots,
    start,
    map,
    lower,
    upper,
    censored_upper,
    multiphase,
    parallel,
    get_joint_precision,
    collapse_spatial_variance,
    collapse_threshold
  )
  c(out, list(...))
}

set_par_value <- function(opt, par) {
  as.numeric(opt$par[par == names(opt$par)])
}

get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  pdHess <- isTRUE(sd_report$pdHess)
  invisible(named_list(final_grads, bad_eig, pdHess))
}

make_year_i <- function(x) {
  x <- as.integer(as.factor(x))
  x - min(x)
}

make_time_lu <- function(time_vec_from_data, full_time_vec = sort(unique(time_vec_from_data))) {
  if (!all(time_vec_from_data %in% full_time_vec)) {
    stop("All time elements not in full time vector.")
  }
  lu <- unique(
    data.frame(
      year_i = make_year_i(full_time_vec),
      time_from_data = full_time_vec,
      stringsAsFactors = FALSE
    )
  )
  lu$extra_time <- !lu$time_from_data %in% time_vec_from_data
  lu <- lu[order(lu$time_from_data),]
  row.names(lu) <- NULL
  lu
}

check_offset <- function(formula) {
  .check <- any(grepl("^offset$",
    gsub(" ", "", unlist(strsplit(as.character(formula), "\\+")))))
  if (.check)
    cli_abort("Contains offset in formula. This is deprecated. Please use the `offset` argument.")
}

check_and_parse_thresh_params <- function(formula, data) {
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  if (any(grepl("linear_thresh", terms_labels)) && any(grepl("logistic_thresh", terms_labels))) {
    cli_abort("Please include only a linear (`breakpt`) *or* a logistic threshold.")
  }
  if (sum(grepl("linear_thresh", terms_labels)) > 1 || sum(grepl("logistic_thresh", terms_labels)) > 1) {
    cli_abort("Please include only a *single* threshold variable.")
  }
  threshold_parameter <- NULL
  if (any(grepl("breakpt", terms_labels))) {
    out <- parse_threshold_formula(formula, "breakpt", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "linear"
  }
  if (any(grepl("logistic", terms_labels))) {
    out <- parse_threshold_formula(formula, "logistic", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "logistic"
  }
  if (!is.null(threshold_parameter)) {
    if (length(threshold_parameter) > 1) {
      cli_abort("`threshold_parameter` must be a single variable name.")
    }
    if (!threshold_parameter %in% names(data)) {
      cli_abort("`threshold_parameter` is not a column in the `data` data frame.")
    }
  }

  if (is.null(threshold_parameter)) {
    X_threshold <- rep(0, nrow(data)) # just placeholder
    threshold_func <- 0L
  } else {
    X_threshold <- data[, names(data) == threshold_parameter, drop = TRUE]
    # indexed 1, 2 because 0 will tell TMB not to estimate this:
    threshold_func <- match(threshold_function, c("linear", "logistic"))
  }
  X_threshold <- as.numeric(unlist(X_threshold))
  list(
    formula = formula, threshold_parameter = threshold_parameter,
    threshold_func = threshold_func, X_threshold = X_threshold
  )
}

parse_threshold_formula <- function(formula, thresh_type_short = "lin_thresh",
  terms_labels) {
  which_thresh <- grep(thresh_type_short, terms_labels)
  temp <- gsub(paste0("^", thresh_type_short, "\\("), "", terms_labels[which_thresh])
  threshold_parameter <- gsub("\\)$", "", temp)
  formula <- stats::update(formula, paste("~ . -", terms_labels[which_thresh]))
  list(formula = formula, threshold_parameter = threshold_parameter)
}

expand_time <- function(df, time_slices, time_column, weights, offset, upr) {
  df[["__weight_sdmTMB__"]] <- if (!is.null(weights)) weights else  1
  df[["__sdmTMB_offset__"]] <- if (!is.null(offset)) offset else 0
  df[["__dcens_upr__"]] <- if (!is.null(upr)) upr else NA_real_
  df[["__fake_data__"]] <- FALSE
  fake_df <- df[1L, , drop = FALSE]
  fake_df[["__fake_data__"]] <- TRUE
  fake_df[["__weight_sdmTMB__"]] <- 0 # IMPORTANT: this turns off these data in the likelihood
  missing_years <- time_slices[!time_slices %in% df[[time_column]]]
  fake_df <- do.call("rbind", replicate(length(missing_years), fake_df, simplify = FALSE))
  fake_df[[time_column]] <- missing_years
  rbind(df, fake_df)
}

safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}

# Extract a short name from a call element without expensive deparse
extract_call_name <- function(call_element, max_width = 80) {
  if (is.null(call_element)) {
    return(NULL)
  }

  # Fast path for symbols (most common case: mesh_obj)
  if (is.symbol(call_element)) {
    return(as.character(call_element))
  }

  # Fast path for function calls (e.g., make_mesh(...))
  if (is.call(call_element)) {
    # Use deparse1 (available in R >= 4.0.0, package requires >= 4.1.0)
    deparsed <- deparse1(call_element, width.cutoff = max_width)
    # If still too long after width cutoff, truncate and add ellipsis
    if (nchar(deparsed) > max_width) {
      deparsed <- paste0(substr(deparsed, 1, max_width - 3), "...")
    }
    return(deparsed)
  }

  obj_class <- class(call_element)[1]
  return(obj_class)
}

barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

check_valid_factor_levels <- function(x, .name = "") {
  assert_that(is.factor(x),
    msg = sprintf("Random effect group column `%s` is not a factor.", .name))
  lev <- sort(levels(x))
  uni <- sort(unique(as.character(x)))
  assert_that(identical(lev, uni),
    msg = sprintf(
      "Random effect group column `%s` has extra factor levels. Please remove them.", .name))
}

#' Check if ggplot2 installed
#'
#' @export
#' @keywords internal
#' @return Returns `TRUE` or `FALSE`.
ggplot2_installed <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

remove_s_and_t2 <- function(formula) {
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  drop <- grep("s\\(", terms_labels)
  dropt2 <- grep("t2\\(", terms_labels)
  tdrop <- terms_labels[union(drop, dropt2)]
  for (i in seq_along(tdrop)) {
    formula <- stats::update(formula, paste("~ . -", tdrop[i]))
  }
  formula
}

has_no_random_effects <- function(obj) {
  length(obj$tmb_random) == 0L
}

#' Get TMB parameter list
#'
#' @param object Fit from [sdmTMB()]
#'
#' @return A named list of parameter values
#'
#' @examples
#' fit <- sdmTMB(present ~ 1, data = pcod_2011, family = binomial(), spatial = "off")
#' pars <- get_pars(fit)
#' names(pars)
#' @export
get_pars <- function(object) {
  # based on glmmTMB:
  ee <- object$tmb_obj$env
  x <- ee$last.par.best
  if (length(ee$random) > 0) x <- x[-ee$random]
  p <- ee$parList(x = x)
  p
}

#' Replicate a prediction data frame over time
#'
#' Useful for replicating prediction grids across time slices used in model
#' fitting.
#'
#' @param dat Data frame.
#' @param time_name Name of time column in output.
#' @param time_values Time values to replicate `dat` over.
#'
#' @return
#' A data frame replicated over `time_values` with a new column based on
#' `time_name`.
#'
#' @examples
#'
#' df <- data.frame(variable = c("a", "b"))
#' replicate_df(df, time_name = "year", time_values = 1:3)
#'
#' head(qcs_grid)
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#' head(nd)
#' table(nd$year)
#' @export
replicate_df <- function(dat, time_name, time_values) {
  nd <- do.call("rbind", replicate(length(time_values), dat, simplify = FALSE))
  nd[[time_name]] <- rep(time_values, each = nrow(dat))
  nd
}

# work one delta model at a time
# only diff with 2nd part is that the 'fake' values need to be larger than first
# so make a function for just 1 component
# then add an increment for 2nd that's always bigger (e.g. 100s vs. 1000s)
# then as.factor() them down in sequence
# check if share_range = TRUE or if one of spatial or spatiotemporal is 'off',
# if so the values in that column should be identical
# check if share_range = TRUE or if one of spatial or spatiotemporal is 'off',
# if so the values in that column should be identical
map_kappa <- function(spatial, spatiotemporal, share_range, a = 100L) {
  if (share_range) {
    if (!spatial && !spatiotemporal) {
      x <- c(NA_integer_, NA_integer_)
    } else {
      x <- c(a, a)
    }
  } else {
    if (spatial && spatiotemporal) x <- c(a, a + 1L)
    if (!spatial || !spatiotemporal) x <- c(a, a)
    if (!spatial && !spatiotemporal) x <- c(NA_integer_, NA_integer_)
  }
  x
}

get_kappa_map <- function(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)) {
  spatial <- spatial == "on"
  spatiotemporal <- spatiotemporal %in% c("on", "iid", "rw", "ar1")
  k <- map_kappa(spatial[1], spatiotemporal[1], share_range[1], 100L)
  if (n_m > 1) {
    k2 <- map_kappa(spatial[2], spatiotemporal[2], share_range[2], 1000L)
    k <- cbind(k, k2)
  }
  as.factor(as.integer(as.factor(k)))
}


#' Calculate scaling factor for hook competition using censored method
#'
#' @param prop_removed A vector (numeric) containing the proportion of baits
#'   removed in each fishing event.
#' @param n_hooks A vector (integers) containing the number of hooks
#'   deployed in each fishing event.
#' @param pstar A single value (numeric) specifying the breakdown point of
#'   observed catch counts as a result of hook competition.
#' @return A vector containing the scale factor used to calculate the
#'   upper bound on catch counts of target species to improve convergence of
#'   censored method. See \doi{10.1139/cjfas-2022-0159}, including supplementary
#'   materials section S.3 for more details.
#' @noRd
get_scale_factor <- function(prop_removed, n_hooks, pstar) {
  # Apply pstar correction
  prop_hook <- signif(((prop_removed - pstar) / (1 - pstar)), 5)
  n_hooks <- round((1 - pstar) * n_hooks)
  # Calculate competition adjustment factor
  prop <- 1 - prop_hook
  prop[prop == 0] <- 1 / n_hooks[prop == 0] # if all hooks saturated - map to 1 hook
  -log(prop) / (1 - prop)
}

#' Calculate an upper bound on catch counts for the censored Poisson family
#'
#' @param prop_removed The proportion of baits removed in each fishing event
#'   from *any* species. I.e., the proportion of hooks returning without bait
#'   for any reason.
#' @param n_catch The observed catch counts on each fishing event of the target
#'   species.
#' @param n_hooks The number of hooks deployed on each fishing event.
#' @param pstar A single value between `0 <= pstar <= 1` specifying the
#'   breakdown point of observed catch counts as a result of hook competition.
#'
#' @details `pstar` could be obtained via inspecting a GAM or other smoother fit
#'   with catch counts as the response, an offset for log(hook count), and
#'   proportion of baits removed for each fishing event as the predictor.
#'   Check when the curve drops off as the proportion bait removed increases.
#'
#' The `lwr` limit for [sdmTMB::censored_poisson()] should be the observed catch
#' counts, i.e., `n_catch` here.
#'
#' If `upr` in [sdmTMB::censored_poisson()] is set to NA, the full
#' right-censored Poisson likelihood is used without any upper bound.
#'
#' The right-censored Poisson density can be written as:
#'
#' ```
#'   dcens_pois <- function(x, lambda) {
#'       1 - ppois(x - 1, lambda)
#'    }
#' ```
#'
#' and the right-censored Poisson density with an upper limit can be written as:
#'
#' ```
#'   dcens_pois_upper <- function(x, lambda, upper) {
#'     ppois(upper, lambda) - ppois(x - 1, lambda)
#'   }
#' ```
#'
#' In practice, these computations are done in log space for numerical
#' stability.
#'
#' @return A numeric vector of upper bound catch counts of the target species to
#'   improve convergence of the censored method.
#'
#' @references See \doi{10.1139/cjfas-2022-0159} for more details.
#' @noRd
#'
#' @examples
#' dat <- structure(
#'   list(
#'     n_catch = c(
#'       78L, 63L, 15L, 6L, 7L, 11L, 37L, 99L, 34L, 100L, 77L, 79L,
#'       98L, 30L, 49L, 33L, 6L, 28L, 99L, 33L
#'     ),
#'     prop_removed = c(
#'       0.61, 0.81, 0.96, 0.69, 0.99, 0.98, 0.25, 0.95, 0.89, 1, 0.95, 0.95,
#'       0.94, 1, 0.95, 1, 0.84, 0.3, 1, 0.99
#'     ), n_hooks = c(
#'       140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L,
#'       140L, 140L, 140L, 140L, 140L, 140L, 140L, 140L
#'     )
#'   ),
#'   class = "data.frame", row.names = c(NA, -20L)
#' )
#' upr <- get_censored_upper(dat$prop_removed, dat$n_catch, dat$n_hooks, pstar = 0.9)
#' upr
#'
#' plot(dat$n_catch, upr, xlab = "N catch", ylab = "N catch upper limit")
#' abline(0, 1, lty = 2)
#' above_pstar <- dat[dat$prop_removed > 0.9, ]
#' upr_pstar <- upr[dat$prop_removed > 0.9]
#' points(above_pstar$n_catch, upr_pstar, col = "red", pch = 20)
#' txt <- paste0(
#'   "Red indicates catch events\n",
#'   "with >= pstar proportion of hooks\n",
#'   "coming back without bait.\n\n",
#'   "The rest are fishing events < pstar\n",
#'   "('high-quality' events)."
#' )
#' text(10, 120, txt, adj = 0)
get_censored_upper <- function(
    prop_removed,
  n_catch,
  n_hooks,
  pstar = 0.95) {
  assertthat::assert_that(
    is.numeric(prop_removed),
    is.numeric(n_catch),
    is.numeric(n_hooks)
  )
  assertthat::assert_that(length(prop_removed) == length(n_catch))
  assertthat::assert_that(length(prop_removed) == length(n_hooks))
  assertthat::assert_that(pstar >= 0)
  assertthat::assert_that(pstar <= 1)
  assertthat::assert_that(all(prop_removed <= 1))
  assertthat::assert_that(all(prop_removed >= 0))
  assertthat::assert_that(all(n_catch >= 0))
  assertthat::assert_that(all(n_hooks >= 0))
  assertthat::assert_that(sum(is.na(c(prop_removed, n_catch, n_hooks))) == 0)

  removed_ind <- prop_removed > pstar
  # FIXME: add cli::cli_abort()?
  # probably need to throw error if length(removed_ind) == 0L
  upper_bound <- rep(0, length(prop_removed))
  scale_fac <- rep(0, length(prop_removed))

  # Calculate part of scaling factor used to get upper bound (part of eq. S.20)
  scale_fac[removed_ind] <-
    get_scale_factor(
      pstar = pstar,
      prop_removed = prop_removed[removed_ind],
      n_hooks = n_hooks[removed_ind]
    )

  # Upper bound for a species (part of eq. S.22)
  upper_bound[removed_ind] <- (prop_removed[removed_ind] - pstar) *
    n_hooks[removed_ind] * scale_fac[removed_ind]

  high <- n_catch
  high[prop_removed >= pstar] <- high[prop_removed >= pstar] +
    upper_bound[prop_removed >= pstar]
  round(high)
}

#' Set delta model for [ggeffects::ggpredict()]
#'
#' Set a delta model component to predict from with [ggeffects::ggpredict()].
#'
#' @param x An [sdmTMB::sdmTMB()] model fit with a delta family such as
#'   [sdmTMB::delta_gamma()].
#' @param model Which delta/hurdle linear predictor to predict/plot with.
#'   `NA` does the combined prediction, `1` does the binomial part, and `2`
#'   does the positive part.
#'
#' @details
#' A complete version of the examples below would be:
#'
#' ```
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2), data = pcod_2011,
#'   spatial = "off", family = delta_gamma())
#'
#' # binomial part:
#' set_delta_model(fit, model = 1) |>
#'   ggeffects::ggpredict("depth_scaled [all]")
#'
#' # gamma part:
#' set_delta_model(fit, model = 2) |>
#'   ggeffects::ggpredict("depth_scaled [all]")
#'
#' # combined:
#' set_delta_model(fit, model = NA) |>
#'   ggeffects::ggpredict("depth_scaled [all]")
#' ```
#'
#' But cannot be run on CRAN until the next version of \pkg{ggeffects} is available
#' on CRAN. For now, you can install the GitHub version of \pkg{ggeffects}.
#' <https://github.com/strengejacke/ggeffects>.
#'
#' @returns
#' The fitted model with a new attribute named `delta_model_predict`.
#' We suggest you use `set_delta_model()` in a pipe (as in the examples)
#' so that this attribute does not persist. Otherwise, [predict.sdmTMB()]
#' will choose this model component by default. You can also remove the
#' attribute yourself after:
#'
#' ```
#' attr(fit, "delta_model_predict") <- NULL
#' ```
#'
#' @examplesIf require("ggeffects", quietly = TRUE)
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2), data = pcod_2011,
#'   spatial = "off", family = delta_gamma())
#'
#' # binomial part:
#' set_delta_model(fit, model = 1)
#'
#' # gamma part:
#' set_delta_model(fit, model = 2)
#'
#' # combined:
#' set_delta_model(fit, model = NA)
#' @export
set_delta_model <- function(x, model = c(NA, 1, 2)) {
  assertthat::assert_that(model[[1]] %in% c(NA, 1, 2),
    msg = "`model` argument not valid; should be one of NA, 1, 2")
  attr(x, "delta_model_predict") <- model[[1]]
  x
}

get_fitted_time <- function(x) {
  if (!"fitted_time" %in% names(x))
    cli_abort("Missing 'fitted_time' element in fitted object. Please refit the model with a current version of sdmTMB.")
  x$fitted_time
}

reload_model <- function(object) {
  if ("parlist" %in% names(object)) {
    # tinyVAST does this to be extra sure... I've found one case where it was needed
    obj <- TMB::MakeADFun(
      data = object$tmb_data,
      parameters = object$parlist, #!! important part
      map = object$tmb_map,
      random = object$tmb_random,
      DLL = "sdmTMB",
      profile = object$control$profile
    )
    obj$env$beSilent()
    nll_new <- obj$fn(object$model$par) #!! important: need to eval once (restores last.par.best etc.)
    if (abs(nll_new - object$model$objective) > 0.01) {
      cli_abort(c("Model fit is not identical to recorded value:", "
        something is not working as expected"))
    }
    object$tmb_obj <- obj
    object
  } else {
    cli_abort("`reload_model()` only works with models fit with sdmTMB 0.5.0.9006 and higher.")
  }
}

reinitialize <- function(x) {
  # replacement for TMB:::isNullPointer; modified from glmmTMB source
  # https://github.com/glmmTMB/glmmTMB/issues/651#issuecomment-912920255
  # https://github.com/glmmTMB/glmmTMB/issues/651#issuecomment-914542795
  is_null_pointer <- function(x) {
    x <- x$tmb_obj$env$ADFun$ptr
    attributes(x) <- NULL
    identical(x, new("externalptr"))
  }
  if (is_null_pointer(x)) {
    x$tmb_obj$env$beSilent()
    x$tmb_obj$fn(x$model$par)
    # x$tmb_obj$retape()
    if ("parlist" %in% names(x) && "last.par.best" %in% names(x)) {
      if (!identical(x$tmb_obj$env$last.par.best, x$last.par.best)) {
        cli_warn(c("Detected a potential issue reloading a saved sdmTMB model.",
          "Please run `fit <- sdmTMB:::reload_model(fit)`,",
          "where `fit` is your fitted model."))
      }
    }
  }
}

#' Set up spatially varying coefficients for category composition models
#'
#' This function helps set up the data structure, formula, and mapping needed
#' for fitting spatially varying coefficient models with categories (e.g., ages,
#' length bins, species) that have both spatial and spatiotemporal random
#' fields. It's particularly useful for age or length composition
#' standardization models.
#'
#' @param data Data frame containing the composition data.
#' @param category_column Character. Name of the category column (e.g., "Age",
#'   "length_bin", "species").
#' @param time_column Character. Name of the time column (e.g., "Year").
#' @param share_spatial_sd Logical. If `TRUE`, all categories share the same
#'   spatial SD. If `FALSE`, each category gets its own spatial SD.
#' @param share_spatiotemporal_sd Logical. If `TRUE`, all category-time
#'   combinations share the same spatiotemporal SD. If `FALSE`, each gets its
#'   own.
#'
#' @return A list containing:
#' \itemize{
#'   \item `data_expanded`: Data frame with added model matrix columns for use in [sdmTMB()].
#'   \item `svc_formula`: Formula for the `spatial_varying` argument in [sdmTMB()].
#'   \item `svc_map`: Map list for the `map` argument in [sdmTMBcontrol()].
#'   \item `info`: List with summary information about the model structure.
#' }
#'
#' @details This function creates spatially varying coefficient structures for
#' composition models by setting up:
#'
#' 1. **Spatial fields**: One field per category (e.g., age-specific spatial
#' fields) 2. **Spatiotemporal fields**: One field per category-time combination
#' (e.g., age-year fields)
#'
#' The sharing of variance parameters is controlled by `share_spatial_sd` and
#' `share_spatiotemporal_sd`. When `TRUE`, all fields of that type share the
#' same variance parameter, which is more parsimonious but assumes similar
#' variance magnitudes across categories.
#'
#' The resulting model structure allows each category to have its own spatial
#' pattern and temporal variation while controlling parameter sharing for
#' identifiability and computational efficiency.
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   age = factor(rep(1:3, each = 20)),
#'   year = rep(2020:2022, 20),
#'   abundance = rnorm(60),
#'   x = runif(60), y = runif(60)
#' )
#'
#' # Set up model components
#' setup <- make_category_svc(
#'   data = data,
#'   category_column = "age",
#'   time_column = "year",
#'   share_spatial_sd = TRUE,
#'   share_spatiotemporal_sd = TRUE
#' )
#'
#' # Check the setup
#' setup$info
#'
#' # See the age composition standardization vignette for more details
#' @export
make_category_svc <- function(data,
                               category_column,
                               time_column,
                               share_spatial_sd = TRUE,
                               share_spatiotemporal_sd = TRUE) {
  # Input validation
  assert_that(is.data.frame(data))
  assert_that(is.character(category_column), length(category_column) == 1L)
  assert_that(is.character(time_column), length(time_column) == 1L)
  assert_that(is.logical(share_spatial_sd), length(share_spatial_sd) == 1L)
  assert_that(is.logical(share_spatiotemporal_sd), length(share_spatiotemporal_sd) == 1L)

  if (!category_column %in% names(data)) {
    cli_abort(paste0("Column '", category_column, "' not found in data."))
  }
  if (!time_column %in% names(data)) {
    cli_abort(paste0("Column '", time_column, "' not found in data."))
  }

  # Create model matrices
  spatial_formula <- as.formula(paste0("~ 0 + ", category_column))
  spatiotemporal_formula <- as.formula(paste0("~ 0 + factor(", time_column, "):", category_column))

  spatial_terms <- model.matrix(spatial_formula, data = data)
  spatiotemporal_terms <- model.matrix(spatiotemporal_formula, data = data)

  # Combine matrices
  combined_matrix <- cbind(spatial_terms, spatiotemporal_terms)
  colnames_combined <- colnames(combined_matrix)

  # Create mapping based on sharing preferences
  n_spatial <- ncol(spatial_terms)
  n_spatiotemporal <- ncol(spatiotemporal_terms)

  if (share_spatial_sd && share_spatiotemporal_sd) {
    # Two shared SDs: one for spatial, one for spatiotemporal
    svc_map <- factor(c(rep(1, n_spatial), rep(2, n_spatiotemporal)))
  } else if (share_spatial_sd && !share_spatiotemporal_sd) {
    # Shared spatial SD, separate spatiotemporal SDs
    svc_map <- factor(c(rep(1, n_spatial), 2:(n_spatiotemporal + 1)))
  } else if (!share_spatial_sd && share_spatiotemporal_sd) {
    # Separate spatial SDs, shared spatiotemporal SD
    svc_map <- factor(c(1:n_spatial, rep(n_spatial + 1, n_spatiotemporal)))
  } else {
    # All separate SDs
    svc_map <- factor(1:(n_spatial + n_spatiotemporal))
  }

  # Expand data
  data_expanded <- cbind(data, combined_matrix)

  # Create formula
  svc_formula <- as.formula(paste0("~ `", paste(colnames_combined, collapse = "` + `"), "`"))

  # Create info summary
  info <- list(
    n_categories = length(unique(data[[category_column]])),
    n_times = length(unique(data[[time_column]])),
    n_spatial_terms = n_spatial,
    n_spatiotemporal_terms = n_spatiotemporal,
    n_variance_parameters = length(unique(svc_map)),
    spatial_terms = colnames(spatial_terms),
    spatiotemporal_terms = colnames(spatiotemporal_terms),
    share_spatial_sd = share_spatial_sd,
    share_spatiotemporal_sd = share_spatiotemporal_sd
  )

  return(list(
    data_expanded = data_expanded,
    svc_formula = svc_formula,
    svc_map = list(ln_tau_Z = svc_map),
    info = info
  ))
}
