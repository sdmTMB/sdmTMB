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
#' @param newton_loops How many Newton optimization steps to try with
#'   [stats::optimHess()] after running [stats::nlminb()]. Sometimes aids
#'   convergence.
#' @param mgcv **Deprecated** Parse the formula with [mgcv::gam()]?
#' @param map_rf **Deprecated** use `spatial = 'off', spatiotemporal = 'off'` in
#'   [sdmTMB()].
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
#'   `lower = list(b_j = c(-5, -5))`.
#' @param upper An optional named list of upper bounds within the optimization.
#' @param get_joint_precision Logical. Passed to `getJointPrecision` in
#'   [TMB::sdreport()]. Must be `TRUE` to use simulation-based methods in
#'   [predict.sdmTMB()] or `[get_index_sims()]`. If not needed, setting this
#'   `FALSE` will reduce object size.
#' @param parallel Argument currently ignored. For parallel processing with 3
#'   cores, as an example, use `TMB::openmp(n = 3, DLL = "sdmTMB")`. But be
#'   careful, because it's not always faster with more cores and there is
#'   definitely an upper limit.
# Number of cores to use.
# Best set with the option
#   `options(sdmTMB.cores = n)`, where `n` is the number of cores.
#   Currently only works on Mac/Linux. Note that the models in sdmTMB
#   often slow down with too many cores. Ideal numbers appear to be
#   a bit less than half the available cores or ~3-4 on the machines
#   we have tested. Also propogates to prediction and index calculation.
#   Can we tweaked after the fact in a fitted model by modifying
#   `fit$control$parallel`.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @return A list of control arguments
#' @export
#' @details
#' Usually used within [sdmTMB()]. For example:
#'
#' ```
#' sdmTMB(..., control = sdmTMBcontrol(newton_loops = 1))
#' ```
#' @examples
#' sdmTMBcontrol()
sdmTMBcontrol <- function(
  eval.max = 2e3L,
  iter.max = 1e3L,
  normalize = FALSE,
  nlminb_loops = 1L,
  newton_loops = 0L,
  mgcv = deprecated(),
  quadratic_roots = FALSE,
  start = NULL,
  map_rf = deprecated(),
  map = NULL,
  lower = NULL,
  upper = NULL,
  multiphase = TRUE,
  profile = FALSE,
  get_joint_precision = TRUE,
  parallel = getOption("sdmTMB.cores", 1L),
  ...) {

  if (is_present(mgcv)) {
    deprecate_warn("0.0.20", "sdmTMBcontrol(mgcv)",
      details = "`mgcv` argument no longer does anything.")
  }

  if (is_present(map_rf)) {
    deprecate_stop("0.0.22", "sdmTMBcontrol(map_rf)", "sdmTMB(spatial = 'off', spatiotemporal = 'off')")
  }

  if (!is.null(parallel)) {
    assert_that(!is.na(parallel), parallel > 0)
    parallel <- as.integer(parallel)
  }

  out <- named_list(
    eval.max,
    iter.max,
    normalize,
    nlminb_loops,
    newton_loops,
    profile,
    quadratic_roots,
    start,
    map,
    lower,
    upper,
    multiphase,
    parallel,
    get_joint_precision
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

  n_par <- 0
  if (any(grepl("breakpt", terms_labels))) {
    out <- parse_threshold_formula(formula, "breakpt", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "linear"
    n_par <- 2
  }
  if (any(grepl("logistic", terms_labels))) {
    out <- parse_threshold_formula(formula, "logistic", terms_labels)
    threshold_parameter <- out$threshold_parameter
    formula <- out$formula
    threshold_function <- "logistic"
    n_par <- 3
  }
  if (!is.null(threshold_parameter)) {
    if (length(threshold_parameter) > 1) {
      cli_abort("`threshold_parameter` must be a single variable name.")
    }
    if (!threshold_parameter %in% c(names(data), "mi")) {
      cli_abort("`threshold_parameter` is not a column in the `data` data frame.")
    }
  }

  mi_est <- isTRUE(threshold_parameter == "mi")
  if (is.null(threshold_parameter)) {
    X_threshold <- matrix(0, nrow = nrow(data), ncol = 0L) # just placeholder
    threshold_func <- 0L
  } else {
    if (mi_est) {
      found_cols <- names(data) %in% c("po2", "invtemp")
      if (sum(found_cols) != 2)
        cli_abort("Metabolic index calcuations require columns named 'po2' and 'invtemp'.")
      X_threshold <- as.matrix(data[, c("po2", "invtemp")])
      n_par <- n_par + 1
    } else {
      X_threshold <- matrix(as.numeric(data[, names(data) == threshold_parameter, drop = TRUE], ncol = 1L))
    }
    # indexed 1, 2 because 0 will tell TMB not to estimate this:
    threshold_func <- match(threshold_function, c("linear", "logistic"))
  }

  list(
    formula = formula, threshold_parameter = threshold_parameter,
    mi_est = mi_est, n_par = n_par,
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

expand_time <- function(df, time_slices, time_column) {
  df[["weight_sdmTMB"]] <- 1
  fake_df <- df[1L, , drop = FALSE]
  fake_df[["weight_sdmTMB"]] <- 0
  missing_years <- time_slices[!time_slices %in% df[[time_column]]]
  fake_df <- do.call("rbind", replicate(length(missing_years), fake_df, simplify = FALSE))
  fake_df[[time_column]] <- missing_years
  rbind(df, fake_df)
}

safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
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

#' Check if INLA installed (i.e., not on CRAN)
#'
#' @export
#' @rdname inla_installed
#' @keywords internal
inla_installed <- function() {
  r1 <- requireNamespace("INLA", quietly = TRUE)
  r2 <- requireNamespace("rgdal", quietly = TRUE)
  r1 && r2
}

#' Check if ggplot2 installed
#'
#' @export
#' @keywords internal
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
  "omega_s" %in% names(obj$tmb_map) &&
    "epsilon_st" %in% names(obj$tmb_map) &&
    "b_rw_t" %in% names(obj$tmb_map) &&
    !"RE" %in% obj$tmb_random
}

#' Get TMB parameter list
#'
#' @param object Fit from [sdmTMB()]
#'
#' @return A named list of parameter values
#'
#' @examplesIf inla_installed()
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
