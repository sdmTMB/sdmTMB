ll_gaussian <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  stats::dnorm(x = withheld_y, mean = withheld_mu, sd = .sd, log = TRUE)
}

ll_tweedie <- function(object, withheld_y, withheld_mu) {
  p <- stats::plogis(object$model$par[["thetaf"]]) + 1
  phi <- exp(object$model$par[["ln_phi"]])
  fishMod::dTweedie(y = withheld_y, mu = withheld_mu, p = p, phi = phi, LOG = TRUE)
}

ll_binomial <- function(object, withheld_y, withheld_mu) {
  stats::dbinom(x = withheld_y, size = 1, prob = withheld_mu, log = TRUE)
}

ll_gamma <- function(object, withheld_y, withheld_mu) {
  .shape <- exp(object$model$par[["ln_phi"]])
  stats::dgamma(x = withheld_y, shape = .shape, scale = withheld_mu / .shape, log = TRUE)
}

ll_lognormal <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  stats::dlnorm(x = withheld_y, meanlog = withheld_mu - 0.5 * (.sd)^2, sdlog = .sd, log = TRUE)
}

dstudent <- function(x, df, mean, sd, ncp, log = FALSE) {
  # from metRology::dt.scaled()
  if (!log) {
    return(stats::dt((x - mean) / sd, df, ncp = ncp, log = FALSE) / sd)
  } else {
    return(stats::dt((x - mean) / sd, df, ncp = ncp, log = TRUE) - log(sd))
  }
}

ll_student <- function(object, withheld_y, withheld_mu) {
  .sd <- exp(object$model$par[["ln_phi"]])
  dstudent(x = withheld_y, df = object$tmb_data$df, mean = withheld_mu, sd = .sd, log = TRUE)
}

ll_nbinom1 <- function(object, withheld_y, withheld_mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  stats::dnbinom(x = withheld_y, size = withheld_mu / phi, mu = withheld_mu, log = TRUE)
}

ll_nbinom2 <- function(object, withheld_y, withheld_mu) {
  phi <- exp(object$model$par[["ln_phi"]])
  stats::dnbinom(x = withheld_y, size = phi, mu = withheld_mu, log = TRUE)
}

# no longer used within sdmTMB_cv(); uses TMB report() instead
ll_sdmTMB <- function(object, withheld_y, withheld_mu) {
  family_func <- switch(object$family$family,
    gaussian = ll_gaussian,
    tweedie = ll_tweedie,
    binomial = ll_binomial,
    lognormal = ll_lognormal,
    student = ll_student,
    Gamma = ll_gamma,
    nbinom1 = ll_nbinom1,
    nbinom2 = ll_nbinom2,
    cli_abort(paste0(
      object$family$family, " not yet implemented. ",
      "Please file an issue on GitHub."
    ))
  )
  family_func(object, withheld_y, withheld_mu)
}

#' Cross validation with sdmTMB models
#'
#' Performs k-fold or leave-future-out cross validation with sdmTMB models.
#' Returns the sum of log likelihoods of held-out data (log predictive density),
#' which can be used to compare models—higher values indicate better
#' out-of-sample prediction. By default, creates folds randomly and stratified
#' by time (set a seed for reproducibility), but folds can be manually assigned
#' via `fold_ids`. See Ward and Anderson (2025) in the References and the
#' [cross-validation vignette](https://sdmTMB.github.io/sdmTMB/articles/cross-validation.html).
#'
#' @param formula Model formula.
#' @param data A data frame.
#' @param mesh Output from [make_mesh()]. If supplied, the same mesh will be
#'   used for all folds. This is faster and usually what you want.
#' @param mesh_args Arguments for [make_mesh()]. If supplied, the mesh will be
#'   reconstructed for each fold.
#' @param time The name of the time column. Leave as `NULL` if this is only
#'   spatial data.
#' @param k_folds Number of folds.
#' @param fold_ids Optional vector containing user fold IDs. Can also be a
#'   single string, e.g. `"fold_id"` representing the name of the variable in
#'   `data`. Ignored if `lfo` is TRUE
#' @param lfo Logical. Use leave-future-out (LFO) cross validation? If `TRUE`,
#'   data from earlier time steps are used to predict future time steps. The
#'   `time` argument must be specified. See Details section below.
#' @param lfo_forecast If `lfo = TRUE`, number of time steps ahead to forecast.
#'   For example, `lfo_forecast = 1` means fitting to time steps 1 to T and
#'   validating on T + 1. See Details section below.
#' @param lfo_validations If `lfo = TRUE`, number of times to step through the
#'   LFO process (i.e., number of validation folds). Defaults to 5. See Details
#'   section below.
#' @param parallel If `TRUE` and a [future::plan()] is supplied, will be run in
#'   parallel.
#' @param use_initial_fit Fit the first fold and use those parameter values
#'   as starting values for subsequent folds? Can be faster with many folds.
#' @param save_models Logical. If `TRUE` (default), the fitted model object for
#'   each fold is stored in the output. If `FALSE`, models are not saved, which
#'   can substantially reduce memory usage for large datasets or many folds.
#'   When `FALSE`, functions that require access to the fitted models (e.g.,
#'   [tidy()], [cv_to_waywiser()]) will not work.
#' @param future_globals A character vector of global variables used within
#'   arguments if an error is returned that \pkg{future.apply} can't find an
#'   object. This vector is appended to `TRUE` and passed to the argument
#'   `future.globals` in [future.apply::future_lapply()]. Useful if global
#'   objects are used to specify arguments like priors, families, etc.
#' @param ... All other arguments required to run the [sdmTMB()] model. The
#'   `weights` argument is supported and will be combined with the internal
#'   fold-assignment mechanism (held-out data are assigned weight 0).
#'
#' @export
#' @return
#' A list:
#' * `data`: Original data plus columns for fold ID (`cv_fold`), CV predicted
#'   value (`cv_predicted`), CV log likelihood (`cv_loglik`), and CV deviance
#'   residuals (`cv_deviance_resid`).
#' * `models`: A list of fitted models, one per fold. `NULL` if `save_models = FALSE`.
#' * `fold_loglik`: Sum of log likelihoods of held-out data per fold (log
#'   predictive density per fold). More positive values indicate better
#'   out-of-sample prediction.
#' * `sum_loglik`: Sum of `fold_loglik` across all folds (total log predictive
#'   density). Use this to compare models—more positive values are better.
#' * `pdHess`: Logical vector: was the Hessian positive definite for each fold?
#' * `converged`: Logical: did all folds converge (all `pdHess` `TRUE`)?
#' * `max_gradients`: Maximum absolute gradient for each fold.
#'
#' @details
#' **Parallel processing**
#'
#' Parallel processing can be used by setting a `future::plan()`.
#'
#' For example:
#'
#' ```
#' library(future)
#' plan(multisession)
#' # now use sdmTMB_cv() ...
#' ```
#'
#' **Leave-future-out cross validation (LFOCV)**
#'
#' An example of LFOCV with 9 time steps, `lfo_forecast = 1`, and
#' `lfo_validations = 2`:
#'
#' - Fit data to time steps 1 to 7, predict and validate step 8.
#' - Fit data to time steps 1 to 8, predict and validate step 9.
#'
#' An example of LFOCV with 9 time steps, `lfo_forecast = 2`, and
#' `lfo_validations = 3`:
#'
#' - Fit data to time steps 1 to 5, predict and validate step 7.
#' - Fit data to time steps 1 to 6, predict and validate step 8.
#' - Fit data to time steps 1 to 7, predict and validate step 9.
#'
#' Note these are time steps as they are presented in order in the data.
#' For example, in the `pcod` data example below steps between data points
#' are not always one year but an `lfo_forecast = 2` forecasts 2 time
#' steps as presented not two years.
#'
#' See example below.
#'
#' @references
#'
#' Ward, E.J., and S.C. Anderson. 2025. Approximating spatial processes with
#' too many knots degrades the quality of probabilistic predictions.
#' bioRxiv 2025.11.14.688354. \doi{10.1101/2025.11.14.688354}.
#'
#' @examples
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
#'
#' # Set parallel processing first if desired with the future package.
#' # See the Details section above.
#'
#' m_cv <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh = mesh, spatial = "off",
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' m_cv$fold_loglik
#' m_cv$sum_loglik
#'
#' head(m_cv$data)
#' m_cv$models[[1]]
#' m_cv$max_gradients
#'
#'
#' \donttest{
#' # Create mesh each fold:
#' m_cv2 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh_args = list(xy_cols = c("X", "Y"), cutoff = 20),
#'   family = tweedie(link = "log"), k_folds = 2
#' )
#'
#' # Use fold_ids:
#' m_cv3 <- sdmTMB_cv(
#'   density ~ 0 + depth_scaled + depth_scaled2,
#'   data = pcod, mesh = mesh,
#'   family = tweedie(link = "log"),
#'   fold_ids = rep(seq(1, 3), nrow(pcod))[seq(1, nrow(pcod))]
#' )
#'
#' # LFOCV:
#' m_lfocv <- sdmTMB_cv(
#'   present ~ s(year, k = 4),
#'   data = pcod,
#'   lfo = TRUE,
#'   lfo_forecast = 2,
#'   lfo_validations = 3,
#'   family = binomial(),
#'   mesh = mesh,
#'   spatial = "off", # fast example
#'   spatiotemporal = "off", # fast example
#'   time = "year" # must be specified
#' )
#'
#' # See how the LFOCV folds were assigned:
#' fold_table <- table(m_lfocv$data$cv_fold, m_lfocv$data$year)
#' fold_table
#' }
sdmTMB_cv <- function(
    formula, data, mesh_args, mesh = NULL, time = NULL,
    k_folds = 8, fold_ids = NULL,
    lfo = FALSE,
    lfo_forecast = 1,
    lfo_validations = 5,
    parallel = TRUE,
    use_initial_fit = FALSE,
    save_models = TRUE,
    future_globals = NULL,
    ...) {
  if (k_folds < 1) cli_abort("`k_folds` must be >= 1.")

  spde <- mesh
  data[["_sdm_order_"]] <- seq_len(nrow(data))
  constant_mesh <- missing(mesh_args)
  if (missing(mesh_args)) mesh_args <- NULL
  if (missing(spde)) spde <- NULL
  if (lfo) fold_ids <- NULL
  # add column of fold_ids stratified across time steps
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }
  if (is.null(fold_ids)) {
    if (lfo) {
      if (length(unique(data[[time]])) < (lfo_validations + lfo_forecast)) {
        cli_abort("Not enough time steps for the desired validation period. Either decrease `lfo_validations` or add more data.")
      }
      # Create lfo_validations + 1 folds, ordered sequentially
      data$cv_fold <- 1
      t_validate <- sort(unique(data[[time]]), decreasing = TRUE)
      for (t in seq(1, lfo_validations+lfo_forecast)) {
        # fold id increasing order + forecast
        data$cv_fold[data[[time]] == t_validate[t]] <- lfo_validations - t + 1 + lfo_forecast
      }
    } else {
      dd <- lapply(split(data, data[[time]]), function(x) {
        x$cv_fold <- sample(rep(seq(1L, k_folds), length.out = nrow(x)), size = nrow(x))
        x
      })
      data <- do.call(rbind, dd)
    }
    fold_ids <- "cv_fold"
  } else {
    # fold_ids passed in; can be numeric, or a named column in `data`
    data$cv_fold <- NA
    if (length(fold_ids) == nrow(data)) {
      data$cv_fold <- fold_ids
    }
    if (length(fold_ids) == 1L && is.character(fold_ids)) {
      if (!fold_ids %in% names(data)) {
        cli_abort("Name of fold identifier not found in data.")
      }
      data$cv_fold <- data[[fold_ids]]
    }
    if (length(fold_ids) > 1 && length(fold_ids) < nrow(data)) {
      cli_abort("Dimension of `fold_ids` doesn't match data and is not a named variable.")
    }
    if (length(which(is.na(data$cv_fold))) > 0) {
      cli_abort("NAs found in `fold_ids`; please check `fold_ids` are specified correctly.")
    }
    k_folds <- length(unique(data$cv_fold))
  }
  if (time == "_sdmTMB_time") { # undo changes above, make time NULL
    data[["_sdmTMB_time"]] <- NULL
    time <- NULL
  }

  dot_args <- as.list(substitute(list(...)))[-1L]

  # Extract user-supplied weights if provided
  if ("weights" %in% names(dot_args)) {
    user_weights <- eval(dot_args$weights, envir = parent.frame())
    if (length(user_weights) != nrow(data)) {
      cli_abort("`weights` must have the same length as the number of rows in `data`.")
    }
    if (any(user_weights <= 0)) {
      cli_abort("`weights` must be positive (> 0).")
    }
  } else {
    user_weights <- rep(1, nrow(data))
  }

  if ("offset" %in% names(dot_args)) {
    if (!is.character(dot_args$offset)) {
      cli_abort("Please use a character value for 'offset' (indicating the column name) for cross validation.")
    }
    .offset <- eval(dot_args$offset)
  } else {
    .offset <- NULL
  }

  if (k_folds > 1) {
    # data in kth fold get weight of 0:
    fold_weights <- ifelse(data$cv_fold == 1L, 0, 1)
  } else {
    fold_weights <- rep(1, nrow(data))
  }
  if (lfo) fold_weights <- ifelse(data$cv_fold == 1L, 1, 0)

  # Combine user weights with fold weights
  weights <- user_weights * fold_weights

  if (use_initial_fit) {
    # run model on first fold to get starting values:

    if (!constant_mesh) {
      if (lfo) {
        dat_fit <- data[data$cv_fold == 1L, , drop = FALSE]
      } else {
        dat_fit <- data[data$cv_fold != 1L, , drop = FALSE]
      }
      # Create mesh on training data, then update with full data
      mesh_args[["data"]] <- dat_fit
      mesh_train <- do.call(make_mesh, mesh_args)
      mesh_args[["data"]] <- data
      mesh_args[["mesh"]] <- mesh_train$mesh
      mesh <- do.call(make_mesh, mesh_args)
    } else {
      mesh <- spde
    }
    dot_args <- list(dot_args)[[1]]
    dot_args$offset <- NULL
    dot_args$weights <- NULL
    .args <- c(list(
      data = data, formula = formula, time = time, mesh = mesh,
      weights = weights, offset = .offset
    ), dot_args)
    fit1 <- do.call(sdmTMB, .args)
  }

  fit_func <- function(k) {
    if (lfo) {
      fold_weights <- ifelse(data$cv_fold <= k, 1, 0)
    } else {
      # data in kth fold get weight of 0:
      fold_weights <- ifelse(data$cv_fold == k, 0, 1)
    }
    # Combine user weights with fold weights
    weights <- user_weights * fold_weights

    if (k == 1L && use_initial_fit) {
      object <- fit1
    } else {
      if (!constant_mesh) {
        if (lfo) {
          dat_fit <- data[data$cv_fold <= k, , drop = FALSE]
        } else {
          dat_fit <- data[data$cv_fold != k, , drop = FALSE]
        }
        # Create mesh on training data, then update with full data
        mesh_args[["data"]] <- dat_fit
        mesh_train <- do.call(make_mesh, mesh_args)
        mesh_args[["data"]] <- data
        mesh_args[["mesh"]] <- mesh_train$mesh
        mesh <- do.call(make_mesh, mesh_args)
      } else {
        mesh <- spde
      }
      dot_args <- as.list(substitute(list(...)))[-1L] # re-evaluate here! issue #54
      dot_args <- list(...)
      dot_args$offset <- NULL
      dot_args$weights <- NULL
      args <- c(list(
        data = data, formula = formula, time = time, mesh = mesh, offset = .offset,
        weights = weights, previous_fit = if (use_initial_fit) fit1 else NULL
      ), dot_args)
      object <- do.call(sdmTMB, args)
    }

    if (lfo) {
      cv_data <- data[data$cv_fold == (k + lfo_forecast), , drop = FALSE]
    } else {
      cv_data <- data[data$cv_fold == k, , drop = FALSE]
    }

    temporal_covariate_diffusion <- !is.null(object$covariate_diffusion_parsed) &&
      isTRUE(object$covariate_diffusion_parsed$needs_time)

    # FIXME: only use TMB report() below to be faster!
    # predict for withheld data:
    # cli_inform("Testing on data fold {k}.")
    if (temporal_covariate_diffusion) {
      predicted_full <- predict(object, newdata = data, type = "response",
        offset = if (!is.null(.offset)) data[[.offset]] else rep(0, nrow(data)))
      match_idx <- match(cv_data[["_sdm_order_"]], predicted_full[["_sdm_order_"]])
      predicted <- predicted_full[match_idx, , drop = FALSE]
    } else {
      predicted <- predict(object, newdata = cv_data, type = "response",
        offset = if (!is.null(.offset)) cv_data[[.offset]] else rep(0, nrow(cv_data)))
    }

    cv_data$cv_predicted <- predicted$est
    response <- get_response(object$formula[[1]])
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted

    # Calculate deviance residuals for held-out data
    dev_resids <- tryCatch({
      residuals(object, type = "deviance")
    }, error = function(e) NULL)

    # Match deviance residuals to held-out data
    if (!is.null(dev_resids) && !all(dev_resids == 0)) {
      # Get the row indices from cv_data (held-out data)
      cv_order <- cv_data[["_sdm_order_"]]
      # Get the row indices from object$data (fitted data)
      obj_order <- object$data[["_sdm_order_"]]
      # Match held-out indices to fitted data indices
      match_idx <- match(cv_order, obj_order)
      # Extract corresponding deviance residuals
      cv_data$cv_deviance_resid <- dev_resids[match_idx]
    } else {
      # Deviance residuals not available for this family
      cv_data$cv_deviance_resid <- NA_real_
    }

    # FIXME: get LFO working with the TMB report() option below!
    # calculate log likelihood for each withheld observation:
    # trickery to get the log likelihood of the withheld data directly
    # from the TMB report():
    if (!lfo) {
      tmb_data <- object$tmb_data

      # Reverse weights: training (weight != 0) → 0, held-out (weight == 0) → user_weight
      tmb_data$weights_i <- ifelse(tmb_data$weights_i == 0, user_weights, 0)

      new_tmb_obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = get_pars(object),
        map = object$tmb_map,
        random = object$tmb_random,
        DLL = "sdmTMB",
        silent = TRUE
      )
      lp <- object$tmb_obj$env$last.par.best
      r <- new_tmb_obj$report(lp)
      cv_loglik <- -1 * r$jnll_obs
      # Extract log-likelihoods for held-out observations (where reversed weights > 0)
      cv_data$cv_loglik <- cv_loglik[tmb_data$weights_i > 0]
    } else { # old method; doesn't work with delta models!
      cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    }

    ## test
    # x2 <- ll_sdmTMB(object, withheld_y, withheld_mu)
    # identical(round(cv_data$cv_loglik, 6), round(x2, 6))
    # cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)

    list(
      data = cv_data,
      model = if (save_models) object else NULL,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  }

  if (requireNamespace("future.apply", quietly = TRUE) && parallel) {
    message(
      "Running fits with `future.apply()`.\n",
      "Set a parallel `future::plan()` to use parallel processing."
    )
    if (!is.null(future_globals)) {
      fg <- structure(TRUE, add = future_globals)
    } else {
      fg <- TRUE
    }
    if (lfo) {
      out <- future.apply::future_lapply(seq_len(lfo_validations), fit_func, future.seed = TRUE, future.globals = fg)
    } else {
      out <- future.apply::future_lapply(seq_len(k_folds), fit_func, future.seed = TRUE, future.globals = fg)
    }
  } else {
    message(
      "Running fits sequentially.\n",
      "Install the future and future.apply packages,\n",
      "set a parallel `future::plan()`, and set `parallel = TRUE` to use parallel processing."
    )
    if (lfo) {
      out <- lapply(seq_len(lfo_validations), fit_func)
    } else {
      out <- lapply(seq_len(k_folds), fit_func)
    }
  }

  models <- if (save_models) lapply(out, `[[`, "model") else NULL
  data <- lapply(out, `[[`, "data")
  fold_cv_ll <- vapply(data, function(.x) sum(.x$cv_loglik), FUN.VALUE = numeric(1L))
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  data[["_sdmTMB_time"]] <- NULL
  row.names(data) <- NULL
  pdHess <- vapply(out, `[[`, "pdHess", FUN.VALUE = logical(1L))
  max_grad <- vapply(out, `[[`, "max_gradient", FUN.VALUE = numeric(1L))
  converged <- all(pdHess)
  out <- list(
    data = data,
    models = models,
    fold_loglik = fold_cv_ll,
    sum_loglik = sum(data$cv_loglik),
    converged = converged,
    pdHess = pdHess,
    max_gradients = max_grad
  )
  `class<-`(out, "sdmTMB_cv")
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

#' @export
#' @import methods
print.sdmTMB_cv <- function(x, ...) {
  nconverged <- sum(x$pdHess)
  nfolds <- length(x$pdHess)
  cat(paste0("Cross validation of sdmTMB models with ", nfolds, " folds.\n"))
  cat("\n")
  if (is.null(x$models)) {
    cat("Models were not saved (save_models = FALSE).\n")
    cat("Only prediction and log likelihood results are available.\n")
    cat("\n")
  } else {
    cat("Summary of the first fold model fit:\n")
    cat("\n")
    print(x$models[[1]])
    cat("\n")
    cat("Access the rest of the models in a list element named `models`.\n")
    cat("E.g. `object$models[[2]]` for the 2nd fold model fit.\n")
    cat("\n")
  }
  cat(paste0(nconverged, " out of ", nfolds, " models are consistent with convergence.\n"))
  cat("Figure out which folds these are in the `converged` list element.\n")
  cat("\n")
  cat(paste0("Out-of-sample log likelihood for each fold: ", paste(round(x$fold_loglik, 2), collapse = ", "), ".\n"))
  cat("Access these values in the `fold_loglik` list element.\n")
  cat("\n")
  cat("Sum of out-of-sample log likelihoods:", round(x$sum_loglik, 2), "\n")
  cat("More positive values imply better out-of-sample prediction.\n")
  cat("Access this value in the `sum_loglik` list element.\n")
}

#' Convert [sdmTMB_cv()] objects to \pkg{sf} format for spatial assessment with \pkg{waywiser}
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Converts cross-validation results to an [sf::sf()] object for use with
#' spatial model assessment tools such as those in the \pkg{waywiser} package.
#' This enables multi-scale spatial assessment of model predictions.
#'
#' @param object An object of class `sdmTMB_cv` from [sdmTMB_cv()].
#' @param ll_names Column names for longitude and latitude in the original data.
#'   **Note the order: longitude first, then latitude.**
#' @param ll_crs The coordinate reference system (CRS) for the `ll_names` columns.
#'   Defaults to `4326` (WGS84 lon/lat).
#' @param utm_crs The projected coordinate reference system (CRS) for the output
#'   sf object. By default (if you're feeling lucky!) automatically detected
#'   using [get_crs()] based on `ll_names`. Can be manually specified as an
#'   EPSG code (e.g., `32609`) or any format accepted by [sf::st_crs()].
#'
#' @return An [sf::sf()] object with POINT geometry containing:
#' \describe{
#'   \item{truth}{The observed response values}
#'   \item{estimate}{The cross-validated predictions}
#'   \item{geometry}{Spatial point locations}
#' }
#'
#' @details
#' This function is particularly useful for assessing spatial models at multiple
#' scales using the \pkg{waywiser} package. After converting to sf format, you
#' can use functions like [waywiser::ww_multi_scale()] to evaluate how model
#' performance changes when predictions are aggregated to different spatial
#' scales.
#'
#' For delta/hurdle models, the combined predictions are returned (e.g., the
#' product of the encounter probability and positive catch rate).
#'
#' @seealso [sdmTMB_cv()], [get_crs()], \url{https://sdmTMB.github.io/sdmTMB/articles/cross-validation.html}
#'
#' @export
#' @examplesIf require("sf", quietly = TRUE) && require("waywiser", quietly = TRUE)
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 12)
#'
#' # Run cross-validation
#' set.seed(123)
#' m_cv <- sdmTMB_cv(
#'   density ~ depth_scaled,
#'   data = pcod_2011,
#'   mesh = mesh,
#'   family = tweedie(),
#'   k_folds = 2
#' )
#'
#' # Convert with default auto-detected CRS based on lon/lat columns:
#' cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"))
#'
#' # Or manually specify the desired UTM CRS:
#' cv_sf <- cv_to_waywiser(m_cv, ll_names = c("lon", "lat"), utm_crs = 32609)
#'
#' # Use with waywiser for multi-scale assessment
#' waywiser::ww_multi_scale(
#'   cv_sf,
#'   truth,    # column name (unquoted)
#'   estimate, # column name (unquoted)
#'   n = list(c(5, 5), c(2, 2)) # 5x5 and 2x2 grid blocks
#' )
cv_to_waywiser <- function(object,
                           ll_names = c("longitude", "latitude"),
                           ll_crs = 4326,
                           utm_crs = get_crs(object$data, ll_names)) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli_abort("The sf package must be installed to use this function.")
  }

  if (is.null(object$models)) {
    cli_abort(c(
      "Models were not saved during cross-validation.",
      "i" = "Set `save_models = TRUE` in `sdmTMB_cv()` to use `cv_to_waywiser()`."
    ))
  }

  # Get response variable name from the first model's formula
  formulas <- object$models[[1]]$formula
  if (is.list(formulas)) {
    # Delta model - use first formula
    formula <- formulas[[1]]
  } else {
    # Non-delta model
    formula <- formulas
  }
  response <- get_response(formula)

  # Extract data
  cv_data <- object$data

  # Check that required columns exist
  if (!response %in% names(cv_data)) {
    cli_abort(c(
      "Response variable '{response}' not found in CV data.",
      "i" = "This may indicate an issue with the sdmTMB_cv object."
    ))
  }
  if (!"cv_predicted" %in% names(cv_data)) {
    cli_abort(c(
      "'cv_predicted' column not found in CV data.",
      "i" = "This may indicate an issue with the sdmTMB_cv object."
    ))
  }
  if (!all(ll_names %in% names(cv_data))) {
    cli_abort(c(
      "Longitude/latitude columns '{paste(ll_names, collapse = ', ')}' not found in CV data.",
      "i" = "Check that `ll_names` matches your data column names."
    ))
  }

  # Create data frame with required columns including lon/lat coordinates
  point_df <- data.frame(
    truth = cv_data[[response]],
    estimate = cv_data$cv_predicted
  )
  # Add the lon/lat coordinate columns
  point_df[[ll_names[1]]] <- cv_data[[ll_names[1]]]
  point_df[[ll_names[2]]] <- cv_data[[ll_names[2]]]

  # Convert to sf object with lon/lat CRS, then transform to UTM
  points <- sf::st_as_sf(point_df,
    coords = ll_names,
    crs = ll_crs
  )
  points <- sf::st_transform(points, utm_crs)

  points
}
