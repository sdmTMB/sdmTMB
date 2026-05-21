#' Turn sdmTMB model output into a tidy data frame
#'
#' @param x Output from [sdmTMB()].
#' @param effects A character value. One of `"fixed"` ('fixed' or main-effect
#'   parameters), `"ran_pars"` (standard deviations, spatial range, and other
#'   random effect and dispersion-related terms), `"ran_vals"` (individual
#'   random intercepts or slopes, if included; behaves like `ranef()`), or `"ran_vcov"` (list
#'   of variance covariance matrices for the random effects, by model and group).
#' @param conf.int Include a confidence interval?
#' @param conf.level Confidence level for CI.
#' @param exponentiate Whether to exponentiate the fixed-effect coefficient
#'   estimates and confidence intervals.
#' @param model Which model to tidy if a delta model (1 or 2). The `model` will be
#'   ignored when effects is `"ran_vals"` (all returned in a single dataframe)
#'
#' @param silent Omit any messages?
#' @param ... Extra arguments (not used).
#'
#' @return A data frame
#' @details
#' Follows the conventions of the \pkg{broom} and \pkg{broom.mixed} packages.
#'
#' Currently, `effects = "ran_pars"` also includes dispersion-related terms
#' (e.g., `phi`), which are not actually associated with random effects.
#'
#' Standard errors for spatial variance terms fit in log space (e.g., variance
#' terms, range, or parameters associated with the observation error) are
#' omitted to avoid confusion. Confidence intervals are still available.
#'
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom stats plogis
#' @examples
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#' tidy(fit, "ran_pars", conf.int = TRUE)
#'
#' pcod_2011$fyear <- as.factor(pcod_2011$year)
#' fit <- sdmTMB(density ~ poly(depth_scaled, 2, raw = TRUE) + (1 | fyear),
#'   data = pcod_2011, mesh = pcod_mesh_2011,
#'   family = tweedie()
#' )
#' tidy(fit, "ran_vals")

.tidy_multi_family_param_rows <- function(est, se, family_spec, crit, conf.int = TRUE) {
  rows <- list()

  add_rows <- function(slot_name, term_name, estimate_fun, std_error_fun = NULL, conf_fun = NULL) {
    slot <- family_spec$param_slot[[slot_name]]
    used <- which(!is.na(slot))
    if (!length(used)) {
      return(NULL)
    }
    out <- lapply(used, function(f) {
      idx <- slot[[f]]
      estimate <- estimate_fun(idx)
      std.error <- if (is.null(std_error_fun)) NA_real_ else std_error_fun(idx)
      row <- data.frame(
        group_name = family_spec$family_labels[[f]],
        term = term_name,
        estimate = estimate,
        std.error = std.error,
        stringsAsFactors = FALSE
      )
      if (conf.int) {
        bounds <- if (is.null(conf_fun)) c(NA_real_, NA_real_) else conf_fun(idx)
        row$conf.low <- bounds[[1]]
        row$conf.high <- bounds[[2]]
      }
      row
    })
    do.call(rbind, out)
  }

  rows$phi <- add_rows(
    slot_name = "ln_phi",
    term_name = "phi",
    estimate_fun = function(idx) as.numeric(est$phi[[idx]]),
    std_error_fun = function(idx) {
      if (is.null(se$phi) || length(se$phi) < idx) NA_real_ else as.numeric(se$phi[[idx]])
    },
    conf_fun = function(idx) {
      if (is.null(est$ln_phi) || is.null(se$ln_phi) || length(se$ln_phi) < idx) {
        return(c(NA_real_, NA_real_))
      }
      c(
        exp(est$ln_phi[[idx]] - crit * se$ln_phi[[idx]]),
        exp(est$ln_phi[[idx]] + crit * se$ln_phi[[idx]])
      )
    }
  )
  rows$tweedie_p <- add_rows(
    slot_name = "thetaf",
    term_name = "tweedie_p",
    estimate_fun = function(idx) plogis(est$thetaf[[idx]]) + 1,
    std_error_fun = function(idx) {
      if (is.null(se$tweedie_p) || length(se$tweedie_p) < idx) NA_real_ else as.numeric(se$tweedie_p[[idx]])
    },
    conf_fun = function(idx) {
      if (is.null(est$thetaf) || is.null(se$thetaf) || length(se$thetaf) < idx) {
        return(c(NA_real_, NA_real_))
      }
      c(
        plogis(est$thetaf[[idx]] - crit * se$thetaf[[idx]]) + 1,
        plogis(est$thetaf[[idx]] + crit * se$thetaf[[idx]]) + 1
      )
    }
  )
  rows$student_df <- add_rows(
    slot_name = "ln_student_df",
    term_name = "student_df",
    estimate_fun = function(idx) exp(est$ln_student_df[[idx]]) + 1,
    std_error_fun = function(idx) {
      if (is.null(se$student_df) || length(se$student_df) < idx) NA_real_ else as.numeric(se$student_df[[idx]])
    },
    conf_fun = function(idx) {
      if (is.null(est$ln_student_df) || is.null(se$ln_student_df) || length(se$ln_student_df) < idx) {
        return(c(NA_real_, NA_real_))
      }
      c(
        exp(est$ln_student_df[[idx]] - crit * se$ln_student_df[[idx]]) + 1,
        exp(est$ln_student_df[[idx]] + crit * se$ln_student_df[[idx]]) + 1
      )
    }
  )
  rows$gengamma_Q <- add_rows(
    slot_name = "gengamma_Q",
    term_name = "gengamma_Q",
    estimate_fun = function(idx) as.numeric(est$gengamma_Q[[idx]]),
    std_error_fun = function(idx) {
      if (is.null(se$gengamma_Q) || length(se$gengamma_Q) < idx) NA_real_ else as.numeric(se$gengamma_Q[[idx]])
    },
    conf_fun = function(idx) {
      if (is.null(est$gengamma_Q) || is.null(se$gengamma_Q) || length(se$gengamma_Q) < idx) {
        return(c(NA_real_, NA_real_))
      }
      c(
        est$gengamma_Q[[idx]] - crit * se$gengamma_Q[[idx]],
        est$gengamma_Q[[idx]] + crit * se$gengamma_Q[[idx]]
      )
    }
  )

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(NULL)
  }
  do.call(rbind, rows)
}

tidy.sdmTMB <- function(x, effects = c("fixed", "ran_pars", "ran_vals", "ran_vcov"), model = 1,
                 conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE,
                 silent = FALSE, ...) {
  effects <- match.arg(effects)
  assert_that(is.logical(exponentiate))
  assert_that(is.logical(conf.int))
  if (conf.int) {
    assert_that(is.numeric(conf.level),
      conf.level > 0, conf.level < 1,
      length(conf.level) == 1,
      msg = "`conf.level` must be length 1 and between 0 and 1")
  }

  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  if (exponentiate) trans <- exp else trans <- I

  reinitialize(x)
  is_areal <- is_areal_fit(x)
  is_car <- is_car_fit(x)
  family_spec <- .object_family_spec(x, caller = "`tidy()`")
  multi_family <- family_spec$n_f > 1L

  delta <- family_spec$n_m == 2L
  assert_that(is.numeric(model))
  assert_that(length(model) == 1L)
  if (delta) assert_that(model %in% c(1, 2), msg = "`model` must be 1 or 2.")
  if (!delta) assert_that(model == 1, msg = "Only one model: `model` must be 1.")

  se_rep <- as.list(x$sd_report, "Std. Error", report = TRUE)
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  se <- as.list(x$sd_report, "Std. Error", report = FALSE)
  est <- as.list(x$sd_report, "Estimate", report = FALSE)

  se <- c(se, se_rep)
  est <- c(est, est_rep)
  # cleanup:
  est$epsilon_st <- NULL
  est$zeta_s <- NULL
  est$omega_s <- NULL
  est$ln_H_input <- NULL

  se$epsilon_st <- NULL
  se$zeta_s <- NULL
  se$omega_s <- NULL
  se$ln_H_input <- NULL

  subset_pars <- function(p, model) {
    .subset_model <- function(x) {
      if (is.null(x) || !length(x)) {
        return(NULL)
      }
      if (is.null(dim(x)) || length(dim(x)) == 1L) {
        if (length(x) >= model) {
          return(as.numeric(x[model]))
        }
        return(as.numeric(x))
      }
      if (length(dim(x)) >= 2L) {
        return(as.numeric(x[, model, drop = TRUE]))
      }
      as.numeric(x)
    }
    p$b_j <- if (model == 1) p$b_j else p$b_j2
    p$ln_tau_O <- .subset_model(p$ln_tau_O)
    p$ln_tau_Z <- .subset_model(p$ln_tau_Z)
    p$ln_tau_E <- .subset_model(p$ln_tau_E)
    p$ln_kappa <- .subset_model(p$ln_kappa)
    p$ln_tau_V <- .subset_model(p$ln_tau_V)
    p$ar1_phi <- .subset_model(p$ar1_phi)
    p$log_sigma_O <- .subset_model(p$log_sigma_O)
    p$log_sigma_E <- .subset_model(p$log_sigma_E)
    p$log_sigma_Z <- .subset_model(p$log_sigma_Z)
    p$log_range <- .subset_model(p$log_range)
    p$logit_rho_sar <- .subset_model(p$logit_rho_sar)

    p$range <- .subset_model(p$range)
    p$sigma_E <- .subset_model(p$sigma_E)
    p$sigma_O <- .subset_model(p$sigma_O)
    p$sigma_Z <- .subset_model(p$sigma_Z)
    p$rho_sar <- .subset_model(p$rho_sar)
    p$alpha_car <- .subset_model(p$alpha_car)
    if (!multi_family) {
      p$ln_phi <- .subset_model(p$ln_phi)
      p$phi <- .subset_model(p$phi)
      p$thetaf <- .subset_model(p$thetaf)
      p$tweedie_p <- .subset_model(p$tweedie_p)
      p$ln_student_df <- .subset_model(p$ln_student_df)
      p$student_df <- .subset_model(p$student_df)
      p$gengamma_Q <- .subset_model(p$gengamma_Q)
    }

    if (!is.null(p$rho_time_unscaled) && length(p$rho_time_unscaled)) {
      p$rho_time_unscaled <- .subset_model(p$rho_time_unscaled)
      p$rho_time <- 2 * plogis(p$rho_time_unscaled) - 1
    } else {
      p$rho_time <- numeric(0)
    }
    if (!is.null(p$logit_rho_sar) && length(p$logit_rho_sar) &&
        (is.null(p$rho_sar) || !length(p$rho_sar)) &&
        (is.null(p$alpha_car) || !length(p$alpha_car))) {
      if (is_car) {
        p$alpha_car <- plogis(p$logit_rho_sar)
      } else {
        p$rho_sar <- 2 * plogis(p$logit_rho_sar) - 1
      }
    }
    p
  }

  est <- subset_pars(est, model)
  se <- subset_pars(se, model)

  if (!multi_family && x$family$family[[model]] %in% c("binomial", "poisson")) {
    se$ln_phi <- NULL
    est$ln_phi <- NULL
    se$phi <- NULL
    est$phi <- NULL
  }

  ii <- 1

  # grab fixed effects:
  .formula <- x$split_formula[[model]]$form_no_bars
  .formula <- remove_s_and_t2(.formula)
  if (!"mgcv" %in% names(x)) x[["mgcv"]] <- FALSE
  fe_names <- colnames(model.matrix(.formula, x$data))
  if (!is.null(x$covariate_diffusion_data)) {
    fe_names <- c(fe_names, x$covariate_diffusion_data$term_coef_name)
  }

  b_j <- est$b_j[!fe_names == "offset", drop = TRUE]
  b_j_se <- se$b_j[!fe_names == "offset", drop = TRUE]
  fe_names <- fe_names[!fe_names == "offset"]
  out <- data.frame(term = fe_names, estimate = b_j, std.error = b_j_se, stringsAsFactors = FALSE)

  if (x$tmb_data$threshold_func > 0) {
    if (x$threshold_function == 1L) {
      par_name <- paste0(x$threshold_parameter, c("-slope", "-breakpt"))
      estimates <- est$b_threshold[,model,drop=TRUE]
      ses <- se$b_threshold[,model,drop=TRUE]
    } else {
      par_name <- paste0(x$threshold_parameter, c("-s50", "-s95", "-smax"))
      estimates <- c(est$s50[model], est$s95[model], est$s_max[model])
      ses <- c(se$s50[model], se$s95[model], se$s_max[model])
    }
    out <- rbind(
      out,
      data.frame(
        term = par_name, estimate = estimates,
        std.error = ses, stringsAsFactors = FALSE
      )
    )
  }

  if (x$tmb_data$has_smooths) {
    p <- print_smooth_effects(x, silent = FALSE, m = model)
    mm <- p$smooth_effects
    out <- rbind(
      out,
      data.frame(
        term = rownames(mm),
        estimate = mm[, "bs"],
        std.error = mm[, "bs_se"],
        stringsAsFactors = FALSE
      )
    )
  }

  if (conf.int) {
    out$conf.low <- as.numeric(trans(out$estimate - crit * out$std.error))
    out$conf.high <- as.numeric(trans(out$estimate + crit * out$std.error))
  }
  # must wrap in as.numeric() otherwise I() leaves 'AsIs' class that affects emmeans package
  out$estimate <- as.numeric(trans(out$estimate))
  if (exponentiate) out$std.error <- NULL

  out_re <- list()
  log_name <- character(0)
  name <- character(0)
  if (!is_areal) {
    log_name <- c(log_name, "log_range")
    name <- c(name, "range")
  }
  if (!multi_family && !isTRUE(is.na(x$tmb_map$ln_phi))) {
    log_name <- c(log_name, "ln_phi")
    name <- c(name, "phi")
  }
  if (x$tmb_data$include_spatial[model]) {
    log_name <- c(log_name, "log_sigma_O")
    name <- c(name, "sigma_O")
  }
  if (!x$tmb_data$spatial_only[model]) {
    log_name <- c(log_name, "log_sigma_E")
    name <- c(name, "sigma_E")
  }
  if (x$tmb_data$spatial_covariate) {
    log_name <- c(log_name, "log_sigma_Z")
    name <- c(name, "sigma_Z")
  }
  if (x$tmb_data$random_walk) {
    log_name <- c(log_name, "ln_tau_V")
    name <- c(name, "tau_V")
  }
  if (!all(est$rho_time == 0)) {
    log_name <- c(log_name, "rho_time_unscaled")
    name <- c(name, "rho_time")
  }
  if (is_areal && "logit_rho_sar" %in% names(est) && length(est$logit_rho_sar)) {
    log_name <- c(log_name, "logit_rho_sar")
    name <- c(name, if (is_car) "alpha_car" else "rho_sar")
  }

  j <- 0
  expects_range <- !is_areal &&
    (x$tmb_data$include_spatial[model] || !x$tmb_data$spatial_only[model])
  if (expects_range && !"log_range" %in% names(est)) {
    cli_warn("This model was fit with an old version of sdmTMB. Some parameters may not be available to the tidy() method. Re-fit the model with the current version of sdmTMB if you need access to any missing parameters.")
  }

  for (i in name) {
    j <- j + 1
    if (i %in% names(est)) {
      .e <- est[[log_name[j]]]
      .se <- se[[log_name[j]]]
      .e <- if (is.null(.e)) NA else .e
      .se <- if (is.null(.se)) NA else .se

      non_log_name <- gsub("ln_", "", gsub("log_", "", log_name))
      this <- non_log_name[j]
      if (this == "tau_V") this <- "sigma_V"
      if (this == "rho_time_unscaled") this <- "rho_time"
      if (this == "logit_rho_sar") this <- if (is_car) "alpha_car" else "rho_sar"

      this_se <- as.numeric(se[[this]])
      this_est <- as.numeric(est[[this]])
      if (!length(this_se) && length(this_est)) {
        this_se <- rep(NA_real_, length(this_est))
      }
      if (length(this_est) && !(all(this_se == 0) && all(this_est == 0))) {
        out_re[[i]] <- data.frame(
          term = i, estimate = this_est, std.error = this_se,
          conf.low = exp(.e - crit * .se),
          conf.high = exp(.e + crit * .se),
          stringsAsFactors = FALSE
        )
        if (this %in% c("rho_time", "rho_sar")) {
          out_re[[i]] <- data.frame(
            term = i,
            estimate = 2 * plogis(.e) - 1,
            # use delta method to get SE in normal space
            std.error = 2 * plogis(.e) * (1 - plogis(.e)) * .se,
            # don't use delta-method for CIs, because they can be outside (-1,1)
            conf.low = 2 * plogis(.e - crit * .se) - 1,
            conf.high = 2 * plogis(.e + crit * .se) - 1,
            stringsAsFactors = FALSE
          )
        }
        if (this == "alpha_car") {
          out_re[[i]] <- data.frame(
            term = i,
            estimate = plogis(.e),
            std.error = plogis(.e) * (1 - plogis(.e)) * .se,
            conf.low = plogis(.e - crit * .se),
            conf.high = plogis(.e + crit * .se),
            stringsAsFactors = FALSE
          )
        }
      }
      ii <- ii + 1
    }
  }
  discard <- unlist(lapply(out_re, function(x) length(x) == 1L)) # e.g. old models and phi
  out_re[discard] <- NULL

  if (!multi_family && "tweedie" %in% x$family$family) {
    out_re$tweedie_p <- data.frame(
      term = "tweedie_p", estimate = plogis(est$thetaf) + 1,
      std.error = se$tweedie_p, stringsAsFactors = FALSE)
    out_re$tweedie_p$conf.low <- plogis(est$thetaf - crit * se$thetaf) + 1
    out_re$tweedie_p$conf.high <- plogis(est$thetaf + crit * se$thetaf) + 1
    ii <- ii + 1
  }

  if (!multi_family && "student" %in% x$family$family) {
    # Check if df was fixed (mapped to NA) or estimated
    df_fixed <- !is.null(x$tmb_map$ln_student_df) && is.na(x$tmb_map$ln_student_df[1])
    if (df_fixed) {
      # Fixed df: report the fixed value with NA std.error
      df_value <- exp(est$ln_student_df) + 1
      out_re$student_df <- data.frame(
        term = "student_df", estimate = df_value,
        std.error = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        stringsAsFactors = FALSE)
    } else {
      # Estimated df: report with confidence intervals
      out_re$student_df <- data.frame(
        term = "student_df", estimate = exp(est$ln_student_df) + 1,
        std.error = se$student_df, stringsAsFactors = FALSE)
      out_re$student_df$conf.low <- exp(est$ln_student_df - crit * se$ln_student_df) + 1
      out_re$student_df$conf.high <- exp(est$ln_student_df + crit * se$ln_student_df) + 1
    }
    ii <- ii + 1
  }

  if ("ar1_phi" %in% names(est)) {
    ar_phi <- est$ar1_phi
    ar_phi_se <- se$ar1_phi
    rho_est <- 2 * stats::plogis(ar_phi) - 1
    rho_lwr <- 2 * stats::plogis(ar_phi - crit * ar_phi_se) - 1
    rho_upr <- 2 * stats::plogis(ar_phi + crit * ar_phi_se) - 1
    out_re[[ii]] <- data.frame(
      term = "rho", estimate = rho_est, std.error = NA,
      conf.low = rho_lwr, conf.high = rho_upr, stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  add_covariate_diffusion_parameter <- function(term_name, covariate_mask_name,
                                                display_name = term_name) {
    if (is.null(x$covariate_diffusion_data) || !length(x$covariate_diffusion_data$covariates)) {
      return(NULL)
    }
    covariates <- x$covariate_diffusion_data$covariates
    keep <- x$covariate_diffusion_data[[covariate_mask_name]]
    if (is.null(keep) || length(keep) != length(covariates)) {
      return(NULL)
    }
    keep_idx <- which(as.logical(keep))
    if (!length(keep_idx) || is.null(est[[term_name]])) {
      return(NULL)
    }
    estimates <- rep(NA_real_, length(covariates))
    ses <- rep(NA_real_, length(covariates))
    estimates[keep_idx] <- est[[term_name]]
    if (!is.null(se[[term_name]])) {
      ses[keep_idx] <- se[[term_name]]
    }
    out <- data.frame(
      term = paste0(display_name, "[", covariates, "]"),
      estimate = estimates,
      std.error = ses,
      conf.low = estimates - crit * ses,
      conf.high = estimates + crit * ses,
      stringsAsFactors = FALSE
    )
    out[keep_idx, , drop = FALSE]
  }
  covariate_diffusion_term_masks <- list(
    list(term_name = "kappaS_dl", display_name = "kappaS_cov_diff", covariate_mask_name = "covariate_has_spatial"),
    list(term_name = "kappaT_dl", display_name = "kappaT_cov_diff", covariate_mask_name = "covariate_has_temporal"),
    list(term_name = "kappaST_dl", display_name = "kappaST_cov_diff", covariate_mask_name = "covariate_has_spacetime"),
    list(term_name = "rhoT", display_name = "rhoT", covariate_mask_name = "covariate_has_temporal"),
    list(term_name = "MSD", display_name = "MSD", covariate_mask_name = "covariate_has_spatial"),
    list(term_name = "RMSD", display_name = "RMSD", covariate_mask_name = "covariate_has_spatial")
  )
  for (term_info in covariate_diffusion_term_masks) {
    term_df <- add_covariate_diffusion_parameter(
      term_info$term_name,
      term_info$covariate_mask_name,
      term_info$display_name
    )
    if (!is.null(term_df)) {
      out_re[[paste0("covariate_diffusion_", term_info$display_name)]] <- term_df
    }
  }

  if (all(!x$tmb_data$include_spatial) && all(x$tmb_data$spatial_only)) out_re$range <- NULL
  if (is_areal) out_re$range <- NULL

  out_re <- do.call("rbind", out_re)
  row.names(out_re) <- NULL
  if (multi_family) {
    family_param_rows <- .tidy_multi_family_param_rows(
      est = est,
      se = se,
      family_spec = family_spec,
      crit = crit,
      conf.int = conf.int
    )
    if (!is.null(family_param_rows)) {
      all_cols <- union(names(out_re), names(family_param_rows))
      for (nm in setdiff(all_cols, names(out_re))) out_re[[nm]] <- NA
      for (nm in setdiff(all_cols, names(family_param_rows))) family_param_rows[[nm]] <- NA
      out_re <- rbind(out_re[, all_cols, drop = FALSE], family_param_rows[, all_cols, drop = FALSE])
    }
  }

  if (identical(est$ln_tau_E, 0)) out_re <- out_re[out_re$term != "sigma_E", ]
  if (identical(est$ln_tau_V, 0)) out_re <- out_re[out_re$term != "sigma_V", ]
  if (identical(est$ln_tau_O, 0)) out_re <- out_re[out_re$term != "sigma_O", ]
  if (identical(est$ln_tau_Z, 0)) out_re <- out_re[out_re$term != "sigma_Z", ]
  if (is.na(x$tmb_map$ar1_phi[model])) out_re <- out_re[out_re$term != "rho", ]

  # Handle anisotropic ranges
  aniso_ranges <- get_anisotropic_ranges(x, m = model)
  if (!is.null(aniso_ranges) && "range" %in% out_re$term) {
    # Remove existing range row(s) and add anisotropic ranges
    out_re <- out_re[out_re$term != "range", ]
    aniso_df <- aniso_ranges_to_df(aniso_ranges)
    if (!is.null(aniso_df)) {
      # Add columns to match out_re structure
      aniso_df$std.error <- NA_real_
      aniso_df$conf.low <- NA_real_
      aniso_df$conf.high <- NA_real_
      out_re <- rbind(out_re, aniso_df)
    }
  }

  if (!conf.int) {
    out_re[["conf.low"]] <- NULL
    out_re[["conf.high"]] <- NULL
  }

  if (sum(x$tmb_data$n_re_groups) > 0L) { # we have random intercepts/slopes
    temp <- get_re_tidy_list(x, crit = crit, model = model, delta = delta)
    cov_mat_list <- list(est = temp$cov_matrices)
    if(conf.int) {
      cov_mat_list[["lo"]] <- temp$cov_matrices_lo
      cov_mat_list[["hi"]] <- temp$cov_matrices_hi
    }
    out_ranef <- temp$out_ranef

    # get names
    cnms <- x$split_formula[[model]]$re_cov_terms$cnms
    flattened_cov <- flatten_cov_output(cov_mat_list, cnms)
    flattened_cov$model <- NULL
    flattened_cov$term <- paste0("sd__", flattened_cov$term)
    if (!conf.int) {
      flattened_cov[["conf.low"]] <- NULL
      flattened_cov[["conf.high"]] <- NULL
    }
    out_re$group_name <- rep(NA, nrow(out_re))
    out_re <- rbind(out_re, flattened_cov)
  } else {
    cov_mat_list <- NULL
    out_ranef <- NULL
  }

  # optional smooth SDs models with smooths
  if (x$tmb_data$has_smooths) {
    p <- print_smooth_effects(x, m = model, silent = TRUE)
    ln_sds <- p$ln_sd_estimates
    # Convert from log-scale to natural scale and update term names
    term_names <- gsub("^ln_SD_s\\(", "sd__s(", row.names(ln_sds))
    term_names <- gsub("^ln_SD_t2\\(", "sd__t2(", term_names)
    # add on CIs - calculate in log space then transform
    ln_sds_df <- data.frame(term = term_names,
                            estimate = exp(ln_sds[,1]),
                            std.error = NA_real_,
                            conf.low = exp(ln_sds[,1] - crit*ln_sds[,2]),
                            conf.high = exp(ln_sds[,1] + crit*ln_sds[,2]),
                            group_name = NA)
    if (!conf.int) {
      ln_sds_df[["conf.low"]] <- NULL
      ln_sds_df[["conf.high"]] <- NULL
    }
    row.names(ln_sds_df) <- NULL
    if (!"group_name" %in% names(out_re)) out_re$group_name <- rep(NA, nrow(out_re))
    out_re <- rbind(out_re, ln_sds_df)
  }
  # optional time-varying random components
  if (!is.null(x$time_varying)) {
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    time_slices <- x$time_lu$time_from_data
    yrs <- rep(time_slices, times = length(tv_names))

    if (delta) {
      out_ranef_tv <- data.frame(
        model = model,
        term = paste0(rep(tv_names, each = length(time_slices)), ":", yrs),
        estimate = c(est$b_rw_t),
        std.error = c(se$b_rw_t),
        conf.low = c(est$b_rw_t) - crit * c(se$b_rw_t),
        conf.high = c(est$b_rw_t) + crit * c(se$b_rw_t),
        stringsAsFactors = FALSE
      )
    } else {
      out_ranef_tv <- data.frame(
        term = paste0(rep(tv_names, each = length(time_slices)), ":", yrs),
        estimate = c(est$b_rw_t),
        std.error = c(se$b_rw_t),
        conf.low = c(est$b_rw_t) - crit * c(se$b_rw_t),
        conf.high = c(est$b_rw_t) + crit * c(se$b_rw_t),
        stringsAsFactors = FALSE
      )
    }

    if(is.null(out_ranef)) {
      out_ranef <- out_ranef_tv
    } else {
      out_ranef_tv$group_name <- NA
      out_ranef_tv$level_ids <- NA
      out_ranef <- rbind(out_ranef, out_ranef_tv)
    }
  }

  # re-order names to match those in "ran_vals"
  if (delta) {
    out_re$model <- rep(model, nrow(out_re))
  }
  # group_name and term might not exist
  if (delta) {
    new_names <- c("model", "group_name", "term", "estimate", "std.error")
  } else {
    new_names <- c("group_name", "term", "estimate", "std.error")
  }
  if(conf.int) new_names <- c(new_names, "conf.low", "conf.high")
  new_names <- new_names[new_names %in% names(out_re)]
  out_re <- out_re[, new_names]
  # remove group_name if not used for simplicity:
  if (all(is.na(out_re$group_name))) out_re$group_name <- NULL

  out <- unique(out) # range can be duplicated
  out_re <- unique(out_re)

  if (requireNamespace("tibble", quietly = TRUE)) {
    frm <- tibble::as_tibble
  } else {
    frm <- as.data.frame
  }

  if (effects == "fixed") {
    return(frm(out))
  } else if (effects == "ran_vals") {
    return(frm(out_ranef))
  } else if (effects == "ran_pars") {
    return(frm(out_re))
  } else if (effects == "ran_vcov") {
    return(cov_mat_list)
  } else {
    cli_abort("The specified 'effects' type is not available.")
  }
}

# Convert anisotropic ranges list to a data frame
aniso_ranges_to_df <- function(aniso_ranges) {
  if (is.null(aniso_ranges) || length(aniso_ranges) == 0) {
    return(NULL)
  }

  # Build term-estimate pairs
  terms <- character()
  estimates <- numeric()

  for (name in names(aniso_ranges)) {
    suffix <- if (name == "min") "range_min"
              else if (name == "max") "range_max"
              else paste0("range_", name)
    terms <- c(terms, suffix)
    estimates <- c(estimates, aniso_ranges[[name]])
  }

  data.frame(
    term = terms,
    estimate = estimates,
    stringsAsFactors = FALSE
  )
}

# Internal helper to extract anisotropic range min/max values
# Returns NULL if not anisotropic, otherwise a list with range info
get_anisotropic_ranges <- function(x, m = 1L) {
  # Check if model has anisotropy enabled
  if (!as.logical(x$tmb_data$anisotropy)) {
    return(NULL)
  }

  # Check if H matrix is available
  H_exists <- any(grepl(
    pattern = "ln_H_input",
    x = names(x$sd_report$par.fixed),
    ignore.case = TRUE
  ))
  if (!H_exists) {
    return(NULL)
  }

  # Calculate anisotropy components using shared helper
  comp <- calculate_anisotropy_components(x, m = m)

  # Helper to calculate range from eigenvector
  rss <- function(V) sqrt(sum(V[1]^2 + V[2]^2))

  # Build output list
  result <- list()

  # Add spatial field ranges if spatial is on
  if (x$spatial[m] != "off") {
    result$spatial_min <- rss(comp$min_s)
    result$spatial_max <- rss(comp$maj_s)
  }

  # Add spatiotemporal field ranges if spatiotemporal is on
  if (x$spatiotemporal[m] != "off") {
    result$spatiotemporal_min <- rss(comp$min_st)
    result$spatiotemporal_max <- rss(comp$maj_st)
  }

  # Determine if we should use simpler names (only one field active, or shared range)
  only_spatial <- (x$spatial[m] != "off" && x$spatiotemporal[m] == "off")
  only_spatiotemporal <- (x$spatial[m] == "off" && x$spatiotemporal[m] != "off")
  shared_range <- (x$tmb_data$share_range[m] == 1L &&
                   x$spatial[m] == "on" &&
                   x$spatiotemporal[m] != "off")

  # Use simpler names if only one field is active or ranges are shared
  if (only_spatial || shared_range) {
    result <- list(
      min = result$spatial_min,
      max = result$spatial_max
    )
  } else if (only_spatiotemporal) {
    result <- list(
      min = result$spatiotemporal_min,
      max = result$spatiotemporal_max
    )
  }

  result
}

# Extract and format random effect estimates from sdmTMB model output
#
# This function extracts random effect estimates, including individual random intercepts
# and slopes, as well as covariance matrices, from an `sdmTMB` model output. It formats
# them into a structured list for further analysis.
#
# @param x An `sdmTMB` model object containing estimated random effects.
# @param crit The critical value for confidence interval computation,
#   typically derived from a normal distribution quantile, e.g.
#   crit = stats::qnorm(1 - (1 - conf.level) / 2)
# @importFrom stats aggregate

get_re_tidy_list <- function(x, crit, model = 1, delta = FALSE) {
  re_b_dfs <- add_model_index(x$split_formula, "re_b_df")
  re_b_df <- do.call(rbind, re_b_dfs)
  names(re_b_df)[names(re_b_df) == "group_indices"] <- "group_id"

  # this function just expands each row from start: end
  expand_row <- function(level_id, start, end, group_id, model) {
    seq_len <- end - start + 1
    data.frame(
      level_ids = rep(level_id, seq_len),
      index = seq(from = start, to = end),
      group_id = rep(group_id, seq_len),
      model = rep(model, seq_len)
    )
  }

  # apply to each row and combine the results
  expanded_rows <- Map(
    expand_row, re_b_df$level_ids,
    re_b_df$start, re_b_df$end,
    re_b_df$group_id, re_b_df$model
  )

  re_b_df <- do.call(rbind, expanded_rows) # list to df
  rownames(re_b_df) <- NULL # reset row names

  # this is all as before
  re_indx <- grep("re_b_pars", names(x$sd_report$value), fixed = TRUE)
  non_nas <- !is.na(x$tmb_map$re_b_pars) # parameters that don't get mapped off

  re_b_df$estimate <- x$sd_report$value[re_indx][non_nas]
  re_b_df$std.error <- x$sd_report$sd[re_indx][non_nas]
  re_b_df$conf.low <- re_b_df$estimate - crit * re_b_df$std.error
  re_b_df$conf.high <- re_b_df$estimate + crit * re_b_df$std.error
  re_b_df$index <- NULL
  re_b_df$group_name <- NA
  re_b_df$term <- NA
  for (i in seq_len(length(x$split_formula))) {
    groupnames <- names(x$split_formula[[i]]$re_cov_terms$cnms)
    for (j in seq_len(length(x$split_formula[[i]]$barnames))) {
      model_grp <- which(re_b_df$model == i & re_b_df$group_id == j)
      re_b_df$group_name[model_grp] <- groupnames[j]
      re_b_df$term[model_grp] <- rep(x$split_formula[[i]]$re_cov_terms$cnms[[j]], length.out = length(model_grp))
    }
  }
  group_key <- stats::aggregate(group_name ~ model + group_id, data = re_b_df, FUN = function(x) x[1])

  # more sensible re-ordering
  re_b_df$group_id <- NULL
  if (delta) {
    re_b_df <- re_b_df[, c("model", "group_name", "term", "level_ids", "estimate", "std.error", "conf.low", "conf.high")]
  } else {
    re_b_df <- re_b_df[, c("group_name", "term", "level_ids", "estimate", "std.error", "conf.low", "conf.high")]
  }
  # remove ":" in the level_ids
  re_b_df$level_ids <- sapply(strsplit(re_b_df$level_ids, ":"), function(x) x[2])
  out_ranef <- re_b_df
  row.names(out_ranef) <- NULL

  re_cov_dfs <- add_model_index(x$split_formula, "re_df")
  # add group names
  re_cov_dfs <- lapply(re_cov_dfs, function(df) {
    df$group <- row.names(df)
    df
  })
  re_cov_df <- do.call(rbind, re_cov_dfs)
  re_cov_df$rows <- re_cov_df$rows + 1 # increement rows/cols from 1, not 0
  re_cov_df$cols <- re_cov_df$cols + 1
  row.names(re_cov_df) <- NULL

  # make sure group index is included correctly
  for (i in seq_len(nrow(re_cov_df))) {
    indx <- which(group_key$model == re_cov_df$model[i] & group_key$group_id == (re_cov_df$group_indices[i] + 1))
    re_cov_df$group[i] <- group_key$group_name[indx]
  }

  re_indx <- grep("re_cov_pars", names(x$sd_report$value), fixed = TRUE)
  non_nas <- which(x$sd_report$value[re_indx] != 0) # remove parameter that gets mapped off
  re_cov_df$estimate <- x$sd_report$value[re_indx][non_nas]
  re_cov_df$std.error <- x$sd_report$sd[re_indx][non_nas]
  re_cov_df$conf.low <- re_cov_df$estimate - crit * re_cov_df$std.error
  re_cov_df$conf.high <- re_cov_df$estimate + crit * re_cov_df$std.error
  # the SD parameters are returned in log space -- calculate CIs in log space then transform
  sds <- which(re_cov_df$is_sd == 1) # index which elements are SDs
  # Calculate CIs in log space first
  re_cov_df$conf.low[sds] <- exp(re_cov_df$estimate[sds] - crit * re_cov_df$std.error[sds])
  re_cov_df$conf.high[sds] <- exp(re_cov_df$estimate[sds] + crit * re_cov_df$std.error[sds])
  # Transform estimates to natural space
  re_cov_df$estimate[sds] <- exp(re_cov_df$estimate[sds])

  re_cov_df <- re_cov_df[, c("rows", "cols", "model", "group", "estimate", "std.error", "conf.low", "conf.high")]
  cov_matrices_lo = create_cov_matrices(re_cov_df, col_name = "conf.low", model = model)
  cov_matrices_hi = create_cov_matrices(re_cov_df, col_name = "conf.high", model = model)
  cov_matrices_est = create_cov_matrices(re_cov_df, model = model)
  list(
    out_ranef = out_ranef,
    cov_matrices = cov_matrices_est,
    cov_matrices_lo = cov_matrices_lo,
    cov_matrices_hi = cov_matrices_hi
  )
}

create_cov_matrices <- function(df, col_name = "estimate", model = 1) {
  # Initialize an empty list to store the covariance matrices
  cov_matrices <- list()

  # Group by model and group, and then create a covariance matrix for each group
  for (m in model) {
    model_df <- df[df$model == m,,drop=FALSE]
    for (group in unique(model_df$group)) {
      group_df <- model_df[model_df$group == group,,drop=FALSE]
      # Create an empty matrix
      cov_matrix <- matrix(NA_real_, nrow = max(group_df$rows), ncol = max(group_df$cols))
      # Fill the matrix with the estimates
      for (i in seq_len(nrow(group_df))) {
        cov_matrix[group_df$rows[i], group_df$cols[i]] <- group_df[[col_name]][i]
      }
      # Add the matrix to the list
      cov_matrices[[paste("Model", m, "Group", group)]] <- cov_matrix
    }
  }
  cov_matrices
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' @rdname tidy.sdmTMB
#' @export
tidy.sdmTMB_cv <- function(x, ...) {
  if (is.null(x$models)) {
    cli_abort(c(
      "Models were not saved during cross-validation.",
      "i" = "Set `save_models = TRUE` in `sdmTMB_cv()` to use `tidy()`."
    ))
  }
  x <- x$models
  out <- lapply(seq_along(x), function(i) {
    df <- tidy.sdmTMB(x[[i]], ...)
    df$cv_split <- i # add a model index column
    df
  })
  do.call("rbind", out)
}

# Flatten a list of estimated random effects into a dataframe
#
# This function extracts random effect estimates from a list of matrices representing
# the mean, lower and upper confidence intervals. Ignoring the off-diagonal elements,
# it returns a dataframe of random effects estimates, their standard errors, and confidence intervals,
# consistent with the other tidy output from `sdmTMB` models.
#
# @param v A list of 3 random effects; for multiple grouping variables this will be a list of lists
# @param cnms A list of vectors of coefficient names for each set of grouping variable
flatten_cov_output <- function(v, cnms) {
  results <- list()
  idx <- 1

  for (name in names(v$est)) {
    est_mat <- v$est[[name]]
    lo_mat  <- v$lo[[name]]
    hi_mat  <- v$hi[[name]]
    # Split_name is the name of the model and group
    split_name <- strsplit(name, " +")[[1]]
    model <- as.integer(split_name[2])
    group_name <- split_name[4]

    term_names <- cnms[[group_name]]
    n <- nrow(est_mat)  # number of coefficients being estimated

    # Extract diagonal values only
    for (i in seq_len(n)) {
      est <- est_mat[i, i]
      if (is.na(est)) next

      conf_low <- if (!is.null(lo_mat)) lo_mat[i, i] else NA
      conf_high <- if (!is.null(hi_mat)) hi_mat[i, i] else NA

      results[[idx]] <- data.frame(
        model = model,
        group_name = group_name,
        term = term_names[i],
        estimate = est,
        std.error = NA, # we may want to change this and also report std.error?
        conf.low = conf_low,
        conf.high = conf_high,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  df <- do.call(rbind, results)
  rownames(df) <- NULL
  df
}
