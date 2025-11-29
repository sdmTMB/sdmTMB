# from glmmTMB:
named_list <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(nm <- names(L))) {
    nm <- snm
  }
  if (any(nonames <- nm == "")) {
    nm[nonames] <- snm[nonames]
  }
  stats::setNames(L, nm)
}

print_model_info <- function(x) {
  delta <- isTRUE(x$family$delta)
  spatial_only <- as.logical(x$tmb_data$spatial_only)
  fit_by <- if (isTRUE(x$reml)) "REML" else "ML"
  if (all(spatial_only)) {
    title <- paste0("Spatial model fit by ", fit_by, " ['sdmTMB']\n")
  } else {
    title <- paste0("Spatiotemporal model fit by ", fit_by, " ['sdmTMB']\n")
  }
  if (all(x$spatial == "off") && all(x$spatiotemporal == "off")) {
    title <- paste0("Model fit by ", fit_by, " ['sdmTMB']\n")
  }

  aniso <- as.logical(x$tmb_data$anisotropy)
  # covariance <- paste0("Covariance: ", if (aniso) "anisotropic" else "isotropic", "\n")
  covariance <- paste0(if (aniso) "anisotropic" else "isotropic")

  formula <- paste0("Formula: ", extract_call_name(x$call$formula), "\n")

  time_name <- extract_call_name(x$call$time)
  if (!is.null(time_name) && time_name != "NULL") {
    time <- paste0("Time column: ", time_name, "\n")
    time <- gsub('\\"', "", time)
    time <- gsub("\\'", "", time)
  } else {
    time <- NULL
  }

  mesh <- paste0("Mesh: ", extract_call_name(x$call$mesh), " (", covariance, " covariance)\n")
  data <- paste0("Data: ", extract_call_name(x$call$data), "\n")

  if ("clean_name" %in% names(x$family)) {
    overall_family <- x$family$clean_name
  } else {
    overall_family <- paste0(x$family$family[1], "(link = '", x$family$link[1], "')")
  }
  overall_family <- paste0("Family: ", overall_family, "\n")

  if (delta) {
    family1 <- paste0("Family: ", x$family$family[1], "(link = '", x$family$link[1], "')")
    family2 <- paste0("Family: ", x$family$family[2], "(link = '", x$family$link[2], "')")
  } else {
    family1 <- NULL
    family2 <- NULL
  }

  criterion <- paste0(fit_by, " criterion at convergence: ", mround(x$model$objective, 3), "\n")

  named_list(
    delta,
    spatial_only,
    title,
    time,
    mesh,
    formula,
    data,
    family1,
    family2,
    overall_family,
    criterion,
    covariance
  )
}

print_main_effects <- function(x, m = 1) {
  b <- as.data.frame(tidy(x, model = m, silent = TRUE))
  b$estimate <- round(b$estimate, 2L)
  b$std.error <- round(b$std.error, 2L)
  mm <- cbind(b$estimate, b$std.error)
  colnames(mm) <- c("coef.est", "coef.se")
  row.names(mm) <- b$term
  mm
}

print_smooth_effects <- function(x, m = 1, edf = NULL, silent = FALSE) {
  sr <- x$sd_report
  sr_se <- as.list(sr, "Std. Error")
  sr_est <- as.list(sr, "Estimate")
  if (x$tmb_data$has_smooths) {
    bs_se <- round(sr_se$bs[, m], 2L)
    bs <- round(sr_est$bs[, m], 2L)
    sm <- parse_smoothers(formula = x$formula[[m]], data = x$data)
    sm_names_orig <- unlist(lapply(sm$Zs, function(x) attr(x, "s.label")))
    sm_names <- gsub("s\\(", "", sm_names_orig)
    sm_names <- gsub("t2\\(", "", sm_names)
    sm_names <- gsub("\\)$", "", sm_names)

    sm_names_bs <- sapply(seq_along(sm_names), function(i) {
      if (grepl("^t2\\(", sm_names_orig[i])) {
        paste0("t2(", sm_names[i], ")")
      } else {
        paste0("s", sm_names[i])
      }
    })
    sm_classes <- unlist(sm$classes)

    xx <- lapply(seq_along(sm_names_bs), function(i) {
      .x <- sm_names_bs[i]
      # For t2() smooths, don't expand - just use the base name
      if (grepl("^t2\\(", sm_names_orig[i])) {
        return(.x)
      }
      # original approach for regular univariate smooths
      n_sm <- grep(",", .x) + 1
      if (length(n_sm)) {
        .x <- paste(.x, seq_len(n_sm), sep = "_")
        gsub(",", "", .x)
      } else {
        .x
      }
    })

    # Expand sm_classes to match the expanded names
    sm_classes_expanded <- unlist(mapply(function(cls, n_expanded) {
      rep(cls, n_expanded)
    }, sm_classes, lengths(xx), SIMPLIFY = FALSE))

    # For coefficient names (sd_name = FALSE) adds indices,
    #   "t2(X,Y)" to "t2(X,Y)_1"
    # For SD names (is_sd_name = TRUE) adding indices after parenthesis
    #   "sds(X,Y)" to "sds(X,Y)_1"
    index_sm_names <- function(x, orig_names, is_sd_name = FALSE) {
      names_indexed <- character(length(x))
      i <- 1
      while (i <= length(x)) {
        if (grepl("^t2\\(", orig_names[i])) { # t2 case
          current_name <- x[i]
          j <- i
          while (j <= length(x) && x[j] == current_name) {
            j <- j + 1
          }
          count <- j - i
          if (count > 1) {
            if (is_sd_name) {
              # SD case: add index after closing parenthesis
              names_indexed[i:(j-1)] <- paste0(current_name, "_", seq_len(count))
            } else {
              # For coef names keep commas and add indices
              names_indexed[i:(j-1)] <- paste0(current_name, "_", seq_len(count))
            }
          } else {
            names_indexed[i] <- current_name
          }
          i <- j
        } else {
          names_indexed[i] <- x[i]
          i <- i + 1
        }
      }
      return(names_indexed)
    }

    # Helper function for smooth SD names - just returns one name per component
    # For regular s(): "depth" becomes sds(depth) or ln_SD_s(depth)
    # For smooths with by: depth):year2003 becomes sds(depth):year2003
    # For t2 smooths: returns one name per component, e.g., [sdt2(X,Y), sdt2(X,Y), sdt2(X,Y)]
    # These will be indexed later to become ..._1, ..._2, ..._3
    expand_smooth_sd_names <- function(prefix, sm_names, sm_names_orig) {
      lapply(seq_along(sm_names), function(i) {
        has_colon <- grepl(":", sm_names[i])
        # Check if original was t2() smoother and modify prefix accordingly
        if (grepl("^t2\\(", sm_names_orig[i])) {
          # Replace 's' pattern in prefix with 't2' for t2 smoothers
          # Handle both "sds(" and "ln_SD_s(" patterns
          prefix_modified <- gsub("^sds\\(", "sdt2(", prefix)
          prefix_modified <- gsub("_s\\(", "_t2(", prefix_modified)
          paste0(prefix_modified, sm_names[i], ifelse(has_colon, "", ")"))
        } else {
          paste0(prefix, sm_names[i], ifelse(has_colon, "", ")"))
        }
      })
    }

    # Create the SD name lists
    sm_names_sds_list <- expand_smooth_sd_names("sds(", sm_names, sm_names_orig)
    ln_sm_names_sds_list <- expand_smooth_sd_names("ln_SD_s(", sm_names, sm_names_orig)

    # Now index all the names
    sm_names_bs <- index_sm_names(unlist(xx),
                                  sm_names_orig,
                                  is_sd_name = FALSE)
    sm_names_sds <- index_sm_names(unlist(sm_names_sds_list),
                                   sm_names_orig,
                                   is_sd_name = TRUE)
    ln_sm_names_sds <- index_sm_names(unlist(ln_sm_names_sds_list),
                                      sm_names_orig,
                                      is_sd_name = TRUE)

    mm_sm <- cbind(bs, bs_se)

    .sm_names_bs <- sm_names_bs[!sm_classes_expanded %in% c("cc.smooth.spec", "cs.smooth.spec")]
    if (length(.sm_names_bs) != nrow(mm_sm)) {
      if (nrow(mm_sm) > 0L) {
        msg <- c(
          "Smoother fixed effect matrix names could not be retrieved.",
          "This could be from using t2(), setting m != 2 in s(), or other unanticipated s() arguments.",
          "This does not affect model fitting.",
          "We'll use generic covariate names ('scovariate') here instead.",
          ""
        )
        if (!silent) cli_inform(msg)
        row.names(mm_sm) <- paste0("scovariate-", seq_len(nrow(mm_sm)))
      } else {
        mm_sm <- NULL
      }
    } else {
      row.names(mm_sm) <- .sm_names_bs
    }

    smooth_sds <- round(exp(sr_est$ln_smooth_sigma[, m]), 2L)
    re_sm_mat <- matrix(NA_real_, nrow = length(smooth_sds), ncol = 1L)
    re_sm_mat[, 1] <- smooth_sds
    # Convert sds() and sdt2() to sd__s() and sd__t2() to match tidy() convention
    print_names <- gsub("^sds\\(", "sd__s(", sm_names_sds)
    print_names <- gsub("^sdt2\\(", "sd__t2(", print_names)
    rownames(re_sm_mat) <- print_names
    colnames(re_sm_mat) <- "Std. Dev."

    # second matrix for ran_pars
    ln_re_sm_mat <- matrix(NA_real_, nrow = length(smooth_sds), ncol = 2L)
    ln_re_sm_mat[, 1] <- sr_est$ln_smooth_sigma[, m]
    ln_re_sm_mat[, 2] <- sr_se$ln_smooth_sigma[, m]
    rownames(ln_re_sm_mat) <- ln_sm_names_sds
    colnames(ln_re_sm_mat) <- c("estimate", "std.error")

    if (!is.null(edf)) {
      if (is_delta(x)) {
        lp_regex <- paste0("^", m, "LP-s\\(")
        edf <- edf[grepl(lp_regex, names(edf))]
      } else {
        edf <- edf[grepl("^s\\(", names(edf))]
      }
      edf <- round(edf, 2)
      re_sm_mat <- cbind(re_sm_mat, matrix(edf, ncol = 1))
      colnames(re_sm_mat)[2] <- "EDF"
    }
  } else {
    re_sm_mat <- NULL
    mm_sm <- NULL
    ln_re_sm_mat <- NULL
  }
  list(smooth_effects = mm_sm, smooth_sds = re_sm_mat, ln_sd_estimates = ln_re_sm_mat)
}

print_int_slope_re <- function(x, m = 1) {
  if (sum(x$tmb_data$n_re_groups)) {
    v <- tidy(x, effects = "ran_vcov", model = m)
    cnms <- x$split_formula[[m]]$re_cov_terms$cnms
    ll <- vapply(cnms, length, FUN.VALUE = 1L)
    mmc2 <- unlist(cnms, use.names = FALSE)
    mmc1 <- lapply(names(ll), \(na) rep(na, ll[[na]])) |> unlist()
    mmc1[duplicated(mmc1)] <- ""
    mmsd <- lapply(v[[1]], \(x) {
      if (ncol(x) == 2L) {
        c(x[1,1], x[2,2])
      } else {
        c(x[1,1])
      }
    })
    mmsd <- unlist(mmsd, use.names = FALSE)
    mmvar <- mmsd^2
    mmcor <- lapply(v[[1]], \(x) {
      if (ncol(x) == 2L) {
        c("", mround(x[2,1], 2))
      } else {
        ""
      }
    })
    mmcor <- unlist(mmcor, use.names = FALSE)
    mm <- cbind(mmc1, mmc2, mround(mmvar, 2), mround(mmsd, 2), mmcor)
    colnames(mm) <- c("Groups", "Name", "Variance", "Std.Dev.", "Corr")
    rownames(mm) <- rep("", nrow(mm))
    if (all(mm[,"Corr"] == "")) mm <- mm[,-5]
    mm
  } else {
    NULL
  }
}

print_iid_re <- function(x, m = 1) {
  .tidy <- tidy(x, "ran_pars", model = m, silent = TRUE)
  if ("sigma_G" %in% .tidy$term) {
    re_int_names <- x$split_formula[[1]]$barnames # TODO DELTA HARDCODED TO 1
    re_int_mat <- matrix(NA_real_, nrow = length(re_int_names), ncol = 1L)
    re_int_mat[, 1L] <- round(.tidy$estimate[.tidy$term == "sigma_G"], 2L)
    rownames(re_int_mat) <- re_int_names
    colnames(re_int_mat) <- "Std. Dev."
  } else {
    re_int_mat <- NULL
  }
  re_int_mat
}

print_time_varying <- function(x, m = 1) {
  if (!is.null(x$time_varying)) {
    sr <- x$sd_report
    sr_se <- as.list(sr, "Std. Error")
    sr_est <- as.list(sr, "Estimate")
    b_rw_t_est <- sr_est$b_rw_t[, , m]
    b_rw_t_se <- sr_se$b_rw_t[, , m]
    tv_names <- colnames(model.matrix(x$time_varying, x$data))
    mm_tv <- cbind(round(as.numeric(b_rw_t_est), 2L), round(as.numeric(b_rw_t_se), 2L))
    colnames(mm_tv) <- c("coef.est", "coef.se")
    time_slices <- x$time_lu$time_from_data
    row.names(mm_tv) <- paste(rep(tv_names, each = length(time_slices)), time_slices, sep = "-")

    # append AR(1) parameters if they are estimated
    p <- tidy(x, effects = "ran_pars", model = m, silent = TRUE)
    if(any(p$term == "rho_time")) {
      rho_tv <- as.matrix(p[p$term=="rho_time",c("estimate","std.error")])
      rho_tv[] <- round(rho_tv, 2L)
      colnames(rho_tv) <- colnames(mm_tv)
      row.names(rho_tv) <- paste("rho", tv_names, sep = "-")
      mm_tv <- rbind(mm_tv, rho_tv)
    }
  } else {
    mm_tv <- NULL
  }
  mm_tv
}

print_range <- function(x, m = 1L, digits = 2L) {
  b <- tidy(x, effects = "ran_pars", model = m, silent = TRUE)
  range <- b$estimate[b$term == "range"]
  if (is.null(range)) {
    return(NULL)
  }

  range <- mround(range, digits)
  range_text <- if (x$tmb_data$share_range[m]) {
    paste0("Mat\u00e9rn range: ", range[1], "\n")
  } else {
    paste0(
      "Mat\u00e9rn range (spatial): ", range[1], "\n",
      "Mat\u00e9rn range (spatiotemporal): ", range[2], "\n"
    )
  }

  if (as.logical(x$tmb_data$anisotropy)) {
    range_text <- print_anisotropy(x = x, m = m)
  }

  if (x$spatial[m] == "off" && x$spatiotemporal[m] == "off") {
    range_text <- NULL
  }

  range_text
}

print_anisotropy <- function(x, m = 1L, digits = 1L, return_dat = FALSE) {
  aniso_df <- plot_anisotropy(x, return_data = TRUE)
  aniso_df$degree <- aniso_df$angle * 180 / pi
  aniso_df_st <- aniso_df_sp <- NULL

  if (isTRUE(x$family$delta)) {
    aniso_df_sp <- aniso_df[aniso_df$random_field == "spatial" &
        aniso_df$model_num == m, ][1, c("a", "b", "degree")]
    aniso_df_st <- aniso_df[aniso_df$random_field == "spatiotemporal" &
        aniso_df$model_num == m, ][1L, c("a", "b", "degree")]
  } else {
    if (x$spatial[m] != "off") {
      aniso_df_sp <- aniso_df[aniso_df$random_field == "spatial", ][1L, c("a", "b", "degree")]
    }
    if (x$spatiotemporal[m] != "off") {
      aniso_df_st <- aniso_df[aniso_df$random_field == "spatiotemporal", ][1L, c("a", "b", "degree")]
    }
  }

  if (x$spatial[m] != "off") {
    aniso_df_sp[1:2] <- mround(aniso_df_sp[1:2], digits)
    aniso_df_sp[3] <- mround(aniso_df_sp[3], 0L)
    aniso_df_sp[4] <- paste0(
      "Mat\u00e9rn anisotropic range (spatial): ",
      aniso_df_sp[2], " to ", aniso_df_sp[1], " at ", aniso_df_sp[3], " deg.", "\n"
    )
  }
  if (x$spatiotemporal[m] != "off") {
    aniso_df_st[1:2] <- mround(aniso_df_st[1:2], digits)
    aniso_df_st[3] <- mround(aniso_df_st[3], 0L)
    aniso_df_st[4] <- paste0(
      "Mat\u00e9rn anisotropic range (spatiotemporal): ",
      aniso_df_st[2], " to ", aniso_df_st[1], " at ", aniso_df_st[3], " deg.", "\n"
    )
  }

  # Spatial only
  if (x$spatial[m] == "on" && x$spatiotemporal[m] == "off") {
    range_text <- aniso_df_sp[[4]]
  }
  # Spatiotemporal only
  if (x$spatial[m] == "off" && x$spatiotemporal[m] != "off") {
    range_text <- aniso_df_st[[4]]
  }
  # Spatial and spatiotemporal shared
  if (x$tmb_data$share_range[m] == 1L && x$spatial[m] == "on" && x$spatiotemporal[m] != "off") {
    range_text <- aniso_df_sp[[4]]
  }
  # Spatial and spatiotemporal NOT shared
  if (x$tmb_data$share_range[m] == 0L && x$spatial[m] == "on" && x$spatiotemporal[m] != "off") {
    range_text <- paste0(aniso_df_sp[[4]], aniso_df_st[[4]])
  }
  if (return_dat) return(list(sp = aniso_df_sp, st = aniso_df_st))
  range_text
}

print_other_parameters <- function(x, m = 1L) {
  b <- tidy(x, "ran_pars", model = m, silent = TRUE)

  get_term_text <- function(term_name = "", pretext = "") {
    b2 <- as.list(x$sd_report, what = "Estimate")
    if (term_name %in% b$term) {
      a <- mround(b$estimate[b$term == term_name], 2L)
      a <- paste0(pretext, ": ", a, "\n")
    } else if (term_name %in% names(b2)) {
      a <- mround(b2[[term_name]], 2L)
      a <- paste0(pretext, ": ", a, "\n")
    } else {
      a <- ""
    }
    a
  }

  phi <- get_term_text("phi", "Dispersion parameter")
  tweedie_p <- get_term_text("tweedie_p", "Tweedie p")
  student_df <- get_term_text("student_df", "Student-t df")
  gengamma_par <- if ('gengamma' %in% family(x)[[m]]) {
    get_term_text("gengamma_Q", "Generalized gamma Q")
    } else ""
  sigma_O <- get_term_text("sigma_O", "Spatial SD")
  xtra <- if (x$spatiotemporal[m] == "ar1") "marginal " else ""
  sigma_E <- get_term_text("sigma_E",
    paste0("Spatiotemporal ", xtra, toupper(x$spatiotemporal[m]), " SD"))
  rho <- get_term_text("rho", "Spatiotemporal AR1 correlation (rho)")

  if ("sigma_Z" %in% b$term) {
    # tidy() takes sigma_Z from the sdreport,
    # which condenses them to unique values if mapped, so:
    sigma_Z <- x$tmb_obj$report(x$tmb_obj$env$last.par.best)$sigma_Z
    sigma_Z <- sigma_Z[,m,drop=TRUE]
    a <- mround(sigma_Z, 2L)
    sigma_Z <- paste0("Spatially varying coefficient SD (", x$spatial_varying,  "): ", a, "\n")
    sigma_Z <- gsub("\\(\\(", "\\(", sigma_Z) # ((Intercept))
    sigma_Z <- gsub("\\)\\)", "\\)", sigma_Z) # ((Intercept))
    sigma_Z <- paste(sigma_Z, collapse = "")
  } else {
    sigma_Z <- ""
  }

  named_list(phi, tweedie_p, student_df, sigma_O, sigma_E, sigma_Z, rho, gengamma_par)
}

print_header <- function(x) {
  info <- print_model_info(x)
  cat(info$title)
  cat(info$formula)
  cat(info$mesh)
  cat(info$time)
  cat(info$data)
  cat(info$overall_family)
}

print_one_model <- function(x, m = 1, edf = FALSE, silent = FALSE) {
  if (edf) {
    .edf <- suppressMessages(cAIC(x, what = "EDF"))
  } else {
    .edf <- NULL
  }
  info <- print_model_info(x)
  main <- print_main_effects(x, m = m)
  smooth <- print_smooth_effects(x, m = m, edf = .edf, silent = TRUE)
  # iid_re <- print_iid_re(x, m = m)
  re <- print_int_slope_re(x, m = m)
  tv <- print_time_varying(x, m = m)
  range <- print_range(x, m = m)
  other <- print_other_parameters(x, m = m)

  if (m == 1) cat(info$family1, "\n")
  if (m == 2) cat(info$family2, "\n")

  if (!is.null(re)) {
    cat("Random intercepts and/or slopes:\n\n")
    cat("Conditional model:\n")
    print(re, quote = FALSE)
    cat("\n")
  }

  cat("Conditional model:\n")
  print(main)
  cat("\n")

  if (!is.null(smooth$smooth_sds)) {
    cat("Smooth terms:\n")
    print(smooth$smooth_sds)
    cat("\n")
  }

  if (!is.null(tv)) {
    cat("Time-varying parameters:\n")
    print(tv)
    cat("\n")
  }

  cat(other$phi)
  cat(other$tweedie_p)
  cat(other$student_df)
  cat(other$gengamma_par)
  cat(other$rho)
  cat(range)
  cat(other$sigma_O)
  cat(other$sigma_Z)
  cat(other$sigma_E)
}
print_footer <- function(x) {
  info <- print_model_info(x)
  cat(info$criterion)

  cat("\nSee ?tidy.sdmTMB to extract these values as a data frame.\n")

  if (as.logical(x$tmb_data$anisotropy)) {
    cat("See ?plot_anisotropy to plot the anisotropic range.\n")
  }
  sink(tempfile())
  suppressWarnings(suppressMessages(s <- sanity(x)))
  sink()
  if (!all(unlist(s))) {
    cat("\n**Possible issues detected! Check output of sanity().**\n")
  }
}

#' @export
#' @import methods
print.sdmTMB <- function(x, ...) {
  reinitialize(x)
  lp <- x$tmb_obj$env$last.par.best
  r <- x$tmb_obj$report(lp)

  delta <- isTRUE(x$family$delta)
  print_header(x)
  if (delta) cat("\nDelta/hurdle model 1: -----------------------------------\n")
  print_one_model(x, 1, ...)
  if (delta) {
    cat("\nDelta/hurdle model 2: -----------------------------------\n")
    print_one_model(x, 2, ...)
  }
  if (delta) cat("\n")
  print_footer(x)
}
