#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tinyVAST)
})

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".")
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  stop("Install pkgload or devtools to load sdmTMB from source.")
}

data(red_snapper, package = "tinyVAST")
red_snapper$Data_type <- relevel(red_snapper$Data_type, ref = "Biomass_KG")
red_snapper$var <- "logdens"

Family <- list(
  Encounter = binomial(link = "cloglog"),
  Count = poisson(link = "log"),
  Biomass_KG = tweedie(link = "log")
)

formula <- Response_variable ~ Data_type + factor(Year) + offset(log(AreaSwept_km2))
formula_sdm <- Response_variable ~ Data_type + factor(Year)

ctrl_tv <- tinyVASTcontrol(
  nlminb_loops = 1,
  newton_loops = 0,
  getsd = FALSE,
  silent = TRUE,
  trace = 0,
  verbose = FALSE,
  calculate_deviance_explained = FALSE
)

fit_tv <- tinyVAST(
  data = red_snapper,
  formula = formula,
  family = Family,
  distribution_column = "Data_type",
  spatial_domain = NULL,
  space_term = NULL,
  spacetime_term = NULL,
  time_term = NULL,
  space_columns = c("Lon", "Lat"),
  time_column = "Year",
  variable_column = "var",
  control = ctrl_tv
)

ctrl_sdm <- sdmTMBcontrol(
  nlminb_loops = 1,
  newton_loops = 0,
  getsd = FALSE
)

fit_sdm <- sdmTMB(
  data = red_snapper,
  formula = formula_sdm,
  family = Family,
  distribution_column = "Data_type",
  offset = log(red_snapper$AreaSwept_km2),
  spatial = "off",
  spatiotemporal = "off",
  control = ctrl_sdm
)

tv_beta <- fit_tv$internal$parlist$alpha_j
names(tv_beta) <- colnames(fit_tv$tmb_inputs$tmb_data$X_ij)

sdm_beta <- fit_sdm$parlist$b_j
names(sdm_beta) <- colnames(fit_sdm$tmb_data$X_ij[[1]])

if (!identical(names(tv_beta), names(sdm_beta))) {
  msg <- c(
    "Fixed-effect names differ between tinyVAST and sdmTMB.",
    "tinyVAST: " %+% paste(names(tv_beta), collapse = ", "),
    "sdmTMB: " %+% paste(names(sdm_beta), collapse = ", ")
  )
  stop(paste(msg, collapse = "\n"))
}

extract_tv_params <- function(fit) {
  log_sigma <- fit$internal$parlist$log_sigma
  edims <- fit$tmb_inputs$tmb_data$Edims_ez
  out <- data.frame(
    family = names(fit$internal$family),
    phi = NA_real_,
    tweedie_p = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(fit$internal$family)) {
    start <- edims[i, "start"]
    len <- edims[i, "length"]
    if (len > 0) {
      seg <- log_sigma[(start + 1):(start + len)]
      out$phi[i] <- exp(seg[1])
      if (len > 1) {
        out$tweedie_p[i] <- 1 + plogis(seg[2])
      }
    }
  }
  out
}

extract_sdm_params <- function(fit) {
  ln_phi_e <- fit$parlist$ln_phi_e
  thetaf_e <- fit$parlist$thetaf_e
  starts_phi <- fit$tmb_data$ln_phi_start
  lens_phi <- fit$tmb_data$ln_phi_len
  starts_theta <- fit$tmb_data$thetaf_start
  lens_theta <- fit$tmb_data$thetaf_len

  out <- data.frame(
    family = names(fit$family),
    phi = NA_real_,
    tweedie_p = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(starts_phi)) {
    if (lens_phi[i] > 0) {
      idx <- starts_phi[i] + 1
      out$phi[i] <- exp(ln_phi_e[idx])
    }
    if (lens_theta[i] > 0) {
      idx <- starts_theta[i] + 1
      out$tweedie_p[i] <- 1 + plogis(thetaf_e[idx])
    }
  }
  out
}

tv_params <- extract_tv_params(fit_tv)
sdm_params <- extract_sdm_params(fit_sdm)

if (!identical(tv_params$family, sdm_params$family)) {
  msg <- c(
    "Family ordering differs between tinyVAST and sdmTMB.",
    "tinyVAST: " %+% paste(tv_params$family, collapse = ", "),
    "sdmTMB: " %+% paste(sdm_params$family, collapse = ", ")
  )
  stop(paste(msg, collapse = "\n"))
}

ll_tv <- -fit_tv$opt$objective
ll_sdm <- as.numeric(logLik(fit_sdm))

tol <- 1e-3
compare_vec <- function(a, b, label) {
  if (length(a) != length(b)) {
    stop("Length mismatch for comparison: ", label)
  }
  term <- names(a)
  if (is.null(term) || length(term) == 0L) term <- names(b)
  if (is.null(term) || length(term) == 0L) term <- as.character(seq_along(a))
  diff <- a - b
  data.frame(
    term = term,
    tinyVAST = as.numeric(a),
    sdmTMB = as.numeric(b),
    diff = as.numeric(diff),
    abs_diff = abs(as.numeric(diff)),
    rel_diff = NA_real_,
    within_tol = abs(as.numeric(diff)) <= tol,
    check = label,
    stringsAsFactors = FALSE
  )
}

beta_cmp <- compare_vec(tv_beta, sdm_beta, "fixed_effects")
phi_cmp <- compare_vec(
  setNames(tv_params$phi, tv_params$family),
  setNames(sdm_params$phi, sdm_params$family),
  "phi"
)
tweedie_cmp <- compare_vec(
  setNames(tv_params$tweedie_p, tv_params$family),
  setNames(sdm_params$tweedie_p, sdm_params$family),
  "tweedie_p"
)

ll_cmp <- data.frame(
  term = "logLik",
  tinyVAST = ll_tv,
  sdmTMB = ll_sdm,
  diff = ll_tv - ll_sdm,
  abs_diff = abs(ll_tv - ll_sdm),
  rel_diff = NA_real_,
  within_tol = abs(ll_tv - ll_sdm) <= tol,
  check = "logLik",
  stringsAsFactors = FALSE
)

summarize_pred_diff <- function(tv, sdm, dist, label, tol, use_relative = FALSE) {
  split_idx <- split(seq_along(dist), dist)
  rows <- lapply(names(split_idx), function(name) {
    idx <- split_idx[[name]]
    tv_vals <- as.numeric(tv[idx])
    sdm_vals <- as.numeric(sdm[idx])
    diff <- tv_vals - sdm_vals
    if (all(is.na(diff))) {
      max_abs <- NA_real_
      max_rel <- NA_real_
      mean_tv <- NA_real_
      mean_sdm <- NA_real_
      mean_diff <- NA_real_
    } else {
      max_abs <- max(abs(diff), na.rm = TRUE)
      denom <- pmax(abs(tv_vals), abs(sdm_vals), 1e-12)
      max_rel <- max(abs(diff) / denom, na.rm = TRUE)
      mean_tv <- mean(tv_vals, na.rm = TRUE)
      mean_sdm <- mean(sdm_vals, na.rm = TRUE)
      mean_diff <- mean(diff, na.rm = TRUE)
    }
    data.frame(
      term = paste0(label, "_", name),
      tinyVAST = mean_tv,
      sdmTMB = mean_sdm,
      diff = mean_diff,
      abs_diff = max_abs,
      rel_diff = max_rel,
      within_tol = if (use_relative) isTRUE(max_rel <= tol) else isTRUE(max_abs <= tol),
      check = label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pred_link_tv <- predict(fit_tv, newdata = red_snapper, what = "p_g")
pred_link_sdm <- predict(
  fit_sdm,
  newdata = red_snapper,
  offset = log(red_snapper$AreaSwept_km2),
  type = "link"
)$est
pred_link_cmp <- summarize_pred_diff(
  pred_link_tv,
  pred_link_sdm,
  red_snapper$Data_type,
  "pred_p_g",
  tol
)

pred_mu_tv <- predict(fit_tv, newdata = red_snapper, what = "mu_g")
pred_mu_sdm <- predict(
  fit_sdm,
  newdata = red_snapper,
  offset = log(red_snapper$AreaSwept_km2),
  type = "response"
)$est
pred_mu_cmp <- summarize_pred_diff(
  pred_mu_tv,
  pred_mu_sdm,
  red_snapper$Data_type,
  "pred_mu_g",
  tol,
  use_relative = TRUE
)

comparison <- rbind(beta_cmp, phi_cmp, tweedie_cmp, ll_cmp, pred_link_cmp, pred_mu_cmp)
print(comparison, row.names = FALSE)

if (all(comparison$within_tol, na.rm = TRUE)) {
  message("All comparisons within tolerance: ", tol)
} else {
  message("Differences exceed tolerance: ", tol)
}
