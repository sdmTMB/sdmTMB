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
  Biomass_KG = tweedie(link = "log"),
  Encounter = binomial(link = "cloglog"),
  Count = poisson(link = "log")
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
    paste0("tinyVAST: ", paste(names(tv_beta), collapse = ", ")),
    paste0("sdmTMB: ", paste(names(sdm_beta), collapse = ", "))
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
      fam <- fit$internal$family[[i]]
      phi_val <- exp(seg[1])
      if ("Gamma" %in% fam$family) {
        phi_val <- exp(-2 * seg[1])
      }
      out$phi[i] <- phi_val
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
    paste0("tinyVAST: ", paste(tv_params$family, collapse = ", ")),
    paste0("sdmTMB: ", paste(sdm_params$family, collapse = ", "))
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

index_cmp <- NULL
if (requireNamespace("sf", quietly = TRUE)) {
  data(red_snapper_shapefile, package = "tinyVAST")
  sf_grid <- sf::st_make_grid(red_snapper_shapefile, cellsize = c(0.2, 0.2))
  sf_grid <- sf::st_intersection(sf_grid, red_snapper_shapefile)
  sf_grid <- sf::st_make_valid(sf_grid)
  if (length(sf_grid) == 0) {
    message("Index comparison skipped: grid intersection returned 0 cells.")
  } else {
    grid_coords <- sf::st_coordinates(sf::st_centroid(sf_grid))
    areas_km2 <- as.numeric(sf::st_area(sf_grid)) / 1e6
    years <- sort(unique(red_snapper$Year))
    tv_idx <- numeric(length(years))
    sdm_idx <- numeric(length(years))
    for (i in seq_along(years)) {
      year <- years[[i]]
      newdata_idx <- data.frame(
        Lat = grid_coords[, "Y"],
        Lon = grid_coords[, "X"],
        Year = year,
        Data_type = "Biomass_KG",
        AreaSwept_km2 = mean(red_snapper$AreaSwept_km2),
        var = "logdens"
      )
      tv_out <- integrate_output(
        fit_tv,
        newdata = newdata_idx,
        area = areas_km2,
        getsd = FALSE,
        bias.correct = FALSE
      )
      tv_idx[i] <- unname(tv_out["Estimate"])
      pred <- predict(
        fit_sdm,
        newdata = newdata_idx,
        offset = log(newdata_idx$AreaSwept_km2),
        return_tmb_object = TRUE
      )
      sdm_out <- get_index(pred, area = areas_km2, bias_correct = FALSE)
      sdm_idx[i] <- sdm_out$est
    }
    index_cmp <- summarize_pred_diff(
      tv_idx,
      sdm_idx,
      as.character(years),
      "index_total",
      tol,
      use_relative = TRUE
    )
  }
} else {
  message("Index comparison skipped: package 'sf' not available.")
}

comparison <- rbind(beta_cmp, phi_cmp, tweedie_cmp, ll_cmp, pred_link_cmp, pred_mu_cmp, index_cmp)

red_snapper_delta <- droplevels(subset(red_snapper, Data_type == "Biomass_KG"))
red_snapper_delta$var <- "logdens"

Family_delta <- delta_gamma(type = "poisson-link")
formula_delta <- Response_variable ~ factor(Year) + offset(log(AreaSwept_km2))
formula_delta_sdm <- Response_variable ~ factor(Year)

fit_tv_delta <- tinyVAST(
  data = red_snapper_delta,
  formula = formula_delta,
  delta_options = list(formula = ~ .),
  family = Family_delta,
  distribution_column = "dist",
  spatial_domain = NULL,
  space_term = NULL,
  spacetime_term = NULL,
  time_term = NULL,
  space_columns = c("Lon", "Lat"),
  time_column = "Year",
  variable_column = "var",
  control = ctrl_tv
)

fit_sdm_delta <- sdmTMB(
  data = red_snapper_delta,
  formula = formula_delta_sdm,
  family = Family_delta,
  offset = log(red_snapper_delta$AreaSwept_km2),
  spatial = "off",
  spatiotemporal = "off",
  control = ctrl_sdm
)

tv_beta_delta1 <- fit_tv_delta$internal$parlist$alpha_j
names(tv_beta_delta1) <- colnames(fit_tv_delta$tmb_inputs$tmb_data$X_ij)

sdm_beta_delta1 <- fit_sdm_delta$parlist$b_j
names(sdm_beta_delta1) <- colnames(fit_sdm_delta$tmb_data$X_ij[[1]])

if (!identical(names(tv_beta_delta1), names(sdm_beta_delta1))) {
  msg <- c(
    "Fixed-effect names differ between tinyVAST and sdmTMB (poisson-link delta, component 1).",
    paste0("tinyVAST: ", paste(names(tv_beta_delta1), collapse = ", ")),
    paste0("sdmTMB: ", paste(names(sdm_beta_delta1), collapse = ", "))
  )
  stop(paste(msg, collapse = "\n"))
}

tv_beta_delta2 <- fit_tv_delta$internal$parlist$alpha2_j
names(tv_beta_delta2) <- colnames(fit_tv_delta$tmb_inputs$tmb_data$X2_ij)

sdm_beta_delta2 <- fit_sdm_delta$parlist$b_j2
names(sdm_beta_delta2) <- colnames(fit_sdm_delta$tmb_data$X_ij[[2]])

if (!identical(names(tv_beta_delta2), names(sdm_beta_delta2))) {
  msg <- c(
    "Fixed-effect names differ between tinyVAST and sdmTMB (poisson-link delta, component 2).",
    paste0("tinyVAST: ", paste(names(tv_beta_delta2), collapse = ", ")),
    paste0("sdmTMB: ", paste(names(sdm_beta_delta2), collapse = ", "))
  )
  stop(paste(msg, collapse = "\n"))
}

beta_delta1_cmp <- compare_vec(tv_beta_delta1, sdm_beta_delta1, "fixed_effects_delta1")
beta_delta2_cmp <- compare_vec(tv_beta_delta2, sdm_beta_delta2, "fixed_effects_delta2")

phi_tv_delta <- exp(-2 * fit_tv_delta$internal$parlist$log_sigma[1])
phi_sdm_delta <- exp(fit_sdm_delta$parlist$ln_phi[2])
phi_delta_cmp <- compare_vec(
  setNames(phi_tv_delta, "phi"),
  setNames(phi_sdm_delta, "phi"),
  "phi_delta"
)

ll_tv_delta <- -fit_tv_delta$opt$objective
ll_sdm_delta <- as.numeric(logLik(fit_sdm_delta))
ll_delta_cmp <- data.frame(
  term = "logLik",
  tinyVAST = ll_tv_delta,
  sdmTMB = ll_sdm_delta,
  diff = ll_tv_delta - ll_sdm_delta,
  abs_diff = abs(ll_tv_delta - ll_sdm_delta),
  rel_diff = NA_real_,
  within_tol = abs(ll_tv_delta - ll_sdm_delta) <= tol,
  check = "logLik_delta",
  stringsAsFactors = FALSE
)

pred_link_tv_delta <- predict(fit_tv_delta, newdata = red_snapper_delta, what = "p_g")
pred_link_sdm_delta_df <- predict(
  fit_sdm_delta,
  newdata = red_snapper_delta,
  offset = log(red_snapper_delta$AreaSwept_km2),
  type = "link"
)
pred_link_sdm_delta <- pred_link_sdm_delta_df$est1 + pred_link_sdm_delta_df$est2
pred_link_delta_cmp <- summarize_pred_diff(
  pred_link_tv_delta,
  pred_link_sdm_delta,
  rep("delta_gamma", nrow(red_snapper_delta)),
  "pred_p_g_delta",
  tol
)

pred_mu_tv_delta <- predict(fit_tv_delta, newdata = red_snapper_delta, what = "mu_g")
pred_mu_sdm_delta <- predict(
  fit_sdm_delta,
  newdata = red_snapper_delta,
  offset = log(red_snapper_delta$AreaSwept_km2),
  type = "response"
)$est
pred_mu_delta_cmp <- summarize_pred_diff(
  pred_mu_tv_delta,
  pred_mu_sdm_delta,
  rep("delta_gamma", nrow(red_snapper_delta)),
  "pred_mu_g_delta",
  tol,
  use_relative = TRUE
)

comparison <- rbind(
  comparison,
  beta_delta1_cmp,
  beta_delta2_cmp,
  phi_delta_cmp,
  ll_delta_cmp,
  pred_link_delta_cmp,
  pred_mu_delta_cmp
)

red_snapper_pl <- droplevels(subset(red_snapper, Data_type %in% c("Encounter", "Biomass_KG")))
red_snapper_pl$Data_type <- relevel(red_snapper_pl$Data_type, ref = "Biomass_KG")
red_snapper_pl$var <- "logdens"

Family_pl <- list(
  Biomass_KG = delta_gamma(type = "poisson-link"),
  Encounter = binomial(link = "cloglog")
)
delta_formula_pl <- ~ .

fit_tv_pl <- tinyVAST(
  data = red_snapper_pl,
  formula = formula,
  delta_options = list(formula = delta_formula_pl),
  family = Family_pl,
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

fit_sdm_pl <- sdmTMB(
  data = red_snapper_pl,
  formula = formula_sdm,
  family = Family_pl,
  distribution_column = "Data_type",
  offset = log(red_snapper_pl$AreaSwept_km2),
  spatial = "off",
  spatiotemporal = "off",
  control = ctrl_sdm
)

tv_beta_pl <- fit_tv_pl$internal$parlist$alpha_j
names(tv_beta_pl) <- colnames(fit_tv_pl$tmb_inputs$tmb_data$X_ij)

sdm_beta_pl <- fit_sdm_pl$parlist$b_j
names(sdm_beta_pl) <- colnames(fit_sdm_pl$tmb_data$X_ij[[1]])

if (!identical(names(tv_beta_pl), names(sdm_beta_pl))) {
  msg <- c(
    "Fixed-effect names differ between tinyVAST and sdmTMB (poisson-link delta).",
    paste0("tinyVAST: ", paste(names(tv_beta_pl), collapse = ", ")),
    paste0("sdmTMB: ", paste(names(sdm_beta_pl), collapse = ", "))
  )
  stop(paste(msg, collapse = "\n"))
}

tv_params_pl <- extract_tv_params(fit_tv_pl)
sdm_params_pl <- extract_sdm_params(fit_sdm_pl)

if (!identical(tv_params_pl$family, sdm_params_pl$family)) {
  msg <- c(
    "Family ordering differs between tinyVAST and sdmTMB (poisson-link delta).",
    paste0("tinyVAST: ", paste(tv_params_pl$family, collapse = ", ")),
    paste0("sdmTMB: ", paste(sdm_params_pl$family, collapse = ", "))
  )
  stop(paste(msg, collapse = "\n"))
}

ll_tv_pl <- -fit_tv_pl$opt$objective
ll_sdm_pl <- as.numeric(logLik(fit_sdm_pl))

beta_cmp_pl <- compare_vec(tv_beta_pl, sdm_beta_pl, "fixed_effects_pl")
phi_cmp_pl <- compare_vec(
  setNames(tv_params_pl$phi, tv_params_pl$family),
  setNames(sdm_params_pl$phi, sdm_params_pl$family),
  "phi_pl"
)
ll_cmp_pl <- data.frame(
  term = "logLik",
  tinyVAST = ll_tv_pl,
  sdmTMB = ll_sdm_pl,
  diff = ll_tv_pl - ll_sdm_pl,
  abs_diff = abs(ll_tv_pl - ll_sdm_pl),
  rel_diff = NA_real_,
  within_tol = abs(ll_tv_pl - ll_sdm_pl) <= tol,
  check = "logLik_pl",
  stringsAsFactors = FALSE
)

pred_link_tv_pl <- predict(fit_tv_pl, newdata = red_snapper_pl, what = "p_g")
pred_link_tv_pl_p1 <- predict(fit_tv_pl, newdata = red_snapper_pl, what = "p1_g")
pred_link_tv_pl_combined <- ifelse(
  red_snapper_pl$Data_type == "Biomass_KG",
  pred_link_tv_pl,
  pred_link_tv_pl_p1
)
pred_link_sdm_pl <- predict(
  fit_sdm_pl,
  newdata = red_snapper_pl,
  offset = log(red_snapper_pl$AreaSwept_km2),
  type = "link"
)$est
pred_link_cmp_pl <- summarize_pred_diff(
  pred_link_tv_pl_combined,
  pred_link_sdm_pl,
  red_snapper_pl$Data_type,
  "pred_p_g_pl",
  tol
)

pred_mu_tv_pl <- predict(fit_tv_pl, newdata = red_snapper_pl, what = "mu_g")
pred_mu_sdm_pl <- predict(
  fit_sdm_pl,
  newdata = red_snapper_pl,
  offset = log(red_snapper_pl$AreaSwept_km2),
  type = "response"
)$est
pred_mu_cmp_pl <- summarize_pred_diff(
  pred_mu_tv_pl,
  pred_mu_sdm_pl,
  red_snapper_pl$Data_type,
  "pred_mu_g_pl",
  tol,
  use_relative = TRUE
)

comparison <- rbind(
  comparison,
  beta_cmp_pl,
  phi_cmp_pl,
  ll_cmp_pl,
  pred_link_cmp_pl,
  pred_mu_cmp_pl
)
print(comparison, row.names = FALSE)

if (all(comparison$within_tol, na.rm = TRUE)) {
  message("All comparisons within tolerance: ", tol)
} else {
  message("Differences exceed tolerance: ", tol)
}
