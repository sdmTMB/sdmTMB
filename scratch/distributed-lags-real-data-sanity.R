suppressPackageStartupMessages({
  library(pkgload)
  library(TMB)
  library(Matrix)
})

# Load the local sdmTMB source tree rather than an installed package version.
pkgload::load_all(".", quiet = TRUE)

# The legacy template used as an external reference implementation.
template_cpp <- "scratch/spacetime_lag_2025_10_24.cpp"
template_dll <- sub("\\.cpp$", "", template_cpp)

# Compile/load the template DLL once for this script.
if (!file.exists(TMB::dynlib(template_dll))) {
  TMB::compile(template_cpp)
}
if (!(basename(template_dll) %in% names(getLoadedDLLs()))) {
  dyn.load(TMB::dynlib(template_dll))
}

# Build the data object used for all checks in this script.
# We convert the response to a transformed Gaussian target for fast fitting.
prepare_dataset <- function() {
  dat <- as.data.frame(sdmTMB::pcod_2011)
  response <- pmax(dat$density, 0)
  keep <- is.finite(dat$depth) &
    is.finite(response) &
    !is.na(dat$year)
  dat <- dat[keep, , drop = FALSE]
  response <- response[keep]
  dat$depth_scaled <- as.numeric(scale(dat$depth))
  dat$y_gaussian <- log1p(response)
  dat
}

# Fit a small real-data distributed-lag model and return the key diagnostics.
fit_real_model <- function() {
  dat <- prepare_dataset()
  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 20)

  fit <- sdmTMB::sdmTMB(
    y_gaussian ~ depth_scaled,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = ~ spatial(depth_scaled),
    control = sdmTMB::sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )

  x_names <- colnames(fit$tmb_data$X_ij[[1]])
  b_hat <- fit$tmb_obj$env$parList(fit$model$par)$b_j
  names(b_hat) <- x_names
  lag_name <- grep("^dl_", x_names, value = TRUE)

  list(
    dataset = "pcod_2011",
    n = nrow(dat),
    years = length(unique(dat$year)),
    objective = fit$model$objective,
    lag_name = lag_name,
    lag_beta = unname(b_hat[lag_name]),
    fit = fit
  )
}

# This matches the lag-template convention: 1 on first sub-diagonal.
make_shift_time_matrix <- function(n_t) {
  Matrix::sparseMatrix(
    i = seq.int(2L, n_t),
    j = seq.int(1L, n_t - 1L),
    x = 1,
    dims = c(n_t, n_t)
  )
}

# Project template output from mesh vertices to observations.
# sdmTMB stores indices as 0-based; R indexing is 1-based.
project_template_to_observations <- function(proto, transformed_vertex_time) {
  A <- proto$tmb_data$A_st
  spatial_idx <- as.integer(proto$tmb_data$A_spatial_index) + 1L
  year_idx <- as.integer(proto$tmb_data$year_i) + 1L
  n_t <- ncol(transformed_vertex_time)

  projected_by_t <- lapply(seq_len(n_t), function(tt) {
    as.numeric(A %*% transformed_vertex_time[, tt])
  })

  vapply(seq_along(year_idx), function(i) {
    projected_by_t[[year_idx[i]]][spatial_idx[i]]
  }, numeric(1))
}

# Build the minimal data list required by the template objective.
# We set n_i = 0 on purpose so we exercise only the lag-transform path (`T_s2t`)
# and avoid differences from unrelated observation-likelihood plumbing.
build_template_data <- function(proto) {
  n_t <- proto$tmb_data$n_t
  n_s <- ncol(proto$tmb_data$A_st)
  zero_A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, n_s))

  list(
    options_z = c(0L, 0L),
    n_t = n_t,
    a_g = 1,
    b_i = numeric(0),
    a_i = numeric(0),
    t_i = integer(0),
    Temp_s2t = proto$distributed_lags_data$covariate_vertex_time[, , 1, drop = TRUE],
    X_ik = matrix(0, nrow = 0, ncol = 1),
    A_is = zero_A,
    A_is2 = zero_A,
    M0_ss = proto$tmb_data$spde$M0,
    M1_ss = proto$tmb_data$spde$M1,
    M2_ss = proto$tmb_data$spde$M2,
    M1_s2s2 = proto$tmb_data$spde$M1,
    invM0_s2s2 = solve(proto$tmb_data$spde$M0),
    M_tt = make_shift_time_matrix(n_t),
    sim_gmrf = 0L
  )
}

# Compare one lag component between:
# 1) current sdmTMB distributed-lag implementation
# 2) legacy template C++ transform + projection done in R
compare_component <- function(component, kappaS = 1.5, kappaT = 0.4, kappaST = -0.3) {
  dat <- prepare_dataset()
  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 20)

  dl_formula <- switch(component,
    spatial = ~ spatial(depth_scaled),
    temporal = ~ temporal(depth_scaled),
    spatiotemporal = ~ spatiotemporal(depth_scaled),
    stop("Unknown component: ", component)
  )

  proto <- sdmTMB::sdmTMB(
    y_gaussian ~ 0,
    data = dat,
    mesh = mesh,
    time = "year",
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = dl_formula,
    do_fit = FALSE
  )

  # Isolate lag contribution in sdmTMB by setting:
  # - all fixed effects to 0
  # - lag coefficient to 1
  lag_idx <- match(proto$distributed_lags_data$term_coef_name, colnames(proto$tmb_data$X_ij[[1]]))
  proto$tmb_params$b_j[] <- 0
  proto$tmb_params$b_j[lag_idx] <- 1

  # Match the lag-scale parameters across both implementations.
  if (length(proto$tmb_params$log_kappaS_dl) > 0) {
    proto$tmb_params$log_kappaS_dl[] <- log(kappaS)
  }
  if (length(proto$tmb_params$log_kappaT_dl) > 0) {
    proto$tmb_params$log_kappaT_dl[] <- log(kappaT)
  }
  if (length(proto$tmb_params$kappaST_dl_unscaled) > 0) {
    proto$tmb_params$kappaST_dl_unscaled[] <- qlogis(-kappaST)
  }

  obj_sdm <- TMB::MakeADFun(
    data = proto$tmb_data,
    parameters = proto$tmb_params,
    map = proto$tmb_map,
    DLL = "sdmTMB",
    silent = TRUE
  )
  eta_sdm <- as.numeric(obj_sdm$report()$eta_i[, 1])

  # Run the template on the same vertex-time covariate input.
  tpl_data <- build_template_data(proto)
  n_s <- ncol(proto$tmb_data$A_st)
  n_t <- proto$tmb_data$n_t

  tpl_par <- list(
    log_kappaS = if (component %in% c("spatial", "spatiotemporal")) log(kappaS) else numeric(0),
    kappaT = if (identical(component, "temporal")) kappaT else numeric(0),
    kappaST = if (identical(component, "spatiotemporal")) kappaST else numeric(0),
    gamma_j = 1,
    beta_k = 0,
    ln_tauO = log(10),
    ln_tauE = log(10),
    ln_kappa = log(sqrt(8) / 100),
    ln_phi = log(1),
    logit_rhoE = 0,
    finv_power = 0,
    omega_s = rep(0, n_s),
    epsilon_st = matrix(0, nrow = n_s, ncol = n_t)
  )

  # `T_s2t` is the transformed lag field on vertices x time.
  obj_tpl <- TMB::MakeADFun(
    data = tpl_data,
    parameters = tpl_par,
    random = c("omega_s", "epsilon_st"),
    DLL = basename(template_dll),
    silent = TRUE
  )
  transformed <- obj_tpl$report()$T_s2t
  term_template <- project_template_to_observations(proto, transformed)

  # Compare sdmTMB lag term vs template lag term at observation level.
  diff <- eta_sdm - term_template
  data.frame(
    dataset = "pcod_2011",
    component = component,
    max_abs_diff = max(abs(diff)),
    mean_abs_diff = mean(abs(diff)),
    stringsAsFactors = FALSE
  )
}

# ----------------------------- Run checks ------------------------------------
cat("== Real-data fitting check ==\n")
fit_result <- fit_real_model()
cat(
  sprintf(
    "%s: n=%d years=%d objective=%.4f %s=%.6f\n",
    fit_result$dataset,
    fit_result$n,
    fit_result$years,
    fit_result$objective,
    fit_result$lag_name,
    fit_result$lag_beta
  )
)

cat("\n== Template sanity comparison (pcod_2011) ==\n")
cmp <- do.call(rbind, lapply(c("spatial", "temporal", "spatiotemporal"), function(comp) {
  compare_component(comp)
}))
print(cmp, row.names = FALSE)

# Component-specific tolerances:
# - spatial/temporal are expected to match at machine precision.
# - spatiotemporal currently allows a small discrepancy.
spatial_ok <- cmp$max_abs_diff[cmp$component == "spatial"] < 1e-8
temporal_ok <- cmp$max_abs_diff[cmp$component == "temporal"] < 1e-8
spatiotemporal_ok <- cmp$max_abs_diff[cmp$component == "spatiotemporal"] < 5e-3

if (!all(spatial_ok, temporal_ok, spatiotemporal_ok)) {
  stop("Template sanity thresholds failed; inspect printed differences.")
}

cat("\nSanity thresholds passed.\n")
