suppressPackageStartupMessages({
  library(pkgload)
  library(TMB)
  library(Matrix)
})

# Use local source tree (current branch state).
pkgload::load_all(".", quiet = TRUE)

# -----------------------------------------------------------------------------
# Data + mesh setup
# -----------------------------------------------------------------------------
prepare_pcod <- function() {
  dat <- as.data.frame(sdmTMB::pcod_2011)
  keep <- is.finite(dat$density) & is.finite(dat$depth) & !is.na(dat$year)
  dat <- dat[keep, , drop = FALSE]
  dat$depth_scaled <- as.numeric(scale(dat$depth))
  dat
}

make_pcod_mesh <- function(dat) {
  sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 20)
}

component_formula <- function(component) {
  switch(component,
    spatial = ~ spatial(depth_scaled),
    temporal = ~ temporal(depth_scaled),
    stop("Unsupported component: ", component)
  )
}

# -----------------------------------------------------------------------------
# sdmTMB optimization (full fit)
# -----------------------------------------------------------------------------
fit_sdmTMB_model <- function(dat, mesh, component) {
  sdmTMB::sdmTMB(
    density ~ 0,
    data = dat,
    mesh = mesh,
    time = "year",
    family = tweedie(link = "log"),
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = component_formula(component),
    control = sdmTMB::sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
  )
}

# -----------------------------------------------------------------------------
# Template optimization (legacy spacetime_lag_2025_10_24.cpp)
# -----------------------------------------------------------------------------
load_template_dll <- function() {
  cpp <- "scratch/spacetime_lag_2025_10_24.cpp"
  dll <- sub("\\.cpp$", "", cpp)
  so <- TMB::dynlib(dll)

  if (!file.exists(so)) {
    TMB::compile(cpp)
  }
  dll_name <- basename(dll)
  if (!(dll_name %in% names(getLoadedDLLs()))) {
    dyn.load(so)
  }
  invisible(dll_name)
}

# Build a proto object just to reuse sdmTMB lag preprocessing/projection internals.
make_proto <- function(dat, mesh, component) {
  sdmTMB::sdmTMB(
    density ~ 0,
    data = dat,
    mesh = mesh,
    time = "year",
    # Gaussian here is just for preprocessing/plumbing; this object is not fitted.
    family = gaussian(),
    spatial = "off",
    spatiotemporal = "off",
    distributed_lags = component_formula(component),
    do_fit = FALSE
  )
}

make_template_data <- function(proto) {
  n_t <- proto$tmb_data$n_t
  # In current sdmTMB plumbing `tmb_data$n_i` can be 0 for `do_fit = FALSE`,
  # so derive observation count from available row-level objects instead.
  n_i <- length(proto$tmb_data$year_i)
  # Subset-copy to materialize a plain sparse matrix object; passing the original
  # object directly has triggered Eigen bounds errors in this template.
  A_obs <- proto$tmb_data$A_st[seq_len(n_i), , drop = FALSE]
  M_tt <- Matrix::sparseMatrix(
    i = seq.int(2L, n_t),
    j = seq.int(1L, n_t - 1L),
    x = 1,
    dims = c(n_t, n_t)
  )

  list(
    options_z = c(0L, 0L),
    n_t = n_t,
    a_g = 1,
    b_i = pmax(proto$data$density, 0),
    a_i = rep(1, n_i),
    t_i = as.integer(proto$tmb_data$year_i), # 0-based in this template
    Temp_s2t = proto$distributed_lags_data$covariate_vertex_time[, , 1, drop = TRUE],
    # Template expects beta_i to have length n_i; keep a zero column and fix beta_k = 0.
    X_ik = matrix(0, n_i, 1),
    A_is = A_obs,
    A_is2 = A_obs,
    M0_ss = proto$tmb_data$spde$M0,
    M1_ss = proto$tmb_data$spde$M1,
    M2_ss = proto$tmb_data$spde$M2,
    M1_s2s2 = proto$tmb_data$spde$M1,
    invM0_s2s2 = solve(proto$tmb_data$spde$M0),
    M_tt = M_tt,
    sim_gmrf = 0L
  )
}

make_template_parameters <- function(proto, component) {
  n_s <- ncol(proto$tmb_data$A_st)
  n_t <- proto$tmb_data$n_t

  list(
    log_kappaS = if (identical(component, "spatial")) log(1.5) else numeric(0),
    kappaT = if (identical(component, "temporal")) 0.4 else numeric(0),
    kappaST = numeric(0),
    gamma_j = 0,
    beta_k = 0,
    # Keep these fixed to avoid fitting unrelated random-field structure.
    ln_tauO = log(10),
    ln_tauE = log(10),
    ln_kappa = log(sqrt(8) / 100),
    ln_phi = log(0.2),
    logit_rhoE = 0,
    finv_power = 0,
    omega_s = rep(0, n_s),
    epsilon_st = matrix(0, nrow = n_s, ncol = n_t)
  )
}

make_template_map <- function(proto) {
  n_s <- ncol(proto$tmb_data$A_st)
  n_t <- proto$tmb_data$n_t
  list(
    beta_k = factor(NA),
    ln_tauO = factor(NA),
    ln_tauE = factor(NA),
    ln_kappa = factor(NA),
    logit_rhoE = factor(NA),
    omega_s = factor(rep(NA, n_s)),
    epsilon_st = factor(rep(NA, n_s * n_t))
  )
}

fit_template_model <- function(proto, dll_name, component) {
  data_tpl <- make_template_data(proto)
  par_tpl <- make_template_parameters(proto, component)
  map_tpl <- make_template_map(proto)

  obj <- TMB::MakeADFun(
    data = data_tpl,
    parameters = par_tpl,
    map = map_tpl,
    DLL = dll_name,
    silent = TRUE
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = 300, eval.max = 400)
  )

  list(obj = obj, opt = opt, report = obj$report(opt$par))
}

# -----------------------------------------------------------------------------
# Extraction + comparison
# -----------------------------------------------------------------------------
extract_sdmTMB_summary <- function(fit, component) {
  par <- fit$tmb_obj$env$parList(fit$model$par)
  lag_idx <- match(fit$distributed_lags_data$term_coef_name, colnames(fit$tmb_data$X_ij[[1]]))

  kappa_name <- if (identical(component, "spatial")) "kappaS" else "kappaT"
  kappa_value <- if (identical(component, "spatial")) {
    exp(par$log_kappaS_dl[1])
  } else {
    exp(par$log_kappaT_dl[1])
  }

  list(
    objective = fit$model$objective,
    kappa_name = kappa_name,
    kappa_value = kappa_value,
    beta_lag = unname(par$b_j[lag_idx]),
    phi = exp(par$ln_phi[1]),
    tweedie_p = 1 + plogis(par$thetaf[1]),
    eta_lag = as.numeric(fit$tmb_obj$report()$eta_i[, 1])
  )
}

extract_template_summary <- function(fit_tpl, component) {
  par <- fit_tpl$opt$par

  kappa_name <- if (identical(component, "spatial")) "kappaS" else "kappaT"
  kappa_value <- if (identical(component, "spatial")) {
    exp(par["log_kappaS"])
  } else {
    par["kappaT"]
  }

  list(
    objective_total = fit_tpl$opt$objective,
    objective_data = fit_tpl$report$jnll_comp[1],
    objective_re = sum(fit_tpl$report$jnll_comp[2:3]),
    kappa_name = kappa_name,
    kappa_value = kappa_value,
    beta_lag = unname(par["gamma_j"]),
    phi = exp(par["ln_phi"]),
    tweedie_p = 1 + plogis(par["finv_power"]),
    eta_lag = as.numeric(fit_tpl$report$gamma_i)
  )
}

compare_effects <- function(eta_a, eta_b) {
  if (length(eta_a) != length(eta_b)) {
    stop("Lag-effect vectors have different lengths: ", length(eta_a), " vs ", length(eta_b))
  }
  data.frame(
    correlation = cor(eta_a, eta_b),
    rmse = sqrt(mean((eta_a - eta_b)^2)),
    mae = mean(abs(eta_a - eta_b)),
    max_abs_diff = max(abs(eta_a - eta_b))
  )
}

run_component <- function(component, dat, mesh, dll_name) {
  cat(sprintf("\n== Component: %s ==\n", component))

  cat("Optimize sdmTMB...\n")
  fit_sdm <- fit_sdmTMB_model(dat, mesh, component)
  sdm_sum <- extract_sdmTMB_summary(fit_sdm, component)

  cat("Optimize template...\n")
  proto <- make_proto(dat, mesh, component)
  fit_tpl <- fit_template_model(proto, dll_name, component)
  tpl_sum <- extract_template_summary(fit_tpl, component)

  param_cmp <- data.frame(
    component = component,
    parameter = c("objective_data", sdm_sum$kappa_name, "lag_beta", "phi", "tweedie_p"),
    sdmTMB = c(sdm_sum$objective, sdm_sum$kappa_value, sdm_sum$beta_lag, sdm_sum$phi, sdm_sum$tweedie_p),
    template = c(tpl_sum$objective_data, tpl_sum$kappa_value, tpl_sum$beta_lag, tpl_sum$phi, tpl_sum$tweedie_p)
  )

  eff_cmp <- cbind(
    data.frame(component = component, stringsAsFactors = FALSE),
    compare_effects(sdm_sum$eta_lag, tpl_sum$eta_lag)
  )

  print(param_cmp, row.names = FALSE)
  cat(
    sprintf(
      "template objective breakdown: total=%.6f data=%.6f re_const=%.6f\n",
      tpl_sum$objective_total,
      tpl_sum$objective_data,
      tpl_sum$objective_re
    )
  )
  print(eff_cmp, row.names = FALSE)

  list(param_cmp = param_cmp, eff_cmp = eff_cmp)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
cat("== Build data ==\n")
dat <- prepare_pcod()
mesh <- make_pcod_mesh(dat)
cat(sprintf("pcod_2011 rows=%d years=%d mesh_vertices=%d\n", nrow(dat), length(unique(dat$year)), mesh$mesh$n))

dll_name <- load_template_dll()
components <- c("spatial", "temporal")
results <- lapply(components, function(comp) run_component(comp, dat, mesh, dll_name))

cat("\n== Summary: lag-effect diagnostics ==\n")
summary_eff <- do.call(rbind, lapply(results, `[[`, "eff_cmp"))
print(summary_eff, row.names = FALSE)

cat("\nDone.\n")
