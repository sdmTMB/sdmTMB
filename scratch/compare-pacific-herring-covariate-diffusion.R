#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TMB)
  library(surveyjoin)
  library(sf)
  library(fmesher)
  library(Matrix)
})

# Compare distributed-lag parameter estimates between legacy TMB and sdmTMB
# for Pacific herring with matched preprocessing and model structure.

root_repo <- normalizePath(".", winslash = "/", mustWork = TRUE)
pkgload::load_all(root_repo, quiet = TRUE)
legacy_dir <- normalizePath("~/src/spacetime-lag", winslash = "/", mustWork = TRUE)
legacy_cpp <- file.path(legacy_dir, "spacetime_lag_2025_10_24.cpp")
out_dir <- file.path(root_repo, "scratch")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

species <- "pacific herring"
rescale <- 1
cutoff <- as.numeric(Sys.getenv("SDMTMB_COMPARE_CUTOFF", "180"))
year_set <- NULL
covariate_mode <- "anomaly"
rel_tol <- 0.01
abs_tol <- 1e-3

message("Preparing data for ", species)

d <- surveyjoin::get_data(species)
d <- as.data.frame(d)
d <- subset(d, survey_name == "eastern Bering Sea" & !is.na(bottom_temp_c))

if (!nrow(d)) stop("No rows after filtering EBS + bottom_temp_c non-missing.")
year_obs <- sort(unique(d$year))
year_set <- seq(min(year_obs) - 1L, max(year_obs))

sf_pts <- sf::st_as_sf(d, coords = c("lon_start", "lat_start"), crs = "+proj=longlat +datum=WGS84")
sf_pts <- sf::st_transform(sf_pts, crs = "+proj=utm +zone=2 +datum=WGS84 +units=km")

coords <- sf::st_coordinates(sf_pts) / rescale
mesh <- fmesher::fm_mesh_2d(coords, cutoff = cutoff / rescale, refine = TRUE)
message("Legacy/shared mesh vertices: ", mesh$n)

spde <- fmesher::fm_fem(mesh, order = 2)

# Legacy-like vertex-year covariate field using nearest-neighbor fill.
Temp_s2t <- array(0, dim = c(mesh$n, length(year_set)), dimnames = list(NULL, year = year_set))
for (ti in seq_along(year_set)) {
  yr <- year_set[ti]
  if (yr %in% sf_pts$year) {
    sf_t <- sf_pts[sf_pts$year == yr, ]
    nn <- RANN::nn2(data = sf::st_coordinates(sf_t) / rescale, query = mesh$loc[, 1:2], k = 1)$nn.idx
    Temp_s2t[, ti] <- sf_t$bottom_temp_c[nn]
  }
}

if ((min(year_obs) - 1L) %in% year_set && all((min(year_obs) + 0:1) %in% year_set)) {
  y0 <- as.character(min(year_obs) - 1L)
  y1 <- as.character(min(year_obs))
  y2 <- as.character(min(year_obs) + 1L)
  Temp_s2t[, y0] <- rowMeans(Temp_s2t[, c(y1, y2), drop = FALSE])
}
missing_interior_years <- setdiff(year_set, year_obs)
missing_interior_years <- missing_interior_years[
  (missing_interior_years - 1L) %in% year_set & (missing_interior_years + 1L) %in% year_set
]
for (yr in missing_interior_years) {
  Temp_s2t[, as.character(yr)] <- rowMeans(
    Temp_s2t[, as.character(c(yr - 1L, yr + 1L)), drop = FALSE]
  )
}

if (identical(covariate_mode, "anomaly")) {
  Temp_s2t <- sweep(Temp_s2t, MARGIN = 1, STATS = rowMeans(Temp_s2t), FUN = "-")
}

analysis_df <- cbind(as.data.frame(sf::st_coordinates(sf_pts)),
  year = sf_pts$year,
  catch_weight = sf_pts$catch_weight,
  effort_km2 = sf_pts$effort / 100,
  lag_x = sf_pts$bottom_temp_c
)

analysis_df$year_num <- as.integer(analysis_df$year)
analysis_df$year_i0 <- match(analysis_df$year_num, year_set) - 1L
analysis_df <- analysis_df[order(analysis_df$year, seq_len(nrow(analysis_df))), ]
rownames(analysis_df) <- NULL

# Build projection after sorting, so A rows match response/offset/year rows.
A_is <- fmesher::fm_evaluator(mesh, loc = as.matrix(analysis_df[, c("X", "Y")]))$proj$A
T_it <- as.matrix(A_is %*% Temp_s2t)
t_i <- match(analysis_df$year_num, year_set)
if (anyNA(t_i)) stop("Observed year outside year_set.")
analysis_df$lag_x <- as.numeric(T_it[cbind(seq_len(nrow(T_it)), t_i)])

sdm_extra_time <- setdiff(year_set, sort(unique(analysis_df$year_num)))
Temp_s2t_sdm <- Temp_s2t
shared_sdm_mesh <- sdmTMB::make_mesh(
  analysis_df,
  xy_cols = c("X", "Y"),
  mesh = mesh
)
if (!isTRUE(all.equal(shared_sdm_mesh$mesh$loc, mesh$loc))) {
  stop("sdmTMB mesh vertices differ from legacy mesh vertices.")
}

M_tt <- Matrix::sparseMatrix(
  i = 2:length(year_set),
  j = 1:(length(year_set) - 1),
  x = 1,
  dims = c(length(year_set), length(year_set))
)

build_legacy_data <- function(component) {
  data_tmb <- list(
    options_z = c(1L, 0L), # stationary version, linear lag
    n_t = length(year_set),
    a_g = rep(1, 1),
    b_i = analysis_df$catch_weight,
    a_i = analysis_df$effort_km2,
    t_i = analysis_df$year_i0,
    Temp_s2t = Temp_s2t,
    X_ik = model.matrix(~ 1, data = analysis_df),
    A_is = A_is,
    A_is2 = A_is,
    A_gs = Matrix::Matrix(0, 1, 1, sparse = TRUE),
    A_gs2 = Matrix::Matrix(0, 1, 1, sparse = TRUE),
    M0_ss = spde$c0,
    M1_ss = spde$g1,
    M2_ss = spde$g2,
    invM0_ss = solve(spde$c0),
    M1_s2s2 = spde$g1,
    invM0_s2s2 = solve(spde$c0),
    M_tt = M_tt,
    sim_gmrf = 0L
  )

  par <- list(
    log_kappaS = log(sqrt(4 / 20)),
    kappaT = 0,
    kappaST = 0,
    gamma_j = 0.1,
    beta_k = 0,
    ln_tauO = log(10) - log(rescale),
    ln_tauE = log(10) - log(rescale),
    ln_kappa = log(sqrt(8) / (100 / rescale)),
    ln_phi = log(1),
    logit_rhoE = 0,
    finv_power = 0,
    omega_s = rep(0, mesh$n),
    epsilon_st = matrix(0, nrow = mesh$n, ncol = length(year_set))
  )

  if (identical(component, "space")) {
    par$kappaT <- numeric(0)
    par$kappaST <- numeric(0)
  }
  if (identical(component, "time")) {
    par$log_kappaS <- numeric(0)
    par$kappaST <- numeric(0)
  }
  if (identical(component, "spacetime")) {
    par$kappaT <- numeric(0)
  }

  list(data = data_tmb, par = par)
}

fit_legacy <- function(component) {
  dat_par <- build_legacy_data(component)
  obj <- TMB::MakeADFun(
    data = dat_par$data,
    parameters = dat_par$par,
    random = c("omega_s", "epsilon_st"),
    silent = TRUE,
    DLL = tools::file_path_sans_ext(basename(legacy_cpp))
  )

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(eval.max = 2000, iter.max = 2000)
  )
  rep <- obj$report()
  pl <- obj$env$parList(opt$par)

  out <- c(
    gamma = as.numeric(pl$gamma_j[1]),
    kappaS = NA_real_,
    kappaT = NA_real_,
    kappaST = NA_real_,
    rhoT = NA_real_,
    obj = opt$objective,
    convergence = opt$convergence
  )

  if ("log_kappaS" %in% names(pl)) out["kappaS"] <- exp(pl$log_kappaS[1])
  if ("kappaT" %in% names(pl)) out["kappaT"] <- pl$kappaT[1]
  if ("kappaST" %in% names(pl)) out["kappaST"] <- pl$kappaST[1]
  if (!is.null(rep$rhoT)) out["rhoT"] <- rep$rhoT[1]
  out
}

fit_sdmTMB <- function(component) {
  cd_form <- switch(component,
    space = ~ space(lag_x),
    time = ~ time(lag_x),
    spacetime = ~ spacetime(lag_x)
  )
  dl_start <- switch(component,
    space = list(log_kappaS_dl = log(sqrt(4 / 20))),
    time = list(log_kappaT_dl = 0),
    spacetime = list(
      log_kappaS_dl = log(sqrt(4 / 20)),
      kappaST_dl_unscaled = 0
    )
  )

  fit <- sdmTMB::sdmTMB(
    catch_weight ~ 1,
    data = analysis_df,
    mesh = shared_sdm_mesh,
    time = "year_num",
    extra_time = sdm_extra_time,
    offset = log(analysis_df$effort_km2),
    family = tweedie(link = "log"),
    spatial = "on",
    spatiotemporal = "ar1",
    covariate_diffusion = cd_form,
    control = sdmTMB::sdmTMBcontrol(
      newton_loops = 0,
      eval.max = 2000,
      iter.max = 2000,
      multiphase = FALSE,
      getsd = FALSE,
      start = c(list(
        b_j = c(0, 0.1),
        ln_tau_O = log(10) - log(rescale),
        ln_tau_E = log(10) - log(rescale),
        ln_kappa = matrix(log(sqrt(8) / (100 / rescale)), nrow = 2L, ncol = 1L),
        ln_phi = log(1),
        thetaf = 0,
        ar1_phi = qlogis((0.5 + 1) / 2)
      ), dl_start)
    ),
    experimental = list(
      covariate_diffusion_covariate_vertex = Temp_s2t_sdm
    ),
    silent = TRUE
  )

  pl <- fit$tmb_obj$env$parList(fit$tmb_obj$env$last.par.best[fit$tmb_obj$env$lfixed()])
  rep <- fit$tmb_obj$report(fit$tmb_obj$env$last.par.best)
  beta_names <- colnames(fit$tmb_data$X_ij[[1]])
  b <- pl$b_j
  names(b) <- beta_names
  dl_name <- grep("^dl_", names(b), value = TRUE)

  c(
    gamma = if (length(dl_name)) as.numeric(b[dl_name[1]]) else NA_real_,
    kappaS = if (!is.null(rep$kappaS_dl) && length(rep$kappaS_dl)) rep$kappaS_dl[1] else NA_real_,
    kappaT = if (!is.null(rep$kappaT_dl) && length(rep$kappaT_dl)) rep$kappaT_dl[1] else NA_real_,
    kappaST = if (!is.null(rep$kappaST_dl) && length(rep$kappaST_dl)) rep$kappaST_dl[1] else NA_real_,
    rhoT = if (!is.null(rep$rhoT) && length(rep$rhoT)) rep$rhoT[1] else NA_real_,
    obj = fit$gradients$objective,
    convergence = fit$convergence
  )
}

# Compile/load legacy DLL
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(legacy_dir)
legacy_dll <- tools::file_path_sans_ext(basename(legacy_cpp))
if (!file.exists(TMB::dynlib(legacy_dll))) {
  TMB::compile(basename(legacy_cpp), framework = "TMBad")
}
dyn.load(TMB::dynlib(legacy_dll))

components <- c("space", "time", "spacetime")
rows <- list()
for (comp in components) {
  message("Fitting component: ", comp)
  leg <- fit_legacy(comp)
  sdm <- fit_sdmTMB(comp)

  pars <- c("gamma", "kappaS", "kappaT", "kappaST", "rhoT")
  for (p in pars) {
    lv <- unname(leg[p])
    sv <- unname(sdm[p])
    abs_diff <- if (is.na(lv) || is.na(sv)) NA_real_ else abs(sv - lv)
    rel <- if (is.na(lv) || is.na(sv) || abs(lv) < 1e-12) NA_real_ else abs_diff / abs(lv)
    rows[[length(rows) + 1L]] <- data.frame(
      component = comp,
      parameter = p,
      legacy = lv,
      sdmTMB = sv,
      abs_diff = abs_diff,
      rel_diff = rel,
      pass = if (is.na(abs_diff)) NA else abs_diff < abs_tol || (!is.na(rel) && rel < rel_tol)
    )
  }
}

out <- do.call(rbind, rows)
out_file <- file.path(out_dir, "pacific-herring-covariate-diffusion-comparison.csv")
write.csv(out, out_file, row.names = FALSE)

message("Saved: ", out_file)
print(out)
