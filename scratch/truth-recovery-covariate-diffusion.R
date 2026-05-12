#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sdmTMB)
})

parse_args <- function(args) {
  defaults <- list(
    n_reps = 10L,
    n_obs = 1500L,
    n_t = 8L,
    cutoff = 0.06,
    n_hotspots = 60L,
    hotspot_sd = 0.035,
    seed = 20260512L,
    noise_sd = 0.15,
    kappaS_true = 6,
    beta_space_true = 0.8,
    beta_time_true = 0.7,
    beta_st_true = 0.6,
    intercept_space_true = 0.5,
    intercept_time_true = -0.2,
    intercept_st_true = 0.1,
    newton_loops = 0L
  )

  if (length(args) == 0L) return(defaults)

  out <- defaults
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L) next
    k <- kv[1]
    v <- kv[2]
    if (!k %in% names(out)) next
    if (is.integer(out[[k]])) {
      out[[k]] <- as.integer(v)
    } else {
      out[[k]] <- as.numeric(v)
    }
  }
  out
}

simulate_one <- function(cfg, rep_id) {
  set.seed(cfg$seed + rep_id)

  dat <- data.frame(
    X = runif(cfg$n_obs),
    Y = runif(cfg$n_obs),
    year = sample(seq_len(cfg$n_t), size = cfg$n_obs, replace = TRUE)
  )

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cfg$cutoff)

  centers <- cbind(runif(cfg$n_hotspots), runif(cfg$n_hotspots))
  weights <- rnorm(cfg$n_hotspots)
  sq_dist <- outer(dat$X, centers[, 1], "-")^2 +
    outer(dat$Y, centers[, 2], "-")^2
  hotspot_field <- as.numeric(exp(-sq_dist / (2 * cfg$hotspot_sd^2)) %*% weights)

  A_st <- mesh$A_st
  M0 <- mesh$spde$c0
  M1 <- mesh$spde$g1

  x_fine <- as.numeric(scale(hotspot_field + rnorm(cfg$n_obs, sd = 0.2)))
  numerator <- as.vector(Matrix::crossprod(A_st, x_fine))
  denominator <- as.vector(Matrix::crossprod(A_st, rep(1, cfg$n_obs)))
  x_vertex <- numeric(ncol(A_st))
  keep <- denominator > 0
  x_vertex[keep] <- numerator[keep] / denominator[keep]

  diffusion_matrix <- M0 + (1 / cfg$kappaS_true^2) * M1
  x_vertex_diffused <- as.numeric(Matrix::solve(diffusion_matrix, M0 %*% x_vertex))
  x_space_truth <- as.numeric(A_st %*% x_vertex_diffused)

  t_scaled <- as.numeric(scale(dat$year))
  x_time <- as.numeric(scale(sin(dat$year / 2) + 0.4 * dat$year + rnorm(cfg$n_obs, sd = 0.25)))
  x_st <- as.numeric(scale(hotspot_field * t_scaled + rnorm(cfg$n_obs, sd = 0.2)))

  dat$x_space <- x_fine
  dat$x_time <- x_time
  dat$x_st <- x_st

  dat$y_space <- cfg$intercept_space_true + cfg$beta_space_true * x_space_truth +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)
  dat$y_time <- cfg$intercept_time_true + cfg$beta_time_true * dat$x_time +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)
  dat$y_st <- cfg$intercept_st_true + cfg$beta_st_true * dat$x_st +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)

  list(data = dat, mesh = mesh)
}

fit_extract <- function(dat, mesh, cfg) {
  ctrl <- sdmTMBcontrol(newton_loops = cfg$newton_loops, getsd = FALSE)

  fit_space <- sdmTMB(
    y_space ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ space(x_space),
    control = ctrl
  )

  fit_time <- sdmTMB(
    y_time ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ time(x_time),
    control = ctrl
  )

  fit_st <- sdmTMB(
    y_st ~ 1,
    data = dat,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = ~ spacetime(x_st),
    control = ctrl
  )

  get_beta <- function(fit) {
    b <- fit$tmb_obj$env$parList()$b_j
    names(b) <- colnames(fit$tmb_data$X_ij[[1]])
    b
  }

  b_space <- get_beta(fit_space)
  b_time <- get_beta(fit_time)
  b_st <- get_beta(fit_st)

  rep_space <- fit_space$tmb_obj$report()
  rep_time <- fit_time$tmb_obj$report()
  rep_st <- fit_st$tmb_obj$report()

  c(
    ll_space = as.numeric(logLik(fit_space)),
    ll_time = as.numeric(logLik(fit_time)),
    ll_st = as.numeric(logLik(fit_st)),
    intercept_space = b_space[["(Intercept)"]],
    beta_space = b_space[["dl_space_x_space"]],
    kappaS_dl = rep_space$kappaS_dl[1],
    intercept_time = b_time[["(Intercept)"]],
    beta_time = b_time[["dl_time_x_time"]],
    kappaT_dl = rep_time$kappaT_dl[1],
    rhoT = rep_time$rhoT[1],
    intercept_st = b_st[["(Intercept)"]],
    beta_st = b_st[["dl_spacetime_x_st"]],
    kappaST_dl = rep_st$kappaST_dl[1]
  )
}

summarize_recovery <- function(res, cfg) {
  truth <- c(
    intercept_space = cfg$intercept_space_true,
    beta_space = cfg$beta_space_true,
    kappaS_dl = cfg$kappaS_true,
    intercept_time = cfg$intercept_time_true,
    beta_time = cfg$beta_time_true,
    intercept_st = cfg$intercept_st_true,
    beta_st = cfg$beta_st_true
  )

  rows <- lapply(names(truth), function(nm) {
    x <- res[[nm]]
    tr <- truth[[nm]]
    c(
      parameter = nm,
      truth = tr,
      mean_est = mean(x),
      median_est = stats::median(x),
      sd_est = stats::sd(x),
      bias = mean(x) - tr,
      rel_bias_pct = 100 * (mean(x) - tr) / ifelse(abs(tr) < 1e-12, NA_real_, tr),
      rmse = sqrt(mean((x - tr)^2))
    )
  })

  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  num_cols <- setdiff(names(out), "parameter")
  out[num_cols] <- lapply(out[num_cols], as.numeric)
  out
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))

  message("Running truth-recovery experiment with config:")
  print(cfg)

  res <- vector("list", cfg$n_reps)
  failed <- integer(0)

  for (i in seq_len(cfg$n_reps)) {
    message(sprintf("rep %d/%d", i, cfg$n_reps))
    sim <- simulate_one(cfg, i)
    fit_res <- try(fit_extract(sim$data, sim$mesh, cfg), silent = TRUE)
    if (inherits(fit_res, "try-error")) {
      failed <- c(failed, i)
      next
    }
    res[[i]] <- fit_res
  }

  ok <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(ok)) {
    stop("All replicates failed.")
  }

  res_df <- as.data.frame(do.call(rbind, res[ok]))
  res_df[] <- lapply(res_df, as.numeric)

  summary_df <- summarize_recovery(res_df, cfg)

  cat("\n=== Truth recovery summary ===\n")
  print(summary_df, row.names = FALSE, digits = 4)

  cat("\n=== Diagnostics ===\n")
  cat("successful reps:", sum(ok), "/", cfg$n_reps, "\n")
  cat("failed reps:", length(failed), "\n")
  if (length(failed) > 0L) {
    cat("failed rep ids:", paste(failed, collapse = ", "), "\n")
  }

  if ("kappaS_dl" %in% names(res_df)) {
    cat(sprintf("kappaS_dl mean (truth %.3f): %.4f\n", cfg$kappaS_true, mean(res_df$kappaS_dl)))
  }
  if ("kappaT_dl" %in% names(res_df)) {
    cat(sprintf("kappaT_dl mean (expected near 0): %.4g\n", mean(res_df$kappaT_dl)))
  }
  if ("rhoT" %in% names(res_df)) {
    cat(sprintf("rhoT mean (expected near 0): %.4g\n", mean(res_df$rhoT)))
  }
}

main()
