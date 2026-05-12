#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sdmTMB)
})

# User-editable defaults:
# Set values directly here (or leave NULL to keep defaults).
USER_SETTINGS <- list(
  n_reps = 10,
  n_obs = 3000
)

parse_args <- function(args) {
  defaults <- list(
    n_reps = 10L,
    n_obs = 1500L,
    n_t = 14L,
    cutoff = 0.06,
    seed = 20260512L,
    noise_sd = 0.08,
    newton_loops = 0L,

    intercept_space_true = 0.5,
    beta_space_true = 0.8,
    kappaS_space_true = 6,

    intercept_time_true = -0.2,
    beta_time_true = 0.7,
    kappaT_time_true = 0.35,

    intercept_st_true = 0.1,
    beta_st_true = 0.8,
    kappaS_st_true = 40,
    kappaST_st_true = -0.1
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
    if (is.integer(out[[k]])) out[[k]] <- as.integer(v) else out[[k]] <- as.numeric(v)
  }
  out
}

apply_user_settings <- function(cfg, user_settings) {
  for (nm in names(user_settings)) {
    val <- user_settings[[nm]]
    if (is.null(val)) next
    if (!nm %in% names(cfg)) next
    cfg[[nm]] <- val
  }
  cfg
}

raw_covariates <- function(dat, seed_offset = 0L) {
  set.seed(seed_offset)
  hotspot_centers <- cbind(runif(50), runif(50))
  hotspot_weights <- rnorm(50)
  sq_dist <- outer(dat$X, hotspot_centers[, 1], "-")^2 +
    outer(dat$Y, hotspot_centers[, 2], "-")^2
  hotspot_field <- as.numeric(exp(-sq_dist / (2 * 0.04^2)) %*% hotspot_weights)

  t_scaled <- as.numeric(scale(dat$year))

  list(
    x_space = as.numeric(scale(hotspot_field + rnorm(nrow(dat), sd = 0.25))),
    x_time = as.numeric(scale(sin(dat$year / 2) + 0.4 * dat$year + rnorm(nrow(dat), sd = 0.25))),
    x_st = as.numeric(scale(hotspot_field * t_scaled + rnorm(nrow(dat), sd = 0.2)))
  )
}

make_proto <- function(dat, mesh, covariate_name, component) {
  d <- dat
  d$y <- 0
  form <- switch(
    component,
    space = stats::as.formula(paste0("~ space(", covariate_name, ")")),
    time = stats::as.formula(paste0("~ time(", covariate_name, ")")),
    spacetime = stats::as.formula(paste0("~ spacetime(", covariate_name, ")")),
    stop("Unknown component")
  )
  sdmTMB(
    y ~ 1,
    data = d,
    mesh = mesh,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    family = gaussian(),
    covariate_diffusion = form,
    do_fit = FALSE
  )
}

compute_term <- function(proto, kappaS = 1, kappaT = 0.1, kappaST = -0.1) {
  # Internal parameterization: kappaST = -plogis(kappaST_dl_unscaled)
  kappaST_unscaled <- stats::qlogis(-kappaST)

  term_mat <- sdmTMB:::.compute_covariate_diffusion_term_values(
    covariate_diffusion_data = proto[["covariate_diffusion_data"]],
    covariate_vertex_time = proto[["covariate_diffusion_data"]][["covariate_vertex_time"]],
    A_st = proto[["tmb_data"]][["A_st"]],
    A_spatial_index = proto[["tmb_data"]][["A_spatial_index"]],
    year_i = proto[["tmb_data"]][["year_i"]],
    n_t = proto[["tmb_data"]][["n_t"]],
    M0 = proto[["tmb_data"]][["spde"]][["M0"]],
    M1 = proto[["tmb_data"]][["spde"]][["M1"]],
    log_kappaS_dl = log(kappaS),
    log_kappaT_dl = log(kappaT),
    kappaST_dl_unscaled = kappaST_unscaled
  )
  as.numeric(term_mat[, 1])
}

simulate_one <- function(cfg, rep_id) {
  set.seed(cfg$seed + rep_id)

  n_sites <- max(10L, as.integer(cfg$n_obs / cfg$n_t))
  site <- seq_len(n_sites)
  site_xy <- data.frame(
    site = site,
    X = runif(n_sites),
    Y = runif(n_sites)
  )
  dat <- merge(
    expand.grid(site = site, year = seq_len(cfg$n_t)),
    site_xy,
    by = "site",
    sort = FALSE
  )
  dat <- dat[, c("X", "Y", "year", "site")]

  # Keep sample size close to cfg$n_obs while preserving repeated sites across time.
  if (nrow(dat) > cfg$n_obs) {
    dat <- dat[sample.int(nrow(dat), cfg$n_obs), , drop = FALSE]
  }
  if (nrow(dat) < cfg$n_obs) {
    add_n <- cfg$n_obs - nrow(dat)
    dat <- rbind(dat, dat[sample.int(nrow(dat), add_n, replace = TRUE), , drop = FALSE])
  }

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cfg$cutoff)

  covs <- raw_covariates(dat, seed_offset = cfg$seed + 1000L + rep_id)
  dat$x_space <- covs$x_space
  dat$x_time <- covs$x_time
  dat$x_st <- covs$x_st

  proto_space <- make_proto(dat, mesh, "x_space", "space")
  proto_time <- make_proto(dat, mesh, "x_time", "time")
  proto_st <- make_proto(dat, mesh, "x_st", "spacetime")

  term_space <- compute_term(proto_space, kappaS = cfg$kappaS_space_true, kappaT = 0.1, kappaST = -0.1)
  term_time <- compute_term(proto_time, kappaS = 1, kappaT = cfg$kappaT_time_true, kappaST = -0.1)
  term_st <- compute_term(proto_st, kappaS = cfg$kappaS_st_true, kappaT = 0.1, kappaST = cfg$kappaST_st_true)

  dat$y_space <- cfg$intercept_space_true + cfg$beta_space_true * term_space +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)
  dat$y_time <- cfg$intercept_time_true + cfg$beta_time_true * term_time +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)
  dat$y_st <- cfg$intercept_st_true + cfg$beta_st_true * term_st +
    rnorm(cfg$n_obs, sd = cfg$noise_sd)

  list(data = dat, mesh = mesh)
}

extract_beta <- function(fit) {
  b <- fit$tmb_obj$env$parList()$b_j
  names(b) <- colnames(fit$tmb_data$X_ij[[1]])
  b
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

  b_space <- extract_beta(fit_space)
  b_time <- extract_beta(fit_time)
  b_st <- extract_beta(fit_st)

  rep_space <- fit_space$tmb_obj$report()
  rep_time <- fit_time$tmb_obj$report()
  rep_st <- fit_st$tmb_obj$report()

  c(
    ll_space = as.numeric(logLik(fit_space)),
    ll_time = as.numeric(logLik(fit_time)),
    ll_st = as.numeric(logLik(fit_st)),

    intercept_space = b_space[["(Intercept)"]],
    beta_space = b_space[["dl_space_x_space"]],
    kappaS_space = rep_space$kappaS_dl[1],
    phi_space = rep_space$phi[1],

    intercept_time = b_time[["(Intercept)"]],
    beta_time = b_time[["dl_time_x_time"]],
    kappaT_time = rep_time$kappaT_dl[1],
    rhoT_time = rep_time$rhoT[1],
    phi_time = rep_time$phi[1],

    intercept_st = b_st[["(Intercept)"]],
    beta_st = b_st[["dl_spacetime_x_st"]],
    kappaS_st = rep_st$kappaS_dl[1],
    kappaST_st = rep_st$kappaST_dl[1],
    phi_st = rep_st$phi[1]
  )
}

summarize_recovery <- function(res, cfg) {
  truth <- c(
    intercept_space = cfg$intercept_space_true,
    beta_space = cfg$beta_space_true,
    kappaS_space = cfg$kappaS_space_true,
    phi_space = cfg$noise_sd,

    intercept_time = cfg$intercept_time_true,
    beta_time = cfg$beta_time_true,
    kappaT_time = cfg$kappaT_time_true,
    rhoT_time = cfg$kappaT_time_true / (1 + cfg$kappaT_time_true),
    phi_time = cfg$noise_sd,

    intercept_st = cfg$intercept_st_true,
    beta_st = cfg$beta_st_true,
    kappaS_st = cfg$kappaS_st_true,
    kappaST_st = cfg$kappaST_st_true,
    phi_st = cfg$noise_sd
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

plot_recovery <- function(res_df, cfg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plot.")
    return(invisible(NULL))
  }

  truth <- c(
    intercept_space = cfg$intercept_space_true,
    beta_space = cfg$beta_space_true,
    kappaS_space = cfg$kappaS_space_true,
    phi_space = cfg$noise_sd,
    intercept_time = cfg$intercept_time_true,
    beta_time = cfg$beta_time_true,
    kappaT_time = cfg$kappaT_time_true,
    rhoT_time = cfg$kappaT_time_true,
    phi_time = cfg$noise_sd,
    intercept_st = cfg$intercept_st_true,
    beta_st = cfg$beta_st_true,
    kappaS_st = cfg$kappaS_st_true,
    kappaST_st = cfg$kappaST_st_true,
    phi_st = cfg$noise_sd
  )

  keep <- intersect(names(truth), names(res_df))
  long <- do.call(rbind, lapply(keep, function(nm) {
    data.frame(
      parameter = nm,
      estimate = res_df[[nm]],
      truth = truth[[nm]],
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = parameter, y = estimate)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = truth), colour = "#d95f02", linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 1.5, colour = "#1b9e77") +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = "Estimate across replicates",
      title = "Covariate diffusion truth recovery",
      subtitle = sprintf("n_reps = %d, n_obs = %d", cfg$n_reps, cfg$n_obs)
    )

  print(p)
  invisible(p)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg <- apply_user_settings(cfg, USER_SETTINGS)

  message("Running truth-recovery experiment with config:")
  print(cfg)

  res <- vector("list", cfg$n_reps)
  failed <- integer(0)

  for (i in seq_len(cfg$n_reps)) {
    message(sprintf("rep %d/%d", i, cfg$n_reps))
    sim <- simulate_one(cfg, i)
    fit_res <- try(fit_extract(sim$data, sim$mesh, cfg), silent = TRUE)
    if (inherits(fit_res, "try-error")) {
      message("rep failed with error: ", as.character(fit_res))
      failed <- c(failed, i)
      next
    }
    res[[i]] <- fit_res
  }

  ok <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(ok)) stop("All replicates failed.")

  res_df <- as.data.frame(do.call(rbind, res[ok]))
  res_df[] <- lapply(res_df, as.numeric)

  summary_df <- summarize_recovery(res_df, cfg)

  cat("\n=== Truth recovery summary ===\n")
  print(summary_df, row.names = FALSE, digits = 4)

  cat("\n=== Diagnostics ===\n")
  cat("successful reps:", sum(ok), "/", cfg$n_reps, "\n")
  cat("failed reps:", length(failed), "\n")
  if (length(failed) > 0L) cat("failed rep ids:", paste(failed, collapse = ", "), "\n")

  plot_recovery(res_df, cfg)
}

if (sys.nframe() == 0L) {
  main()
}
