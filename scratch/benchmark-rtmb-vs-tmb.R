# Benchmark: tmb vs rtmb backends
# Basic spatial model, increasing mesh size / n observations.
# Times and memory (via bench::mark / gc) for: sdmTMB() fit, predict(newdata),
# get_index(bias_correct = TRUE).
#
# Run from package root: Rscript scratch/benchmark-rtmb-vs-tmb.R

library(sdmTMB)
library(dplyr)

set.seed(1)

# ---- simulate data over a range of sizes -----------------------------------
# roughly scale n_obs with mesh size so this represents realistic surveys
sizes <- list(
  list(n_obs = 800,  cutoff = 40),
  list(n_obs = 2000, cutoff = 20),
  list(n_obs = 4000, cutoff = 8),
  list(n_obs = 6000, cutoff = 5),
  list(n_obs = 6000, cutoff = 3),
  list(n_obs = 6000, cutoff = 1.5),
  list(n_obs = 6000, cutoff = 1.1)
)

sim_dat <- function(n_obs, seed = 1) {
  set.seed(seed)
  predictor_dat <- data.frame(
    X = runif(n_obs, 0, 100),
    Y = runif(n_obs, 0, 100)
  )
  predictor_dat
}

bench_one <- function(n_obs, cutoff, backend, seed = 1) {
  dat_base <- sim_dat(n_obs, seed = seed)
  mesh <- make_mesh(dat_base, xy_cols = c("X", "Y"), cutoff = cutoff)
  n_verts <- mesh$mesh$n

  set.seed(seed)
  sim_out <- sdmTMB_simulate(
    formula = ~1,
    data = dat_base,
    mesh = mesh,
    family = gaussian(),
    range = 30,
    phi = 0.3,
    sigma_O = 0.6,
    seed = seed,
    B = 0.5
  )

  gc(full = TRUE)
  t_fit0 <- proc.time()
  mem_fit <- tryCatch({
    gc(reset = TRUE, full = TRUE)
    fit <- sdmTMB(
      observed ~ 1,
      data = sim_out,
      mesh = mesh,
      family = gaussian(),
      control = sdmTMBcontrol(backend = backend)
    )
    m2 <- gc(full = TRUE)
    # "max used" column (Mb) since the reset, summed over Ncells/Vcells rows
    list(fit = fit, mem = sum(m2[, 6]))
  }, error = function(e) list(fit = NULL, error = conditionMessage(e)))
  t_fit <- (proc.time() - t_fit0)[["elapsed"]]

  if (is.null(mem_fit$fit)) {
    return(data.frame(
      backend = backend, n_obs = n_obs, cutoff = cutoff, n_verts = n_verts,
      fit_time = NA, fit_maxmem_mb = NA, predict_time = NA,
      index_time = NA, converged = FALSE, error = mem_fit$error
    ))
  }
  fit <- mem_fit$fit

  # newdata for prediction: grid over the domain, one time value (no time here)
  newdata <- expand.grid(
    X = seq(0, 100, length.out = 50),
    Y = seq(0, 100, length.out = 50)
  )

  gc(full = TRUE)
  t_pred0 <- proc.time()
  pred <- tryCatch(
    predict(fit, newdata = newdata),
    error = function(e) {message("predict error: ", conditionMessage(e)); NULL}
  )
  t_pred <- (proc.time() - t_pred0)[["elapsed"]]

  gc(full = TRUE)
  t_idx0 <- proc.time()
  idx <- tryCatch(
    get_index(fit, newdata = newdata, bias_correct = TRUE),
    error = function(e) {message("get_index error: ", conditionMessage(e)); NULL}
  )
  t_idx <- (proc.time() - t_idx0)[["elapsed"]]

  data.frame(
    backend = backend,
    n_obs = n_obs,
    cutoff = cutoff,
    n_verts = n_verts,
    fit_time = t_fit,
    fit_maxmem_mb = mem_fit$mem,
    predict_time = t_pred,
    index_time = t_idx,
    converged = fit$pos_def_hessian %||% NA,
    max_gradient = tryCatch(max(abs(fit$gradients)), error = function(e) NA),
    error = NA_character_
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

results <- list()
i <- 1
for (s in sizes) {
  for (backend in c("tmb", "rtmb")) {
    cat(sprintf("\n=== backend=%s n_obs=%d cutoff=%g ===\n", backend, s$n_obs, s$cutoff))
    res <- bench_one(s$n_obs, s$cutoff, backend, seed = 1)
    print(res)
    results[[i]] <- res
    i <- i + 1
  }
}

out <- do.call(rbind, results)
print(out)
saveRDS(out, "scratch/benchmark-rtmb-vs-tmb-results.rds")
write.csv(out, "scratch/benchmark-rtmb-vs-tmb-results.csv", row.names = FALSE)

cat("\n\n==== Summary (fit_time in seconds) ====\n")
print(
  out |> select(backend, n_obs, n_verts, fit_time, fit_maxmem_mb, predict_time, index_time)
)

cat("\n\n==== tmb vs rtmb ratio (rtmb / tmb) ====\n")
wide_fit <- reshape(out[, c("backend", "n_obs", "fit_time")], idvar = "n_obs", timevar = "backend", direction = "wide")
wide_fit$ratio_fit <- wide_fit$fit_time.rtmb / wide_fit$fit_time.tmb
print(wide_fit)

wide_pred <- reshape(out[, c("backend", "n_obs", "predict_time")], idvar = "n_obs", timevar = "backend", direction = "wide")
wide_pred$ratio_predict <- wide_pred$predict_time.rtmb / wide_pred$predict_time.tmb
print(wide_pred)

wide_idx <- reshape(out[, c("backend", "n_obs", "index_time")], idvar = "n_obs", timevar = "backend", direction = "wide")
wide_idx$ratio_index <- wide_idx$index_time.rtmb / wide_idx$index_time.tmb
print(wide_idx)
