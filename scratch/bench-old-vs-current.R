# Benchmark sdmTMB: old commit (TMB only) vs current TMB vs current RTMB.
#
# Steps (pcod, mesh cutoff = 8, spatial + spatiotemporal IID fields, year
# factors, qcs_grid newdata): fit, predict, get_index(bias_correct = FALSE/TRUE),
# simulate(nsim = 400), plus project(nsim = 400) for a 30-year AR1 projection
# (dogfish, as in ?project). If project() is missing from a library, it is NA.
#
# Each configuration runs in its own fresh R process against its own installed
# library (non-debug builds). Peak memory is the peak resident set size of that
# process (sampled with `ps` every ~50 ms) during each step, not a delta.
#
# Usage (from package root):
#   Rscript scratch/bench-old-vs-current.R              # full run
#   Rscript scratch/bench-old-vs-current.R --reps=5     # more timing reps
#   Rscript scratch/bench-old-vs-current.R --reinstall  # also redo the old install
#   Rscript scratch/bench-old-vs-current.R --md-only    # rewrite .md/.csv from the last run
#
# Old library is cached; the current version is reinstalled from the working
# tree (copied, so your src/*.o files are untouched) on every run.

old_sha <- "891d8e7f914759338ef91d078032ee85d2ed778b"
cache_dir <- tools::R_user_dir("sdmTMB-bench", "cache")

# ---- worker ----------------------------------------------------------------
run_worker <- function(label, lib, backend, out_file, reps) {
  .libPaths(c(lib, .libPaths()))
  suppressMessages(library(sdmTMB))
  cat(label, ": sdmTMB", as.character(packageVersion("sdmTMB")), "from",
    find.package("sdmTMB"), "\n")
  has_backend <- "backend" %in% names(formals(sdmTMBcontrol))
  ctl <- if (backend == "default") sdmTMBcontrol() else sdmTMBcontrol(backend = backend)

  # background RSS sampler
  pid <- Sys.getpid()
  rss_file <- tempfile()
  file.create(rss_file)
  system(sprintf(
    "while kill -0 %d 2>/dev/null; do ps -o rss= -p %d >> %s; sleep 0.05; done",
    pid, pid, rss_file
  ), wait = FALSE)
  Sys.sleep(0.3)
  rss_now <- function() as.numeric(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE))
  n_lines <- function() length(readLines(rss_file))

  res <- list()
  step <- function(name, expr, n = reps) {
    expr <- substitute(expr)
    env <- parent.frame()
    times <- numeric(n)
    peak <- rss_now()
    for (i in seq_len(n)) {
      gc(full = TRUE)
      n0 <- n_lines()
      times[i] <- system.time(val <- eval(expr, env))[["elapsed"]]
      r <- c(as.numeric(readLines(rss_file)[-seq_len(n0)]), rss_now())
      peak <- max(peak, r)
    }
    cat(sprintf("%-28s %7.2f s (median of %d)  peak %6.0f MB\n",
      name, median(times), n, peak / 1024))
    res[[name]] <<- data.frame(
      config = label, step = name, seconds = median(times), peak_mb = peak / 1024
    )
    invisible(val)
  }
  fail <- function(name) {
    res[[name]] <<- data.frame(config = label, step = name, seconds = NA, peak_mb = NA)
  }

  # pcod: fit / predict / index / simulate
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 8)
  fit <- step("fit", sdmTMB(
    density ~ 0 + as.factor(year), data = pcod, mesh = mesh, time = "year",
    family = delta_gamma(), spatial = "on", spatiotemporal = "iid",
    control = ctl, silent = TRUE
  ), n = 1)
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  step("predict", predict(fit, newdata = nd))
  pred <- predict(fit, newdata = nd, return_tmb_object = TRUE)
  step("get_index(bias_correct = FALSE)", get_index(pred, bias_correct = FALSE))
  step("get_index(bias_correct = TRUE)", get_index(pred, bias_correct = TRUE))
  step("simulate(nsim = 400)", simulate(fit, nsim = 400))
  rm(fit, pred)

  # dogfish: project() with 30 years and AR1 (as in ?project)
  if ("project" %in% getNamespaceExports("sdmTMB")) {
    mesh_d <- make_mesh(dogfish, c("X", "Y"), cutoff = 10)
    hist_years <- 2004:2022
    proj_grid <- replicate_df(wcvi_grid, "year", c(hist_years, max(hist_years) + 1:30))
    fit_d <- sdmTMB(
      catch_weight ~ 1, time = "year", offset = log(dogfish$area_swept),
      extra_time = hist_years, spatial = "on", spatiotemporal = "ar1",
      data = dogfish, mesh = mesh_d, family = tweedie(link = "log"),
      control = ctl, silent = TRUE
    )
    step("project(nsim = 400, 30 yr AR1)", {
      set.seed(1)
      project(fit_d, newdata = proj_grid, nsim = 400, silent = TRUE)
    })
  } else {
    fail("project(nsim = 400, 30 yr AR1)")
  }
  saveRDS(do.call(rbind, res), out_file)
}

# ---- outputs ---------------------------------------------------------------
labs <- c("old TMB", "current TMB", "current RTMB")

write_outputs <- function(res, out_dir) {
  steps <- unique(res$step)
  get <- function(v) sapply(labs, function(l) res[[v]][match(paste(l, steps), paste(res$config, res$step))])
  tt <- get("seconds"); mm <- get("peak_mb")
  tab <- data.frame(
    step = steps,
    old_s = tt[, "old TMB"], tmb_s = tt[, "current TMB"], rtmb_s = tt[, "current RTMB"],
    rtmb_speedup_vs_old = tt[, "old TMB"] / tt[, "current RTMB"],
    old_mb = mm[, "old TMB"], tmb_mb = mm[, "current TMB"], rtmb_mb = mm[, "current RTMB"]
  )
  print(format(tab, digits = 3), row.names = FALSE)
  write.csv(tab, file.path(out_dir, "bench-old-vs-current-results.csv"), row.names = FALSE)

  fmt <- function(x, digits) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
  f1 <- function(x) fmt(x, 1)
  f2 <- function(x) fmt(x, 2)
  md_table <- function(header, rows) c(
    paste0("| ", paste(header, collapse = " | "), " |"),
    paste0("|", paste(c(":--", rep("--:", length(header) - 1)), collapse = "|"), "|"),
    vapply(rows, function(r) paste0("| ", paste(r, collapse = " | "), " |"), "")
  )
  time_rows <- lapply(seq_len(nrow(tab)), function(i) with(tab[i, ], c(
    paste0("`", step, "`"), f2(old_s), f2(tmb_s), f2(rtmb_s),
    ifelse(is.na(rtmb_speedup_vs_old), "NA", paste0(f1(rtmb_speedup_vs_old), "x"))
  )))
  mem_rows <- lapply(seq_len(nrow(tab)), function(i) with(tab[i, ], c(
    paste0("`", step, "`"), f1(old_mb), f1(tmb_mb), f1(rtmb_mb)
  )))
  md <- c(
    "## sdmTMB benchmark: old commit vs current TMB vs current RTMB",
    "",
    paste0("Compares sdmTMB at commit `", old_sha, "` (TMB only) with the current ",
      "version using the TMB backend and the default RTMB backend. Each ",
      "configuration ran in a fresh R process from a non-debug install."),
    "",
    "**Models**",
    "",
    "- `fit`, `predict`, `get_index()`, `simulate()`: `pcod` data, mesh with `cutoff = 8`, ",
    "  `density ~ 0 + as.factor(year)`, delta-gamma family (`delta_gamma()`), `time = \"year\"`, ",
    "  spatial and IID spatiotemporal random fields. `predict()` and `get_index()` use ",
    "  `qcs_grid` replicated across all years; `get_index()` is run with `bias_correct` ",
    "  `FALSE` and `TRUE`; `simulate()` uses `nsim = 400`.",
    "- `project()`: `dogfish` data (as in `?project`), mesh with `cutoff = 10`, Tweedie ",
    "  with an area-swept offset, spatial field and AR1 spatiotemporal field, fitted to ",
    "  2004-2022 and projected 30 further years on `wcvi_grid` with `nsim = 400`.",
    "",
    "Times are elapsed seconds (median over reps; `fit` is run once). Speedup is ",
    "old TMB time divided by current RTMB time.",
    "",
    "### Time (s)",
    "",
    md_table(c("step", "old TMB", "current TMB", "current RTMB", "RTMB speedup vs old"), time_rows),
    "",
    "### Peak memory (MB)",
    "",
    md_table(c("step", "old TMB", "current TMB", "current RTMB"), mem_rows),
    "",
    "Peak memory is the peak process RSS during the step (includes baseline session memory)."
  )
  writeLines(md, file.path(out_dir, "bench-old-vs-current-results.md"))
  cat("\nWrote", file.path(out_dir, "bench-old-vs-current-results.md"), "\n")
}

# ---- orchestrator ----------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) && args[1] == "--worker") {
  run_worker(args[2], args[3], args[4], args[5], as.integer(args[6]))
  quit(save = "no")
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script <- normalizePath(script)
pkg_root <- normalizePath(file.path(dirname(script), ".."))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
results_file <- file.path(cache_dir, "last-results.rds")

# rewrite the outputs from the last run without rerunning anything
if ("--md-only" %in% args) {
  write_outputs(readRDS(results_file), dirname(script))
  quit(save = "no")
}

reps <- 3L
if (any(grepl("^--reps=", args))) reps <- as.integer(sub("^--reps=", "", grep("^--reps=", args, value = TRUE)))
reinstall <- "--reinstall" %in% args


install_lib <- function(src, lib) {
  unlink(lib, recursive = TRUE)
  dir.create(lib, recursive = TRUE)
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
      paste0("--library=", shQuote(lib)), shQuote(src)),
    env = "MAKEFLAGS=-j4")
  if (status != 0) stop("install failed for ", src)
}

old_lib <- file.path(cache_dir, paste0("lib-", substr(old_sha, 1, 8)))
if (reinstall || !dir.exists(file.path(old_lib, "sdmTMB"))) {
  message("Installing old sdmTMB (", substr(old_sha, 1, 8), ")...")
  src <- file.path(tempdir(), "sdmTMB-old")
  dir.create(src)
  system(sprintf("git -C %s archive %s | tar -x -C %s",
    shQuote(pkg_root), old_sha, shQuote(src)))
  install_lib(src, old_lib)
}

message("Installing current sdmTMB from working tree...")
src <- file.path(tempdir(), "sdmTMB-current")
dir.create(src)
system(sprintf(
  "rsync -a --exclude .git --exclude scratch --exclude vignettes --exclude '*.o' --exclude '*.so' %s/ %s/",
  shQuote(pkg_root), shQuote(src)))
cur_lib <- file.path(cache_dir, "lib-current")
install_lib(src, cur_lib)

configs <- list(
  list(label = "old TMB", lib = old_lib, backend = "default"),
  list(label = "current TMB", lib = cur_lib, backend = "tmb"),
  list(label = "current RTMB", lib = cur_lib, backend = "rtmb")
)
out <- lapply(configs, function(cf) {
  out_file <- tempfile(fileext = ".rds")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), "--worker", shQuote(cf$label), shQuote(cf$lib),
      cf$backend, shQuote(out_file), reps))
  if (status != 0 || !file.exists(out_file)) stop("worker failed: ", cf$label)
  readRDS(out_file)
})
res <- do.call(rbind, out)
rownames(res) <- NULL
saveRDS(res, results_file)

write_outputs(res, dirname(script))
