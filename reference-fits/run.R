#!/usr/bin/env Rscript

# Reference-fit regression suite. Run from the package root:
#
#   Rscript reference-fits/run.R check    [--backend tmb|rtmb] [--filter REGEX] [--cores N] [--installed]
#   Rscript reference-fits/run.R validate [--backend tmb|rtmb|both] [--filter REGEX] [--cores N]
#   Rscript reference-fits/run.R record   [--filter REGEX] [--cores N]
#
# `check` refits every case and compares it with ref/reference.csv.
# `validate` fits with standard errors, runs sanity(), and reports timing and
# TMB vs RTMB agreement without writing anything.
# `record` is `validate` for both backends; if every case passes, it writes
# the TMB values to ref/reference.csv (with --filter, only those cases are
# replaced). See README.md.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else ""
if (!mode %in% c("check", "validate", "record")) {
  cat("Usage: Rscript reference-fits/run.R check|validate|record [options]\n")
  quit(status = 1L)
}
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
backend_opt <- opt("backend", if (mode == "check") "tmb" else "both")
backends <- if (mode == "record" || backend_opt == "both") c("tmb", "rtmb") else backend_opt
filter <- opt("filter", NULL)
cores <- as.integer(opt("cores", 1L))

rf_dir <- "reference-fits"
ref_file <- file.path(rf_dir, "ref", "reference.csv")
if (!dir.exists(rf_dir)) stop("Run from the package root.", call. = FALSE)

if ("--installed" %in% args) {
  suppressPackageStartupMessages(library(sdmTMB))
} else {
  pkgload::load_all(quiet = TRUE)
}
options(cli.num_colors = 1L, lifecycle_verbosity = "quiet")

source(file.path(rf_dir, "cases.R"), local = TRUE)
if (!is.null(filter)) cases <- cases[grepl(filter, names(cases))]
if (!length(cases)) stop("No cases match the filter.", call. = FALSE)

# Tolerances --------------------------------------------------------------------

# A metric passes if |new - ref| <= abs + rel * |ref|. Case-level overrides
# live in cases.R (`tol`).
default_tol <- list(
  nll = c(abs = 1e-6, rel = 1e-7),
  loglik = c(abs = 1e-6, rel = 1e-7),
  n_fixed = c(abs = 0, rel = 0),
  se = c(abs = 1e-6, rel = 1e-4),
  other = c(abs = 1e-6, rel = 1e-5)
)
metric_tol <- function(metric, case_tol) {
  key <- if (metric %in% names(default_tol)) metric else
    if (grepl("_se_", metric)) "se" else "other"
  tol <- default_tol[[key]]
  override <- case_tol[[metric]]
  if (is.null(override)) override <- case_tol[[key]]
  if (!is.null(override)) tol[names(override)] <- override
  tol
}
within_tol <- function(ref, new, tol) {
  !is.na(new) & abs(new - ref) <= tol[["abs"]] + tol[["rel"]] * abs(ref)
}

# Fitting and metrics ---------------------------------------------------------------

summarize_pred <- function(p, prefix, delta) {
  cols <- if (delta) c(est1 = "est1", est2 = "est2") else c(est = "est")
  out <- numeric()
  for (nm in names(cols)) {
    v <- p[[cols[[nm]]]]
    out[paste0(prefix, "_", nm, "_mean")] <- mean(v)
    out[paste0(prefix, "_", nm, "_sd")] <- stats::sd(v)
  }
  out
}

case_metrics <- function(fit, case) {
  nd <- case$newdata
  off <- if (!is.null(case$pred_offset)) nd[[case$pred_offset]]
  delta <- isTRUE(fit$family$delta)
  out <- c(
    nll = fit$model$objective,
    loglik = as.numeric(stats::logLik(fit)),
    n_fixed = length(fit$model$par)
  )
  out <- c(out, summarize_pred(predict(fit), "fitted", delta))
  out <- c(out, summarize_pred(predict(fit, newdata = nd, offset = off), "link", delta))
  pr <- predict(fit, newdata = nd, offset = off, type = "response")
  out["response_mean"] <- mean(pr$est)
  if ("re_form_na" %in% case$extras) {
    p <- predict(fit, newdata = nd, offset = off, re_form = NA)
    out <- c(out, summarize_pred(p, "re_form_na", delta))
  }
  if ("index" %in% case$extras) {
    ind <- get_index(fit, newdata = nd, offset = off, area = nd$area,
      bias_correct = FALSE)
    out["index_log_est_mean"] <- mean(ind$log_est)
    out["index_se_mean"] <- mean(ind$se)
  }
  if ("cog" %in% case$extras) {
    cog <- get_cog(fit, newdata = nd, offset = off, area = nd$area,
      bias_correct = FALSE, format = "wide")
    out["cog_x_mean"] <- mean(cog$est_x)
    out["cog_y_mean"] <- mean(cog$est_y)
    out["cog_se_x_mean"] <- mean(cog$se_x)
  }
  out
}

run_case <- function(name, backend, getsd) {
  case <- cases[[name]]
  a <- case$args
  if (is.null(a$control)) a$control <- sdmTMBcontrol()
  a$control$backend <- backend
  a$control$getsd <- getsd
  # The first (no-random-effect) phase can step outside the valid region for
  # identity and inverse links; a single phase is also faster here.
  a$control$multiphase <- FALSE
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    fit <- suppressMessages(suppressWarnings(do.call(sdmTMB, a)))
    metrics <- suppressMessages(suppressWarnings(case_metrics(fit, case)))
    sanity_ok <- NA
    sanity_fail <- ""
    if (getsd) {
      s <- suppressMessages(sanity(fit, silent = TRUE))
      checks <- unlist(s[setdiff(names(s), "all_ok")])
      failed <- names(checks)[!checks]
      sanity_ok <- all(failed %in% case$allow_sanity)
      sanity_fail <- paste(failed, collapse = ",")
    }
    list(metrics = metrics, sanity_ok = sanity_ok, sanity_fail = sanity_fail,
      error = NA_character_)
  }, error = function(e) {
    list(metrics = numeric(), sanity_ok = FALSE, sanity_fail = "",
      error = conditionMessage(e))
  })
  res$seconds <- proc.time()[["elapsed"]] - t0
  res
}

run_all <- function(backend, getsd) {
  cat(sprintf("Fitting %d cases with %s...\n", length(cases), backend))
  out <- parallel::mclapply(names(cases), run_case, backend = backend,
    getsd = getsd, mc.cores = cores, mc.preschedule = FALSE)
  names(out) <- names(cases)
  out
}

to_long <- function(results) {
  rows <- lapply(names(results), function(nm) {
    m <- results[[nm]]$metrics
    if (!length(m)) return(NULL)
    data.frame(case = nm, metric = names(m), value = unname(m))
  })
  do.call(rbind, rows)
}

# Compare two long tables; `ref` supplies the tolerance scale.
compare_long <- function(ref, new) {
  key <- function(x) paste(x$case, x$metric, sep = "::")
  m <- merge(ref, new, by = c("case", "metric"), all = TRUE,
    suffixes = c("_ref", "_new"))
  m$ok <- mapply(function(case, metric, r, n) {
    if (is.na(r) || is.na(n)) return(FALSE)
    within_tol(r, n, metric_tol(metric, cases[[case]]$tol))
  }, m$case, m$metric, m$value_ref, m$value_new)
  m$rel_diff <- abs(m$value_new - m$value_ref) / pmax(abs(m$value_ref), 1e-12)
  m[order(m$case, m$metric), ]
}

# Recorded for information only; never compared.
timing_metrics <- c("seconds_tmb", "seconds_rtmb")

format_value <- function(x) trimws(formatC(x, digits = 12, format = "g"))

print_failures <- function(cmp, label) {
  bad <- cmp[!cmp$ok, , drop = FALSE]
  cat(sprintf("\n%s: %d of %d metrics differ (%d cases).\n", label, nrow(bad),
    nrow(cmp), length(unique(bad$case))))
  if (nrow(bad)) {
    print(data.frame(case = bad$case, metric = bad$metric,
      ref = format_value(bad$value_ref), new = format_value(bad$value_new),
      rel_diff = signif(bad$rel_diff, 3)), row.names = FALSE, right = FALSE)
  }
  invisible(nrow(bad) == 0L)
}

# Local (untracked) timing history: compares each run with the previous run of
# the same mode/backend/case, flags slowdowns, then stores the new timings.
timing_file <- file.path(rf_dir, "timings-local.csv")
slow_rel <- 0.10 # flag if at least 10% slower than the previous run
slow_abs <- 0.02 # ... and at least this many seconds slower (ignore timer noise)

flag_slow_timings <- function(results, backend) {
  new <- data.frame(mode = mode, backend = backend, case = names(results),
    seconds = round(vapply(results, `[[`, numeric(1), "seconds"), 2))
  old <- if (file.exists(timing_file)) utils::read.csv(timing_file,
    stringsAsFactors = FALSE) else new[0, ]
  prev <- old[old$mode == mode & old$backend == backend, ]
  m <- merge(new, prev[c("case", "seconds")], by = "case",
    suffixes = c("", "_prev"))
  slow <- m[m$seconds - m$seconds_prev >= slow_abs &
    m$seconds >= (1 + slow_rel) * m$seconds_prev, , drop = FALSE]
  if (nrow(slow)) {
    slow$change <- sprintf("+%.0f%%", 100 * (slow$seconds / slow$seconds_prev - 1))
    cat(sprintf("\n%s (%s): %d case(s) >= %.0f%% slower than the previous run:\n",
      backend, mode, nrow(slow), 100 * slow_rel))
    print(slow[order(-slow$seconds / slow$seconds_prev),
      c("case", "seconds_prev", "seconds", "change")], row.names = FALSE, right = FALSE)
  } else if (nrow(m)) {
    cat(sprintf("\n%s (%s): no cases >= %.0f%% slower than the previous run.\n",
      backend, mode, 100 * slow_rel))
  }
  keep <- old[!(old$mode == mode & old$backend == backend & old$case %in% new$case), ]
  out <- rbind(keep, new)
  utils::write.csv(out[order(out$mode, out$backend, out$case), ], timing_file,
    row.names = FALSE)
}

report_status <- function(results, backend) {
  tab <- data.frame(
    case = names(results),
    seconds = round(vapply(results, `[[`, numeric(1), "seconds"), 2),
    sanity = vapply(results, function(r) if (is.na(r$sanity_ok)) "" else
      if (r$sanity_ok) "ok" else paste("FAIL", r$sanity_fail), character(1)),
    error = vapply(results, `[[`, character(1), "error")
  )
  utils::write.csv(tab, file.path(rf_dir, paste0("validate-", backend, ".csv")),
    row.names = FALSE)
  cat(sprintf("\n%s: %.2f s total\n", backend, sum(tab$seconds)))
  flag_slow_timings(results, backend)
  bad <- tab[!is.na(tab$error) | tab$sanity != "ok", , drop = FALSE]
  bad$error <- substr(gsub("\\s+", " ", bad$error), 1, 70)
  if (nrow(bad)) print(bad, row.names = FALSE, right = FALSE)
  slow <- utils::head(tab[order(-tab$seconds), c("case", "seconds")], 8)
  cat("Slowest cases:\n")
  print(slow, row.names = FALSE, right = FALSE)
  invisible(nrow(bad) == 0L)
}

# Modes ---------------------------------------------------------------------------

if (mode == "check") {
  if (!file.exists(ref_file)) stop("No reference file; run `record` first.", call. = FALSE)
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
  ref <- ref[ref$case %in% names(cases) & !ref$metric %in% timing_metrics, ]
  ok <- TRUE
  for (b in backends) {
    results <- run_all(b, getsd = FALSE)
    errors <- Filter(function(r) !is.na(r$error), results)
    for (nm in names(errors)) cat(sprintf("ERROR %s: %s\n", nm, errors[[nm]]$error))
    flag_slow_timings(results, b)
    cmp <- compare_long(ref, to_long(results))
    out_file <- file.path(rf_dir, paste0("check-", b, ".csv"))
    utils::write.csv(cmp, out_file, row.names = FALSE)
    ok <- print_failures(cmp, paste(b, "vs reference")) && ok
  }
  quit(status = if (ok) 0L else 1L)
}

# validate / record: fit with standard errors so sanity() can run.
results <- lapply(stats::setNames(backends, backends), run_all, getsd = TRUE)
ok <- TRUE
for (b in backends) ok <- report_status(results[[b]], b) && ok
if (length(backends) == 2L) {
  cmp <- compare_long(to_long(results$tmb), to_long(results$rtmb))
  ok <- print_failures(cmp, "RTMB vs TMB") && ok
  cat("\nLargest TMB vs RTMB relative differences by metric type:\n")
  cmp$type <- sub("_(est|est1|est2)_(mean|sd)$", "", cmp$metric)
  print(stats::aggregate(rel_diff ~ type, cmp, max), row.names = FALSE)
}

if (mode == "record") {
  if (!ok) {
    cat("\nNot recording: fix the failures above first.\n")
    quit(status = 1L)
  }
  new <- to_long(results$tmb)
  for (b in backends) {
    new <- rbind(new, data.frame(case = names(results[[b]]),
      metric = paste0("seconds_", b),
      value = round(vapply(results[[b]], `[[`, numeric(1), "seconds"), 2)))
  }
  if (!is.null(filter) && file.exists(ref_file)) {
    old <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
    new <- rbind(old[!old$case %in% names(cases), ], new)
  }
  new <- new[order(new$case, new$metric), ]
  new$value <- format_value(new$value)
  dir.create(dirname(ref_file), showWarnings = FALSE)
  utils::write.csv(new, ref_file, row.names = FALSE, quote = FALSE)
  pkgs <- c("TMB", "RTMB", "Matrix", "fmesher", "mgcv")
  writeLines(c(
    paste("recorded:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("git:", system("git rev-parse --short HEAD", intern = TRUE)),
    paste("R:", R.version.string),
    paste("platform:", R.version$platform),
    paste0(pkgs, ": ", vapply(pkgs, function(p) as.character(utils::packageVersion(p)), ""))
  ), file.path(rf_dir, "ref", "session.txt"))
  cat(sprintf("\nWrote %d metrics for %d cases to %s\n", nrow(new),
    length(unique(new$case)), ref_file))
}
quit(status = if (ok) 0L else 1L)
