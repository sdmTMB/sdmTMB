# Runs every method on fitted models as the first call on a fresh stripped
# copy (null tape pointer, as after readRDS()) and on an unstripped copy, so a
# method that skips reinitialize() shows up as stripped-only failure.
# A segfault ends the script; the last "running" line names the method.
#
# Rscript scratch/strip-methods-check.R

devtools::load_all(quiet = TRUE)

years <- sort(unique(pcod_2011$year))
fit <- sdmTMB(density ~ depth_scaled + depth_scaled2, data = pcod_2011,
  mesh = pcod_mesh_2011, time = "year", spatiotemporal = "iid",
  anisotropy = TRUE, family = tweedie())
fit_s <- sdmTMB(density ~ s(depth, k = 5), data = pcod_2011,
  mesh = pcod_mesh_2011, family = tweedie())
fit_re <- sdmTMB(density ~ depth_scaled + (1 | fyear),
  data = transform(pcod_2011, fyear = factor(year)), mesh = pcod_mesh_2011,
  family = tweedie())
dfit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = pcod_mesh_2011,
  family = delta_gamma())
nd <- replicate_df(qcs_grid[seq(1, nrow(qcs_grid), by = 50), ], "year", years)
nd_proj <- replicate_df(pcod_2011[seq(1, nrow(pcod_2011), by = 100), ],
  "year", c(years, max(years) + 2))
cv <- sdmTMB_cv(density ~ depth_scaled, data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie(), k_folds = 2, save_models = TRUE)

methods <- list(
  print = function(m) utils::capture.output(print(m)),
  summary = function(m) utils::capture.output(summary(m)),
  tidy = function(m) tidy(m, "ran_pars", conf.int = TRUE),
  coef = coef, confint = confint, vcov = vcov, fixef = fixef,
  sigma = sigma, family = family, formula = formula,
  terms = terms, nobs = nobs, logLik = logLik, AIC = stats::AIC,
  extractAIC = extractAIC, deviance = deviance, df.residual = df.residual,
  fitted = fitted, model.frame = model.frame, model.matrix = model.matrix,
  cAIC = cAIC, cAIC_edf = function(m) cAIC(m, what = "EDF"),
  compare_deviance = function(m) compare_deviance(m, fit),
  sanity = function(m) sanity(m, silent = TRUE),
  get_pars = get_pars,
  get_prior_parameters = get_prior_parameters,
  predict = predict,
  predict_se = function(m) predict(m, newdata = nd[1:5, ], se_fit = TRUE),
  predict_nsim = function(m) predict(m, nsim = 2),
  residuals = residuals,
  residuals_mvn = function(m) residuals(m, type = "mle-mvn"),
  simulate = function(m) simulate(m, nsim = 2),
  simulate_mvn = function(m) simulate(m, nsim = 2, type = "mle-mvn"),
  dharma_residuals = function(m) {
    dharma_residuals(simulate(m, nsim = 20), m, plot = FALSE)
  },
  spread_sims = function(m) spread_sims(m, nsim = 2),
  gather_sims = function(m) gather_sims(m, nsim = 2),
  get_index = function(m) get_index(m, newdata = nd, bias_correct = FALSE),
  get_index_bc = function(m) get_index(m, newdata = nd, bias_correct = TRUE),
  get_index_split = function(m) {
    get_index_split(m, newdata = nd, nsplit = 2, bias_correct = FALSE)
  },
  get_cog = function(m) get_cog(m, newdata = nd),
  get_eao = function(m) get_eao(m, newdata = nd),
  get_weighted_average = function(m) {
    get_weighted_average(m, newdata = nd, vector = nd$depth)
  },
  get_index_sims = function(m) get_index_sims(predict(m, newdata = nd, nsim = 2)),
  get_range_edge = function(m) {
    get_range_edge(predict(m, newdata = nd, nsim = 2), axis = nd$Y)
  },
  project = function(m) project(m, nd_proj, nsim = 1, silent = TRUE),
  run_extra_optimization = function(m) {
    run_extra_optimization(m, nlminb_loops = 0, newton_loops = 1)
  },
  update = function(m) update(m),
  plot_anisotropy = plot_anisotropy, plot_anisotropy2 = plot_anisotropy2,
  emmeans = function(m) emmeans::emmeans(m, ~depth_scaled),
  Effect = function(m) effects::Effect("depth_scaled", m),
  ggpredict = function(m) ggeffects::ggpredict(m, "depth_scaled"),
  visreg = function(m) visreg::visreg(m, "depth_scaled", plot = FALSE)
)
smooth_methods <- list(
  print = function(m) utils::capture.output(print(m)),
  tidy = tidy,
  cAIC_edf = function(m) cAIC(m, what = "EDF"),
  plot_smooth = function(m) plot_smooth(m, return_data = TRUE),
  spread_sims = function(m) spread_sims(m, nsim = 2),
  predict_se = function(m) predict(m, newdata = nd[1:5, ], se_fit = TRUE),
  visreg = function(m) visreg::visreg(m, "depth", plot = FALSE),
  ggpredict = function(m) ggeffects::ggpredict(m, "depth")
)
re_methods <- list(
  ranef = ranef,
  tidy_ran_vals = function(m) tidy(m, "ran_vals"),
  cAIC = cAIC,
  predict_se = function(m) predict(m, newdata = fit_re$data[1:5, ], se_fit = TRUE),
  simulate_mvn = function(m) simulate(m, nsim = 2, type = "mle-mvn")
)
delta_methods <- list(
  print = function(m) utils::capture.output(print(m)),
  tidy = function(m) tidy(m, model = 2),
  predict = predict,
  predict_se = function(m) predict(m, newdata = nd[1:5, ], se_fit = TRUE),
  residuals = residuals,
  simulate = function(m) simulate(m, nsim = 2),
  simulate_mvn = function(m) simulate(m, nsim = 2, type = "mle-mvn"),
  cAIC = cAIC,
  get_index = function(m) get_index(m, newdata = nd, bias_correct = FALSE),
  visreg_delta = function(m) visreg_delta(m, "depth_scaled", model = 2, plot = FALSE)
)
cv_methods <- list(
  tidy = function(x) tidy(x),
  print = function(x) utils::capture.output(print(x)),
  elpd = function(x) loo::elpd(x),
  fold_print = function(x) utils::capture.output(print(x$models[[1]])),
  fold_predict = function(x) predict(x$models[[1]]),
  fold_cAIC = function(x) cAIC(x$models[[1]]),
  fold_residuals = function(x) residuals(x$models[[1]], type = "mle-mvn"),
  fold_get_index = function(x) {
    get_index(x$models[[1]], newdata = nd, bias_correct = FALSE)
  }
)

strip_fit <- function(m) strip_sdmTMB(m, strip_joint_precision = TRUE)
strip_cv <- function(x) {
  x$models <- lapply(x$models, strip_fit)
  x
}
# a serialize round trip copy keeps the original from being retaped
copy <- function(x) unserialize(serialize(x, NULL))

run <- function(label, methods, x, strip) {
  out <- lapply(names(methods), function(nm) {
    message("running ", label, ": ", nm)
    one <- function(y) {
      set.seed(1) # the same draws for both copies
      tryCatch({
        suppressWarnings(suppressMessages(methods[[nm]](y)))
        "ok"
      }, error = function(e) conditionMessage(e))
    }
    data.frame(model = label, method = nm,
      unstripped = one(copy(x)), stripped = one(strip(x)))
  })
  do.call(rbind, out)
}

res <- rbind(
  run("main", methods, fit, strip_fit),
  run("smooth", smooth_methods, fit_s, strip_fit),
  run("ranef", re_methods, fit_re, strip_fit),
  run("delta", delta_methods, dfit, strip_fit),
  run("cv", cv_methods, cv, strip_cv)
)
res$unstripped <- substr(gsub("\n", " ", res$unstripped), 1, 70)
res$stripped <- substr(gsub("\n", " ", res$stripped), 1, 70)
bad <- res[res$unstripped != "ok" | res$stripped != "ok", ]
cat("\n", sum(res$stripped == "ok"), "/", nrow(res), "ok when stripped\n")
if (nrow(bad)) print(bad, row.names = FALSE) else cat("All methods ok\n")
