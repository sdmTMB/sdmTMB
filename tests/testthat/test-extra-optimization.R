test_that("Extra optimization runs and reduces gradients", {
  skip_on_cran()
  d <- subset(pcod, year >= 2013)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"))
  lpb <- m$tmb_obj$env$last.par.best

  m1 <- run_extra_optimization(m, nlminb_loops = 2, newton_loops = 1)
  expect_lt(max(m1$gradients), max(m$gradients))
  expect_lte(m1$model$objective, m$model$objective + 1e-8) # Newton accepts roundoff-level ties
  expect_gt(m1$model$iterations, m$model$iterations)

  # the original fit is untouched
  expect_identical(m$tmb_obj$env$last.par.best, lpb)

  # saved parameters match the new optimum, so predict() uses them
  expect_identical(m1$last.par.best, m1$tmb_obj$env$last.par.best)
  expect_equal(m1$parlist$b_j, unname(m1$model$par[names(m1$model$par) == "b_j"]))
  expect_equal(
    predict(m1)$est,
    predict(m1, newdata = d)$est,
    tolerance = 1e-6
  )
  expect_identical(m1$pos_def_hessian, m1$sd_report$pdHess)
})

test_that("Newton updates shrink steps to stay within parameter limits", {
  obj <- list(
    fn = function(x) sum((x + 1)^2),
    gr = function(x) 2 * (x + 1)
  )
  opt <- list(par = c(x = 0.5), objective = obj$fn(0.5))

  # the full step (to -1) crosses the lower limit; a halved step doesn't
  out <- run_newton_loops(
    newton_loops = 1L,
    opt = opt,
    obj = obj,
    lower = c(x = 0),
    upper = c(x = Inf)
  )
  expect_gte(out$par[["x"]], 0)
  expect_lt(out$objective, opt$objective)

  # no step length stays within the limits
  out <- run_newton_loops(
    newton_loops = 1L,
    opt = opt,
    obj = obj,
    lower = c(x = 0.5),
    upper = c(x = Inf)
  )
  expect_equal(out$par, opt$par)
  expect_equal(out$objective, opt$objective)
})

test_that("Newton updates survive Hessian failures and empty parameters", {
  obj <- list(
    fn = function(x) sum((x + 1)^2),
    gr = function(x) if (length(x) && x[[1]] == 0.5) 1 else NaN
  )
  opt <- list(par = c(x = 0.5), objective = obj$fn(0.5))
  expect_equal(run_newton_loops(1L, opt, obj), opt)

  empty <- list(par = numeric(0), objective = 1)
  expect_equal(run_newton_loops(1L, empty, obj), empty)
})

test_that("nlminb loops accumulate iterations and keep the best fit", {
  obj <- list(
    fn = function(x) sum((x - 3)^2),
    gr = function(x) 2 * (x - 3)
  )
  opt <- stats::nlminb(0, obj$fn, obj$gr)
  out <- run_nlminb_loops(2L, opt, obj, lower = -Inf, upper = Inf,
    control = list())
  expect_lte(out$objective, opt$objective)
  expect_gt(out$iterations, opt$iterations)

  # a failing restart retains the previous fit
  bad <- list(fn = function(x) NaN, gr = function(x) NaN)
  expect_equal(
    run_nlminb_loops(1L, opt, bad, lower = -Inf, upper = Inf, control = list(),
      suppress_warnings = TRUE),
    opt
  )
})

test_that("sdmTMBcontrol() stores suppress_nlminb_warnings", {
  expect_true(sdmTMBcontrol()$suppress_nlminb_warnings)
  expect_false(sdmTMBcontrol(suppress_nlminb_warnings = FALSE)$suppress_nlminb_warnings)
})
