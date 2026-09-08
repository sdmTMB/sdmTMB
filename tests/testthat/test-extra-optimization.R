test_that("Extra optimization runs and reduces gradients", {
  skip_on_cran()
  d <- subset(pcod, year >= 2013)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
    data = d, time = "year", mesh = pcod_spde, family = tweedie(link = "log"))

  m1 <- run_extra_optimization(m, nlminb_loops = 1, newton_loops = 1)
  expect_lt(max(m1$gradients), max(m$gradients))
})

test_that("Newton updates respect parameter limits", {
  obj <- list(
    fn = function(x) sum((x + 1)^2),
    gr = function(x) 2 * (x + 1)
  )
  opt <- list(par = c(x = 0.5), objective = obj$fn(0.5))

  out <- run_newton_loops(
    newton_loops = 1L,
    opt = opt,
    obj = obj,
    lower = c(x = 0),
    upper = c(x = Inf)
  )

  expect_equal(out$par, opt$par)
  expect_equal(out$objective, opt$objective)
})
