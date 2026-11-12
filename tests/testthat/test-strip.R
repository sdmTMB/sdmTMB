test_that("strip_sdmTMB() shrinks a fit without changing it or the original", {
  skip_on_cran()
  fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = pcod_mesh_2011,
    family = tweedie())
  p <- predict(fit)
  small <- strip_sdmTMB(fit)
  no_jp <- strip_sdmTMB(fit, strip_joint_precision = TRUE)
  size <- function(x) length(serialize(x, NULL))
  expect_lt(size(small), size(fit))
  expect_lt(size(no_jp), size(small))
  expect_null(no_jp$sd_report$jointPrecision)
  # the original keeps its live tape and Cholesky factor
  expect_false(is.null(fit$tmb_obj$env$L.created.by.newton))
  expect_false(is.null(attr(fit$tmb_obj$env$ADFun, "rcpp")))

  expect_equal(predict(small)$est, p$est)
  set.seed(1)
  s <- predict(small, nsim = 2)
  expect_equal(dim(s), c(nrow(pcod_2011), 2L))

  expect_message(x <- spread_sims(no_jp, nsim = 3), "getJointPrecision")
  expect_equal(nrow(x), 3L)

  f <- tempfile(fileext = ".rds")
  saveRDS(no_jp, f)
  reloaded <- readRDS(f)
  expect_equal(predict(reloaded)$est, p$est)
  set.seed(1)
  r <- residuals(reloaded, type = "mle-eb")
  set.seed(1)
  expect_equal(r, residuals(fit, type = "mle-eb"))
})
