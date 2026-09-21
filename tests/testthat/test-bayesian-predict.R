test_that("MCMC predictions work with an unfit sdmTMB object", {
  set.seed(1)
  dat <- data.frame(
    X = stats::runif(80),
    Y = stats::runif(80),
    a1 = stats::rnorm(80),
    observed = stats::rnorm(80)
  )
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 0.2)
  fit <- sdmTMB(
    observed ~ a1,
    data = dat,
    mesh = mesh,
    family = gaussian(),
    spatial = "on",
    do_fit = FALSE
  )

  samps <- matrix(fit$tmb_obj$env$last.par, ncol = 1L)
  expect_silent(pred <- predict(fit, mcmc_samples = samps))
  expect_equal(dim(pred), c(nrow(dat), 1L))
})
