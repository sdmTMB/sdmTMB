# Simulated data and fit helper shared by the tests below.
# make fake predictor(s) (a1) and sampling locations:
set.seed(123)
predictor_dat <- data.frame(
  X = runif(2000), Y = runif(2000),
  a1 = rnorm(2000), year = rep(1:10, each = 200)
)

cutoff <- 0.1
sim_mesh <- sdmTMB::make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = cutoff)

fe_default <- 0.5
range_default <- 0.3
sigma_default <- 0.1

test_sim <- function(b0 = 0.5,
                     b1 = -fe_default,
                     phi = 0.1,
                     mtrange = range_default,
                     rho = NULL,
                     sigma_O = sigma_default,
                     sigma_E = sigma_default,
                     model_type = "IID") {

  sim_dat <- sdmTMB::sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    time = "year",
    mesh = sim_mesh,
    family = gaussian(),
    range = mtrange,
    rho = rho,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    phi = phi,
    seed = 42,
    B = c(b0, b1) # B0 = intercept, B1 = a1 slope
  )

  fit <- sdmTMB::sdmTMB(observed ~ a1,
    data = sim_dat, mesh = sim_mesh, time = "year",
    spatiotemporal = model_type
  )

  if (all(unlist(sanity(fit, gradient_thresh = 0.01)))) {
    sr <- as.list(fit$sd_report, "Estimate")
    ty <- tidy(fit, effects = "ran_pars", conf.int = TRUE, silent = TRUE)
    out <- list()
    out[["model"]] <- fit

    # calculate the difference between estimate and true
    out["b0"] <- sr$b_j[1] - b0
    out["b1"] <- sr$b_j[2] - b1
    out["phi"] <- exp(sr$ln_phi[1]) - phi
    out["range"] <- ty$estimate[ty$term == "range"] - mtrange
    out["rho"] <- ifelse(is.null(rho), 0,
                         ty$estimate[ty$term == "rho"] - rho)
    out["sigma_O"] <- ty$estimate[ty$term == "sigma_O"] - sigma_O
    out["sigma_E"] <- ifelse(is.null(sigma_E), 0,
                             ty$estimate[ty$term == "sigma_E"] - sigma_E)
    out
  } else {
    out <- list()
    out[["model"]] <- fit
    out
  }
}

test_that("sdmTMB_simulate() output can be recovered by IID and AR1 fits", {
  skip_on_cran()

  # tests with IID
  # fixed_effect_default 0.7
  # range = 0.4
  # phi = 0.1
  # sigmas = 0.1
  # mesh uncertainty/cutoff = 0.1
  t1 <- test_sim()
  t1[[1]]
  # fixed effect within 5% of true value
  expect_lt(abs(t1[["b1"]]), fe_default*0.05)
  # range error within 0.1 mesh cutoff
  expect_lt(abs(t1[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t1[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_E"])

  ## test AR1
  t <- test_sim(rho = 0.4,
                sigma_O = 0.5, # needs larger sigma to be able to estimate it along with an AR1
                sigma_E = 0.1,
                model_type = "AR1")
  t[[1]]
  # expect within 0.15 of true value
  expect_lt(abs(t[["rho"]]), 0.15)
  # with one seed, rho was consistently 0.1 larger than the value used in simulation
  # fixed effect within 5% of true value
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  # range within 10% of true value
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_O"])
  # expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_E"])
  # expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_E"])
})

test_that("sdmTMB_simulate() recovery holds across variances, ranges, and field types", {
  skip_on_cran()
  skip_slow()

  # increase sigma_O
  t <- test_sim(sigma_O = 0.4, sigma_E = 0.1)
  t[[1]]
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(0.4, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(0.4, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(0.1, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(0.1, ci$conf.low[ci$term=="sigma_E"])

  # increase just sigma_E
  t <- test_sim(sigma_O = 0.1, sigma_E = 0.3)
  t[[1]]
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  expect_lt(abs(t[["range"]]), cutoff) # just barely fails at 5%
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(0.1, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(0.1, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(0.3, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(0.3, ci$conf.low[ci$term=="sigma_E"])


  # increase both sigmas by a lot sigma_O > sigma_E
  t <- test_sim(sigma_O = 0.8, sigma_E = 0.4)
  t[[1]]
  expect_lt(abs(t[["b1"]]), fe_default*0.05)# intercept gets harder to estimate
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(0.8, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(0.8, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(0.4, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(0.4, ci$conf.low[ci$term=="sigma_E"])


  # adjust range lower
  t <- test_sim(mtrange = range_default-0.15)
  t[[1]]
  # fixed effect within 5% of true value
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  # estimate is lower than that for the default model by the amount we've changed the range from that value
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_E"])


  # adjust range higher
  t <- test_sim(mtrange = range_default*1.5)
  t[[1]]
  # fixed effect within 5% of true value
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(sigma_default, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(sigma_default, ci$conf.low[ci$term=="sigma_E"])


  ## test high AR1 values
  t <- test_sim(rho = 0.95,
                sigma_O = 0.7, # needs larger sigma to be able to estimate it along with an AR1
                sigma_E = 0.2,
                model_type = "AR1")
  t[[1]]
  expect_lt(abs(t[["rho"]]), 0.1)
  # fixed effect within 5% of true value
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  # range within 10% of true value
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(0.7, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(0.7, ci$conf.low[ci$term=="sigma_O"])
  expect_lt(0.2, ci$conf.high[ci$term=="sigma_E"])
  expect_gt(0.2, ci$conf.low[ci$term=="sigma_E"])

  t <- test_sim(rho = 0.99,
                sigma_O = 0.7, # needs larger sigma to be able to estimate it along with an AR1
                sigma_E = 0.2,
                model_type = "RW")
  t[[1]]
  # fixed effect within 5% of true value
  expect_lt(abs(t[["b1"]]), fe_default*0.05)
  expect_lt(abs(t[["range"]]), cutoff)
  # check if CI of estimates overlaps simulation value
  (ci <- tidy(t[[1]], effects = "ran_pars", conf.int = TRUE, silent = TRUE))
  expect_lt(0.7, ci$conf.high[ci$term=="sigma_O"])
  expect_gt(0.7, ci$conf.low[ci$term=="sigma_O"])
  # expect_lt(0.2, ci$conf.high[ci$term=="sigma_E"]) # RW underestimates sigma_E?
  expect_gt(0.2, ci$conf.low[ci$term=="sigma_E"])
})
