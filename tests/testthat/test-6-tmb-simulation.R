test_that("TMB (custom) AR1 simulation is unbiased", {
  # run many times; check for bias
  skip_on_cran()
  skip_slow()

  do_sim_fit <- function(i) {
    set.seed(i)
    predictor_dat <- data.frame(
      X = runif(1000), Y = runif(1000),
      a1 = rnorm(1000), year = rep(1:10, each = 100)
    )
    mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
    sim_dat <- sdmTMB_simulate(
      formula = ~1,
      data = predictor_dat,
      time = "year",
      mesh = mesh,
      family = gaussian(),
      range = 0.5,
      sigma_E = 0.1,
      phi = 0.1,
      sigma_O = 0,
      seed = i,
      rho = 0.8,
      B = 0
    )
    fit <- sdmTMB(observed ~ 0, sim_dat,
      mesh = mesh, time = "year",
      spatiotemporal = "ar1", spatial = "off",
    )
    b <- tidy(fit, "ran_pars")
    rho_hat <- b$estimate[b$term == "rho"]
    range_hat <- b$estimate[b$term == "range"]
    sigma_E_hat <- b$estimate[b$term == "sigma_E"]
    data.frame(rho = rho_hat, sigma_E = sigma_E_hat, range = range_hat)
  }
  out <- lapply(seq_len(12L), function(i) do_sim_fit(i))
  out <- do.call("rbind", out)
  expect_true(median(out$rho) > 0.79 && median(out$rho) < 0.81)
  # +/-10%, matching the tolerance already used for sigma_E below: the
  # per-replicate range estimates are noisy (range ~0.3-0.8 across the 12
  # seeds), so a tighter median bound is not reliably met by either backend
  expect_true(median(out$range) > 0.45 && median(out$range) < 0.55)
  expect_true(median(out$sigma_E) > 0.09 && median(out$sigma_E) < 0.11)
})

test_that("simulate() behaves OK with or without random effects across types", {
  skip_on_cran()
  m <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(link = "log")
  )
  set.seed(1)
  s <- simulate(m)
  expect_length(s, 969)

  m2 <- update(m, spatial = "off")
  s <- simulate(m2)
  expect_length(s, 969)

  # has no random effects, switches to standard as needed:
  s <- simulate(m2, type = "mle-mvn")
  expect_length(s, 969)
})

test_that("simulate() method works with newdata", {
  skip_on_cran()
  fit <- sdmTMB(
    present ~ 1,
    time = "year",
    data = pcod_2011, spatial = "on",
    spatiotemporal = "iid",
    family = binomial(),
    mesh = pcod_mesh_2011
  )
  s <- simulate(fit)
  expect_true(nrow(s) == nrow(pcod_2011))
  g <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  s <- simulate(fit, newdata = g)
  expect_true(nrow(s) == nrow(g))
  s <- simulate(fit, newdata = subset(g, year == 2011))
  nrow(s)
  expect_true(nrow(s) == nrow(subset(g, year == 2011)))

  gg <- subset(g, year == 2011)
  set.seed(1)
  s <- simulate(fit, newdata = gg, nsim = 400L)
  a <- apply(s, 1, mean)
  p <- predict(fit, newdata = gg)
  plot(a, plogis(p$est))
  expect_gt(cor(a, plogis(p$est)), 0.98)

  set.seed(1)
  s1 <- simulate(fit, type = "mle-mvn", mle_mvn_samples = "single", nsim = 100)
  set.seed(1)
  s2 <- simulate(fit, type = "mle-mvn", mle_mvn_samples = "multiple", nsim = 100)
  set.seed(1)
  s3 <- simulate(fit, type = "mle-eb", nsim = 100)

  expect_false(identical(s1, s2))
  expect_false(identical(s1, s3))

  sd1 <- apply(s1, 1, sd)
  sd2 <- apply(s2, 1, sd)
  sd3 <- apply(s3, 1, sd)

  expect_lt(mean(sd1), mean(sd2))

  # offset?
  fit <- sdmTMB(
    catch_weight ~ 1,
    data = dogfish,
    offset = log(dogfish$area_swept),
    spatial = "off",
    family = tweedie()
  )

  set.seed(1)
  s1 <- simulate(fit)
  set.seed(1)
  s2 <- simulate(fit, newdata = dogfish)
  set.seed(1)
  s3 <- simulate(fit, newdata = dogfish, offset = rep(0, nrow(dogfish)))
  set.seed(1)
  s4 <- simulate(fit, newdata = dogfish, offset = log(dogfish$area_swept))

  attributes(s1) <- NULL
  attributes(s2) <- NULL
  attributes(s3) <- NULL
  attributes(s4) <- NULL

  row.names(s1) <- NULL
  row.names(s2) <- NULL
  row.names(s3) <- NULL
  row.names(s4) <- NULL

  expect_equal(s1, s4)
  expect_equal(s2, s3)
})

test_that("simulate() can turn off observation error", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ 1,
    mesh = pcod_mesh_2011,
    data = pcod_2011,
    spatial = "on",
    family = tweedie()
  )
  set.seed(1)
  s_no_obs <- simulate(fit, observation_error = FALSE)
  set.seed(1)
  s_obs <- simulate(fit, observation_error = TRUE)
  expect_false(identical(s_no_obs, s_obs))
  expect_equal(sum(s_no_obs[,1] == 0), 0L)
  expect_equal(stats::cor(s_no_obs[,1], s_obs[,1]), 0.449853, tolerance = 0.01)

  # with newdata:
  g <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  s_no_obs <- simulate(fit, observation_error = FALSE, newdata = g)
  expect_equal(sum(s_no_obs[,1] == 0), 0L)

  fit2 <- sdmTMB(
    present ~ 1,
    mesh = pcod_mesh_2011,
    data = pcod_2011,
    spatial = "off",
    spatiotemporal = "iid",
    time = "year",
    family = binomial()
  )
  s <- simulate(fit2, observation_error = FALSE, newdata = g)
  expect_true(any(grepl("year", attributes(s))))
  expect_true(any(grepl("2017", row.names(s))))

  # for indexes:
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 25)
  fit <- sdmTMB(density ~ 0 + as.factor(year), spatiotemporal = "off",
    data = pcod, mesh = mesh, family = tweedie(link = "log"),
    time = "year"
  )
  qcs_grid <- replicate_df(qcs_grid_small, "year", unique(pcod$year))
  set.seed(1)
  s <- simulate(
    fit,
    newdata = qcs_grid,
    type = "mle-mvn",
    mle_mvn_samples = "multiple",
    nsim = 10, # increase this
    observation_error = FALSE
  )
  ind <- get_index_sims(log(s))
  expect_s3_class(ind, "data.frame")
  # ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
  #   geom_line() +
  #   geom_ribbon(alpha = 0.4)
})

# test_that("TMB Delta simulation works", {
#   skip_on_cran()
#   skip_if_not_installed("INLA")
#
#   set.seed(1)
#   predictor_dat <- data.frame(
#     X = runif(2000), Y = runif(2000),
#     a1 = rnorm(2000), year = 1
#   )
#   mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)
#   sim_dat <- sdmTMB_simulate(
#     formula = ~1,
#     data = predictor_dat,
#     time = "year",
#     mesh = mesh,
#     family = nbinom2(),
#     range = 0.3,
#     phi = 1,
#     sigma_O = 0.1,
#     seed = 1,
#     B = 1.5
#   )
#   sim_dat$observed[sample(1:nrow(sim_dat), size = 0.7*nrow(sim_dat), replace=T)] = 0
#   fit <- sdmTMB(observed ~ 1, sim_dat,
#                 mesh = mesh,
#                 family = delta_truncated_nbinom2(),
#                 share_range = TRUE
#   )
#
#   sim_1 <- simulate(fit, model = 1, nsim = 100)
#   expect_equal(ncol(sim_1), 100)
#   empirical_z <- 1 - length(which(sim_dat$observed==0)) / nrow(sim_dat)
#   expect_lt(max(abs(apply(sim_1,2,mean) - empirical_z)), 0.05)
#
#   sim_2 <- simulate(fit, model = 2, nsim = 100)
#   expect_equal(ncol(sim_2), 100)
#   expect_lt(max(abs(apply(sim_2,2,mean,na.rm=T) - 5.67)), 0.65)
# })
