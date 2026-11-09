rtmb_delta_data <- function(n = 90L, seed = 71) {
  set.seed(seed)
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n), w = runif(n),
    time = rep(1:3, each = n / 3), g = factor(rep(1:6, length.out = n)),
    off = rnorm(n, sd = 0.1), positive = rgamma(n, 2, 2))
  d$response <- ifelse(runif(n) < 0.6, rgamma(n, 2, 1), 0)
  d$present <- as.integer(d$response > 0)
  d
}

test_that("RTMB binomial, Gamma, and lognormal likelihoods match frozen TMB values", {
  d <- rtmb_delta_data()
  d$positive[[4L]] <- NA
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
  cases <- list(
    list(present ~ z, binomial(), "on"),
    list(present ~ z, binomial(link = "cloglog"), "off"),
    list(positive ~ z, Gamma(link = "log"), "on"),
    list(positive ~ z, lognormal(), "on")
  )
  for (case in cases) {
    fit <- sdmTMB(case[[1L]], data = d, mesh = mesh, family = case[[2L]],
      spatial = case[[3L]], weights = rep(c(1, 0.5), 45L), do_fit = FALSE)
    expect_frozen_fit_data("delta", fit, d[c(1, 40, 90), ],
      info = paste(case[[2L]]$family, case[[2L]]$link), project = case[[3L]] == "on")
  }
  # The inverse link needs a positive linear predictor.
  fit <- sdmTMB(positive ~ z, data = d, mesh = mesh, spatial = "off",
    family = Gamma(link = "inverse"), do_fit = FALSE)
  expect_frozen_fit_data("delta", fit, d, info = "inverse", project = FALSE,
    edit = function(p) {
      p$b_j[] <- c(1, 0.1)
      p
    })
})

test_that("RTMB delta models match frozen TMB values", {
  skip_on_cran()
  d <- rtmb_delta_data()
  d$response[[5L]] <- NA
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
  newdata <- d[c(1, 40, 90), ]
  cases <- list(
    delta_gamma = list(response ~ z + s(w, k = 4) + (1 | g),
      family = delta_gamma(), time = "time", spatiotemporal = "ar1",
      offset = "off"),
    delta_lognormal = list(response ~ z, family = delta_lognormal(),
      time = "time", spatial = list("on", "off"),
      spatiotemporal = list("iid", "off"), weights = runif(90L)),
    poisson_link = list(response ~ z, offset = "off",
      family = delta_gamma(type = "poisson-link"),
      spatial = list("on", "off")),
    dispersion = list(response ~ z, dispformula = ~w,
      family = delta_gamma(), spatial = "off"),
    priors = list(response ~ z, family = delta_lognormal(), time = "time",
      spatiotemporal = "rw", bayesian = TRUE,
      priors = sdmTMBpriors(b = normal(c(0, 0), c(1, 1)),
        matern_s = pc_matern(1, 1), matern_st = pc_matern(1, 1),
        phi = halfnormal(0, 1))),
    reml = list(response ~ z + s(w, k = 4), family = delta_gamma(),
      reml = TRUE),
    time_varying = list(response ~ 1, family = delta_gamma(),
      time = "time", time_varying = ~ 0 + z, spatial = "off",
      priors = sdmTMBpriors(sigma_V = gamma_cv(0.5, 0.5))),
    svc_threshold = list(response ~ breakpt(z), family = delta_gamma(),
      spatial_varying = ~ 0 + w),
    svc_partial = list(response ~ 1, family = delta_gamma(),
      spatial = list("on", "off"), spatial_varying = ~ 0 + w)
  )
  for (name in names(cases)) {
    fit <- do.call(sdmTMB, c(cases[[name]],
      list(data = d, mesh = mesh, do_fit = FALSE)))
    expect_frozen_fit_data("delta", fit, newdata, info = name,
      project = name %in% c("delta_gamma", "delta_lognormal",
        "poisson_link", "svc_threshold", "svc_partial"))
  }
})

test_that("Delta SVC fields have a density without a spatial intercept field", {
  d <- rtmb_delta_data()
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
  fit <- sdmTMB(response ~ 1, data = d, mesh = mesh, family = delta_gamma(),
    spatial = list("on", "off"), spatial_varying = ~ 0 + w, do_fit = FALSE)
  expect_false(anyNA(fit$tmb_map$ln_kappa[3:4]))
  obj <- fit$tmb_obj
  par <- obj$par
  nll <- obj$fn(par)
  par[names(par) == "ln_tau_Z"][2L] <- 1
  expect_false(isTRUE(all.equal(obj$fn(par), nll)))
})

test_that("RTMB multi-family models match frozen TMB values", {
  d <- rtmb_delta_data()
  d$dist <- rep(c("g", "d", "p"), 30L)
  d$y2 <- ifelse(d$dist == "p", rpois(90L, 2),
    ifelse(d$dist == "g", rnorm(90L), d$response))
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
  newdata <- d[c(1, 2, 3, 90), ]
  fit <- sdmTMB(y2 ~ z, data = d, mesh = mesh, time = "time",
    spatiotemporal = "iid", offset = "off", distribution_column = "dist",
    family = list(g = gaussian(), d = delta_gamma(type = "poisson-link"),
      p = poisson()), do_fit = FALSE)
  expect_frozen_fit_data("delta", fit, newdata, info = "poisson-link mix")
  fit <- sdmTMB(y2 ~ z, data = d, mesh = mesh, distribution_column = "dist",
    family = list(g = gaussian(), d = delta_lognormal(), p = nbinom2()),
    do_fit = FALSE)
  expect_frozen_fit_data("delta", fit, newdata, info = "delta mix")
})

test_that("RTMB delta fits give finite SEs, predictions, and draws", {
  skip_on_cran()
  set.seed(72)
  n <- 150L
  d <- data.frame(x = runif(n), y = runif(n), z = rnorm(n))
  mu <- exp(0.5 + 0.3 * d$z + cos(3 * d$y))
  d$response <- ifelse(runif(n) < plogis(0.3 + 0.5 * d$z),
    rgamma(n, 2, 2 / mu), 0)
  mesh <- make_mesh(d, c("x", "y"), n_knots = 30L, type = "kmeans")
  for (family in list(delta_gamma(), delta_gamma(type = "poisson-link"))) {
    rtmb <- sdmTMB(response ~ z, data = d, mesh = mesh, family = family,
      spatial = list("off", "on"),
      control = sdmTMBcontrol(multiphase = FALSE))
    info <- family$type
    se <- summary(rtmb$sd_report, "fixed")[, "Std. Error"]
    expect_true(all(is.finite(se)), info = info)
    actual <- predict(rtmb, newdata = d[1:6, ], se_fit = TRUE)
    expect_true(all(is.finite(actual$est_se)), info = info)
    set.seed(1)
    draws <- simulate(rtmb, nsim = 300)
    expect_equal(dim(draws), c(n, 300L))
    expect_equal(mean(draws == 0), mean(d$response == 0), tolerance = 0.1,
      info = info)
    expect_equal(rowMeans(draws), predict(rtmb, type = "response")$est,
      tolerance = 0.15, info = info)
  }
})

test_that("RTMB reports response means when observations are not simulated", {
  d <- rtmb_delta_data()
  mesh <- make_mesh(d, c("x", "y"), n_knots = 25L, type = "kmeans")
  cases <- list(
    list(response ~ z, family = delta_gamma(type = "poisson-link")),
    list(response ~ z, family = delta_lognormal()),
    list(present ~ z, family = binomial())
  )
  for (case in cases) {
    fit <- do.call(sdmTMB, c(case, list(data = d, mesh = mesh,
      do_fit = FALSE)))
    p <- rtmb_test_parameters(fit$tmb_params, fit$tmb_map)
    data <- fit$tmb_data
    data$sim_obs <- 0L
    expect_frozen_values(paste("delta response means",
      paste(case$family$family, collapse = "-"), case$family$type),
      as.vector(make_sdmTMB_adfun(data, p, fit$tmb_map,
        fit$tmb_random)$simulate()$y_i), tolerance = 1e-7)
  }
})
