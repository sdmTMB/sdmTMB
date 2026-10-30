# eDNA standard curve (`stdcurve()`) models

sim_stdcurve <- function(seed = 1, n_plates = 20, n_obs = 400,
                         detection = c("linear", "log")) {
  detection <- match.arg(detection)
  set.seed(seed)
  plates <- paste0("p", seq_len(n_plates))
  f <- function(x) if (detection == "log") x else exp(x)
  hyper <- if (detection == "linear") {
    list(mu = c(38, -1.45, -1, 0.5), sigma = c(0.6, 0.05, 0.3, 0.1))
  } else {
    list(mu = c(38, -1.45, 1, 1.5), sigma = c(0.6, 0.05, 0.3, 0.2))
  }
  xi <- vapply(1:4, function(k) {
    stats::rnorm(n_plates, hyper$mu[k], hyper$sigma[k])
  }, numeric(n_plates))
  phi <- 0.3
  sim_ct <- function(x, j) {
    p <- stats::plogis(xi[j, 3] + xi[j, 4] * f(x))
    stats::rbinom(length(x), 1, p) *
      stats::rnorm(length(x), xi[j, 1] + xi[j, 2] * x, phi)
  }
  standards <- expand.grid(
    known_conc_ul = 10^seq(-1, 5), replicate = 1:3, plate = plates,
    stringsAsFactors = FALSE
  )
  j <- match(standards$plate, plates)
  standards$Ct <- sim_ct(log(standards$known_conc_ul), j)
  # some non-detections recorded as NA, as in lab output
  standards$Ct[standards$Ct == 0][1:3] <- NA

  predictor_dat <- data.frame(X = stats::runif(n_obs, 0, 1),
    Y = stats::runif(n_obs, 0, 1))
  mesh <- make_mesh(predictor_dat, c("X", "Y"), cutoff = 0.1)
  sim <- sdmTMB_simulate(~ 1, data = predictor_dat, mesh = mesh,
    family = gaussian(), range = 0.5, sigma_O = 0.5, phi = 0.01, B = 1,
    seed = seed)
  sim$plate <- sample(plates, n_obs, replace = TRUE)
  sim$Ct <- sim_ct(sim$eta, match(sim$plate, plates))
  list(data = sim, standards = standards, mesh = mesh, xi = xi,
    hyper = hyper, phi = phi, plates = plates)
}

test_that("stdcurve() validates the standards", {
  std <- data.frame(plate = c("a", "b"), Ct = c(30, 0), known_conc_ul = c(1, 2))
  expect_s3_class(stdcurve(std), "family")
  expect_identical(stdcurve(std)$family, "stdcurve")
  expect_error(stdcurve(), "data frame")
  expect_error(stdcurve(std[, -1]), "plate")
  expect_error(stdcurve(std[, -2]), "Ct")
  expect_error(stdcurve(std[, -3]), "known_conc_ul")
  std$known_conc_ul[1] <- 0
  expect_error(stdcurve(std), "> 0")
  std$known_conc_ul[1] <- 1
  std$Ct[1] <- -1
  expect_error(stdcurve(std), ">= 0")
  names(std) <- c("pl", "ct", "conc")
  std$ct[1] <- 1
  expect_s3_class(stdcurve(std, plate = "pl", ct = "ct", conc = "conc"),
    "family")
})

test_that("stdcurve models recover simulated standard curves", {
  skip_on_cran()
  local_rtmb_backend()
  s <- sim_stdcurve()
  fit <- sdmTMB(Ct ~ 1, data = s$data, mesh = s$mesh,
    family = stdcurve(standards = s$standards))
  expect_true(all(unlist(sanity(fit, silent = TRUE)[c("hessian_ok",
    "eigen_values_ok", "nlminb_ok")])))

  # latent log concentration intercept
  b <- tidy(fit)
  expect_equal(b$estimate[b$term == "(Intercept)"], 1, tolerance = 0.15)
  expect_equal(sigma(fit), s$phi, tolerance = 0.1)

  sc <- get_stdcurve(fit)
  expect_named(sc, c("plates", "hyper"))
  mu <- sc$hyper$estimate[sc$hyper$parameter == "mu"]
  expect_equal(mu[1:2], s$hyper$mu[1:2], tolerance = 0.02)
  expect_equal(mu[3:4], s$hyper$mu[3:4], tolerance = 0.3)
  est <- matrix(sc$plates$estimate, ncol = 4)
  est <- est[match(s$plates, unique(sc$plates$plate)), ] # sorted plate names
  expect_gt(stats::cor(est[, 1], s$xi[, 1]), 0.95)
  expect_gt(stats::cor(est[, 2], s$xi[, 2]), 0.9)
  expect_gt(stats::cor(est[, 3], s$xi[, 3]), 0.5)
  expect_identical(unique(sc$plates$plate), sort(s$plates))

  # post-fit methods
  p <- predict(fit)
  expect_equal(stats::cor(p$est, s$data$eta), 1, tolerance = 0.1)
  nd <- data.frame(X = c(0.2, 0.5), Y = c(0.3, 0.5))
  expect_equal(nrow(predict(fit, newdata = nd)), 2L)
  r <- residuals(fit, type = "mle-eb")
  expect_length(r, nrow(s$data))
  expect_true(all(is.finite(r)))
  # response residuals and fitted() are on the Ct scale
  pars <- get_pars(fit)
  xi <- pars$std_xi[fit$tmb_data$stdcurve$plate_i, ]
  ct_hat <- stats::plogis(xi[, 3] + xi[, 4] * exp(p$est)) *
    (xi[, 1] + xi[, 2] * p$est)
  expect_equal(fitted(fit), ct_hat)
  expect_equal(residuals(fit, type = "response"), s$data$Ct - ct_hat)
  expect_equal(mean(residuals(fit, type = "response")), 0, tolerance = 0.5)
  expect_error(residuals(fit, type = "pearson"), "Pearson")
  # default mle-mvn: plate coefficients come from the same posterior draw
  # as the linear predictor
  set.seed(1)
  r_mvn <- residuals(fit)
  expect_true(all(is.finite(r_mvn)))
  expect_equal(stats::sd(r_mvn), 1, tolerance = 0.15)
  set.seed(1)
  draw <- .one_sample_posterior(fit)
  mu <- predict(fit, newdata = fit$data, mcmc_samples = matrix(draw, ncol = 1L),
    nsim = 1L, offset = fit$offset)[, 1L]
  y <- fit$response[, 1L]
  r_draw <- qres_stdcurve(fit, y, mu, .par = draw)
  expect_equal(r_mvn, r_draw)
  set.seed(1)
  .one_sample_posterior(fit) # same RNG state as above
  r_mode <- qres_stdcurve(fit, y, mu)
  expect_false(isTRUE(all.equal(r_draw, r_mode)))
  set.seed(1)
  sims <- simulate(fit, nsim = 2)
  expect_equal(dim(sims), c(nrow(s$data), 2L))
  expect_equal(mean(sims[, 1] > 0), mean(s$data$Ct > 0), tolerance = 0.15)
  expect_output(print(fit), "stdcurve")

  # simulate(newdata = ) uses each new row's own plate
  nd <- s$data[c(10:1, 1:5), ] # reordered, with repeats
  set.seed(3)
  sims_nd <- simulate(fit, newdata = nd, nsim = 200, type = "mle-eb")
  expect_equal(dim(sims_nd), c(nrow(nd), 200L))
  pars <- get_pars(fit)
  j <- match(nd$plate, fit$tmb_data$stdcurve$plates)
  eta <- predict(fit, newdata = nd)$est
  pos_mean <- pars$std_xi[j, 1] + pars$std_xi[j, 2] * eta
  sim_pos_mean <- apply(sims_nd, 1, function(x) mean(x[x > 0]))
  expect_lt(max(abs(sim_pos_mean - pos_mean)), 0.15)
  expect_error(simulate(fit, newdata = nd[, names(nd) != "plate"], nsim = 1),
    "plate")
  nd$plate[1] <- "unknown"
  expect_error(simulate(fit, newdata = nd, nsim = 1), "standards")
  # concentration predictions don't need plates
  expect_equal(nrow(predict(fit, newdata = nd[, c("X", "Y")])), nrow(nd))

  # fixing a plate SD with map/start
  fit2 <- update(fit, control = sdmTMBcontrol(
    map = list(ln_std_sigma = factor(c(1, 2, NA, 3))),
    start = list(ln_std_sigma = c(0, -2, log(0.01), -2))))
  h <- get_stdcurve(fit2)$hyper
  expect_equal(h$estimate[7], 0.01)
  expect_true(is.na(h$std.error[7]))
  expect_false(anyNA(h$std.error[-7]))
})

test_that("stdcurve detection on the log scale works", {
  skip_on_cran()
  local_rtmb_backend()
  s <- sim_stdcurve(seed = 2, detection = "log")
  fit <- sdmTMB(Ct ~ 1, data = s$data, mesh = s$mesh,
    family = stdcurve(standards = s$standards, detection = "log"))
  b <- tidy(fit)
  expect_equal(b$estimate[b$term == "(Intercept)"], 1, tolerance = 0.15)
  sc <- get_stdcurve(fit)
  mu <- sc$hyper$estimate[sc$hyper$parameter == "mu"]
  expect_equal(mu, s$hyper$mu, tolerance = 0.3)
})

test_that("stdcurve cross-validation predicts Ct", {
  skip_on_cran()
  local_rtmb_backend()
  s <- sim_stdcurve(n_plates = 6, n_obs = 200)
  mesh <- make_mesh(s$data, c("X", "Y"), cutoff = 0.15)
  set.seed(1)
  cv <- sdmTMB_cv(Ct ~ 1, data = s$data, mesh = mesh,
    family = stdcurve(standards = s$standards), k_folds = 2)
  expect_true(all(cv$converged))
  d <- cv$data
  # fold predictions are expected Ct on each row's plate, not log concentration
  for (k in 1:2) {
    fit_k <- cv$models[[k]]
    rows <- d$cv_fold == k
    expect_equal(nrow(fit_k$data), nrow(d)) # all rows, held-out at weight 0
    expect_equal(d$cv_predicted[rows], fitted(fit_k)[rows])
  }
  expect_gt(min(d$cv_predicted), 5)
  expect_lt(abs(mean(d$Ct - d$cv_predicted)), 2)
  expect_true(is.finite(cv$sum_loglik))
})

test_that("stdcurve models check their inputs", {
  skip_on_cran()
  local_rtmb_backend()
  s <- sim_stdcurve(n_plates = 4, n_obs = 100)
  fam <- stdcurve(standards = s$standards)
  mesh <- make_mesh(s$data, c("X", "Y"), cutoff = 0.2)
  expect_error(sdmTMB(Ct ~ 1, data = s$data[, names(s$data) != "plate"],
    mesh = mesh, family = fam), "plate")
  d <- s$data
  d$plate[1] <- "unknown"
  expect_error(sdmTMB(Ct ~ 1, data = d, mesh = mesh, family = fam),
    "standards")
  d <- s$data
  d$Ct[1] <- -1
  expect_error(sdmTMB(Ct ~ 1, data = d, mesh = mesh, family = fam), ">= 0")
  expect_error(sdmTMB(Ct ~ 1, data = s$data, mesh = mesh, family = fam,
    control = sdmTMBcontrol(backend = "tmb")), "rtmb")
  one <- stdcurve(standards = s$standards[s$standards$plate == "p1", ])
  d <- s$data[s$data$plate == "p1", ]
  expect_error(sdmTMB(Ct ~ 1, data = d, spatial = "off", family = one),
    "at least 2 plates")
  expect_error(sdmTMB(Ct ~ 1, data = s$data, mesh = mesh, family = fam,
    dispformula = ~ X), "dispformula")
  expect_error(get_stdcurve(sdmTMB(Ct ~ 1, data = s$data, spatial = "off")),
    "stdcurve")
})
