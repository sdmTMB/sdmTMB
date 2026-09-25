# One spatial model shared by most tests. Bias correction is only checked
# once, for get_index(), because it is slow for COG and EAO and they share its
# code path.
index_fit <- fit_once(function() {
  sdmTMB(
    density ~ 0 + as.factor(year),
    data = pcod_2011, mesh = pcod_mesh_2011,
    spatiotemporal = "off", time = "year",
    family = tweedie(link = "log")
  )
})
index_grid <- function() replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))

test_that("get_index works", {
  skip_on_cran()
  m <- index_fit()
  nd <- index_grid()

  expect_error(get_index(predict(m, newdata = nd)), regexp = "return_tmb_object")
  lifecycle::expect_deprecated(
    predictions <- predict(m, newdata = nd, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  ind <- get_index(predictions, bias_correct = FALSE)
  expect_s3_class(ind, "data.frame")
  expect_equal(get_index(m, newdata = nd, bias_correct = FALSE), ind)
  expect_equal(get_index(m, nd, bias_correct = FALSE), ind)
  expect_warning(ind_positional <- get_index(predictions, FALSE), "positional")
  expect_equal(ind_positional, ind)
  expect_error(get_index(predictions, FALSE, 0.9), "ambiguous")
  expect_equal(ind$est, c(15514.38579, 18087.5931, 21469.68763, 10818.77689),
    tolerance = 1e-5)
  expect_equal(ind$se_natural, c(2063.19087, 2506.177356, 2880.466712, 1681.308913),
    tolerance = 1e-4)
  expect_equal(ind$lwr, c(11954.6533, 13786.05969, 16505.35323, 7978.033581),
    tolerance = 1e-4)

  ind_bc <- get_index(predictions, bias_correct = TRUE)
  expect_s3_class(ind_bc, "data.frame")
  expect_equal(ind_bc$est, c(17125.71065, 19966.1714, 23699.53043, 11942.41558),
    tolerance = 1e-5)
  expect_gt(mean(ind_bc$est - ind$est), 0)

  indsp <- get_index_split(m, nd, nsplit = 2, bias_correct = FALSE)
  expect_equal(ind, indsp)
  expect_identical(chunk_time(c(1, 2, 3), 2), list(`1` = c(1, 2), `2` = 3))
  expect_error(chunk_time(c(1, 2), 0))
  expect_error(chunk_time(c(1, 2), -1))
  expect_error(chunk_time(c(1, 2), "a"))
  expect_error(chunk_time(c(1, 2), 0.2))
  expect_error(get_index_split(m, nd, nsplit = 2, predict_args = "a"), regexp = "list")

  cog <- get_cog(predictions)
  expect_s3_class(cog, "data.frame")
  expect_equal(get_cog(m, newdata = nd), cog)
  expect_equal(cog$est, c(rep(464.7814192, 4), rep(5752.411089, 4)),
    tolerance = 1e-5)
  cog_wide <- get_cog(predictions, format = "wide")
  expect_s3_class(cog_wide, "data.frame")
  expect_equal(names(cog_wide), c("year", "est_x", "lwr_x", "upr_x", "se_x",
    "est_y", "lwr_y", "upr_y", "se_y", "type"))
  expect_equal(cog$est[cog$coord == "X"], cog_wide$est_x)

  expect_error(get_index(predictions, area = c(1, 2, 3)), regexp = "area")

  # splits work with areas:
  set.seed(1)
  areas <- rlnorm(nrow(nd), meanlog = 0, sdlog = 0.1)
  ind <- get_index(predictions, area = areas, bias_correct = FALSE)
  indsp <- get_index_split(m, nd, nsplit = 2, area = areas, bias_correct = FALSE)
  expect_equal(ind, indsp)
})

test_that("get_index_sims() roughly matches get_index()", {
  skip_on_cran()
  m <- index_fit()
  nd <- index_grid()
  ind <- get_index(m, newdata = nd, bias_correct = FALSE)
  set.seed(1)
  ind_sim <- get_index_sims(predict(m, newdata = nd, nsim = 100L))
  expect_s3_class(ind_sim, "data.frame")
  expect_gt(cor(ind_sim$est, ind$est), 0.9)
  expect_gt(cor(ind_sim$lwr, ind$lwr), 0.9)
  expect_gt(cor(ind_sim$upr, ind$upr), 0.9)
})

test_that("get_index works with offsets and splits", {
  skip_on_cran()
  m2 <- sdmTMB(
    data = dogfish,
    formula = catch_weight ~ 0 + as.factor(year),
    offset = log(dogfish$area_swept),
    spatiotemporal = "off", spatial = "off",
    time = "year",
    family = tweedie(link = "log")
  )
  nd2 <- replicate_df(wcvi_grid, "year", unique(dogfish$year))
  set.seed(1)
  fake_offset <- rnorm(nrow(nd2), 0, 0.1)
  ind <- get_index(m2, newdata = nd2, offset = fake_offset, bias_correct = FALSE)
  indsp <- get_index_split(m2, nd2, nsplit = 2, offset = fake_offset, bias_correct = FALSE)
  expect_equal(ind, indsp)
  expect_error(
    get_index(m2, newdata = nd2, offset = fake_offset[-1], bias_correct = FALSE),
    "one value per row"
  )
  expect_error(
    get_index(m2, newdata = nd2, predict_args = list(offset = fake_offset),
      bias_correct = FALSE),
    "Reserved arguments"
  )
  expect_error(
    get_index(m2, newdata = nd2, predict_args = list(NA), bias_correct = FALSE),
    "must be named"
  )
  expect_warning(
    get_index_split(m2, nd2, nsplit = 2,
      predict_args = list(offset = fake_offset), bias_correct = FALSE),
    "top-level `offset`"
  )
  expect_error(
    get_index_split(m2, nd2, offset = fake_offset,
      predict_args = list(offset = fake_offset), bias_correct = FALSE),
    "not both"
  )
})

test_that("index errors are returned as needed", {
  skip_on_cran()
  g <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  expect_error(
    sdmTMB(
      density ~ 1,
      data = pcod_2011,
      spatial = "off", spatiotemporal = "off",
      family = tweedie(link = "log"),
      time = "year",
      predict_args = list(newdata = g),
      index_args = list(area = 1)
    ), regexp = "do_index" # missing!
  )

  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, spatial = "off", spatiotemporal = "off",
    family = tweedie(link = "log"),
    time = "year"
  )
  lifecycle::expect_deprecated(
    p1 <- predict(fit, newdata = NULL, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  expect_error(get_index(p1), "newdata") # missing!

  suppressMessages(
    i <- get_index(fit, newdata = g, bias_correct = FALSE)
  )
  expect_s3_class(i, "data.frame")
})

test_that("get_index() can override the derived response link for cloglog binomial models", {
  skip_on_cran()

  set.seed(1)
  n_trials <- 10
  d <- expand.grid(time = 1:2, station = seq_len(60))
  eta <- c(-1.2, -0.5)[d$time]
  p <- 1 - exp(-exp(eta))
  d$successes <- stats::rbinom(nrow(d), size = n_trials, prob = p)
  d$prop <- d$successes / n_trials
  d$trials <- n_trials

  m <- sdmTMB(
    prop ~ 0 + as.factor(time),
    data = d,
    family = binomial(link = "cloglog"),
    weights = d$trials,
    spatial = "off",
    spatiotemporal = "off",
    time = "time"
  )

  nd <- data.frame(time = 1:2)
  lifecycle::expect_deprecated(
    pred <- predict(m, newdata = nd, return_tmb_object = TRUE),
    "return_tmb_object"
  )

  idx_default <- get_index(pred, area = n_trials, bias_correct = FALSE)
  idx_log <- get_index(pred, area = n_trials, bias_correct = FALSE, derived_link = "log")
  idx_log_direct <- get_index(m, newdata = nd, area = n_trials,
    bias_correct = FALSE, derived_link = "log")
  expect_equal(idx_default$est, binomial(link = "cloglog")$linkinv(pred$data$est) * n_trials)
  expect_equal(idx_log$est, exp(pred$data$est) * n_trials)
  expect_equal(idx_log_direct, idx_log)

  idx_split <- get_index_split(
    m, nd,
    nsplit = 2,
    area = rep(n_trials, nrow(nd)),
    bias_correct = FALSE,
    derived_link = "log"
  )
  expect_equal(idx_log, idx_split)

  do_index_fit <- function(index_args) {
    sdmTMB(
      prop ~ 0 + as.factor(time),
      data = d,
      family = binomial(link = "cloglog"),
      weights = d$trials,
      spatial = "off",
      spatiotemporal = "off",
      time = "time",
      do_index = TRUE,
      predict_args = list(newdata = nd),
      index_args = index_args
    )
  }
  m_do_index <- do_index_fit(list(area = n_trials, derived_link = "log"))
  expect_equal(get_index(m_do_index, bias_correct = FALSE)$est, idx_log$est)
  m_do_index_default <- do_index_fit(list(area = n_trials))
  idx_fit_override <- get_index(m_do_index_default, bias_correct = FALSE,
    derived_link = "log")
  expect_equal(idx_fit_override$est, idx_log$est)
})

test_that("get_cog works with subsets of years", {
  skip_on_cran()
  d <- pcod_2011[pcod_2011$year %in% c(2011, 2013, 2015), , drop = FALSE]
  mesh <- make_mesh(d, c("X", "Y"), mesh = pcod_mesh_2011$mesh)

  m <- sdmTMB(
    density ~ 0 + as.factor(year),
    data = d,
    time = "year",
    spatiotemporal = "iid",
    spatial = "off",
    mesh = mesh,
    family = tweedie()
  )
  nd <- replicate_df(qcs_grid_small, "year", unique(d$year))
  nd_2011 <- replicate_df(qcs_grid_small, "year", 2011)
  nd_3 <- replicate_df(qcs_grid_small, "year", c(2015, 2011))

  # use get_weighted_average to halve time:
  cog_full <- get_weighted_average(m, newdata = nd, bias_correct = FALSE, vector = nd$X)
  expect_equal(cog_full$est, c(466.9521354, 475.6277198, 462.6024796), tolerance = 1e-5)
  cog_2011 <- get_weighted_average(m, newdata = nd_2011, bias_correct = FALSE, vector = nd_2011$X)
  cog_3 <- get_weighted_average(m, newdata = nd_3, bias_correct = FALSE, vector = nd_3$X)
  expect_equal(cog_2011$est, subset(cog_full, year == 2011)$est)
  expect_equal(cog_3$est, subset(cog_full, year %in% c(2015, 2011))$est)
})

test_that("get_index works with subsets of years", {
  skip_on_cran()

  m <- sdmTMB(
    density ~ 0 + as.factor(year),
    data = pcod_2011,
    time = "year",
    spatiotemporal = "off",
    spatial = "off",
    mesh = pcod_mesh_2011,
    family = delta_gamma()
  )
  nd <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  nd_2011 <- replicate_df(qcs_grid_small, "year", 2011)
  nd_2 <- replicate_df(qcs_grid_small, "year", c(2011, 2013))
  nd_3 <- replicate_df(qcs_grid_small, "year", c(2015, 2011))

  index_full <- get_index(m, newdata = nd, bias_correct = FALSE)
  expect_equal(index_full$est, c(19932.10781, 18159.99112, 24159.95403, 11393.82835),
    tolerance = 1e-4)
  index_2011 <- get_index(m, newdata = nd_2011, bias_correct = FALSE)
  index_2 <- get_index(m, newdata = nd_2, bias_correct = FALSE)
  index_3 <- get_index(m, newdata = nd_3, bias_correct = FALSE)
  expect_equal(index_2011$est, subset(index_full, year == 2011)$est)
  expect_equal(index_2$est, subset(index_full, year %in% c(2011, 2013))$est)
  expect_equal(index_3$est, subset(index_full, year %in% c(2015, 2011))$est)

  cog <- get_cog(m, newdata = nd, bias_correct = FALSE)
  cog2011 <- get_cog(m, newdata = nd_2011, bias_correct = FALSE)
  expect_equal(cog2011$est, cog$est[cog$year == 2011])

  eao <- get_eao(m, newdata = nd, bias_correct = FALSE)
  # no random effects, so every cell is occupied:
  expect_equal(eao$est, rep(nrow(qcs_grid_small), 4), tolerance = 1e-5)
  eao2011 <- get_eao(m, newdata = nd_2011, bias_correct = FALSE)
  expect_equal(eao2011$est, eao$est[eao$year == 2011])
})

test_that("Index integration with area vector works with extra time and possibly not all time elements in prediction data #323", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ s(depth),
    time_varying_type = "ar1",
    time_varying = ~ 1,
    time = "year",
    spatial = "off",
    spatiotemporal = "off",
    extra_time = c(2012, 2014, 2016),
    data = pcod_2011,
    family = tweedie(link = "log")
  )
  # with all years:
  nd <- replicate_df(qcs_grid_small, "year", seq(2011, 2017))
  nd$area <- 4
  ind0 <- get_index(fit, newdata = nd, area = nd$area, bias_correct = FALSE)

  # newdata doesn't have all fitted years:
  nd <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  nd$area <- 4
  ind <- get_index(fit, newdata = nd, area = nd$area, bias_correct = FALSE)
  expect_equal(ind$est - ind0$est[ind0$year %in% seq(2011, 2017, 2)], c(0, 0, 0, 0))
  expect_equal(ind$se - ind0$se[ind0$year %in% seq(2011, 2017, 2)], c(0, 0, 0, 0))
})

test_that("get_index(), get_eao(), and get_cog() take area as a vector or column name", {
  skip_on_cran()
  m <- index_fit()
  set.seed(1)
  g <- qcs_grid_small
  g$area <- runif(nrow(g), 0.9, 1.1)
  nd <- replicate_df(g, "year", unique(pcod_2011$year))

  ind <- get_index(m, newdata = nd, area = nd$area, bias_correct = FALSE)
  expect_equal(get_index(m, newdata = nd, area = "area", bias_correct = FALSE), ind)
  expect_equal(
    get_index_split(m, newdata = nd, nsplit = 2, area = "area",
      bias_correct = FALSE),
    ind
  )
  expect_equal(
    get_eao(m, newdata = nd, area = "area"),
    get_eao(m, newdata = nd, area = nd$area)
  )
  expect_equal(
    get_cog(m, newdata = nd, area = "area"),
    get_cog(m, newdata = nd, area = nd$area)
  )
})

# https://github.com/sdmTMB/sdmTMB/issues/408
test_that("Models error our nicely with Inf or -Inf covariates before get_index()", {
  d <- pcod
  d$depth_scaled[1] <- -Inf
  expect_error(m <- sdmTMB(
    data = d,
    formula = density ~ 0 + as.factor(year) + depth_scaled,
    spatiotemporal = "off", # speed
    spatial = "off", # speed
    time = "year",
    family = delta_gamma(type = "poisson-link")
  ), regexp = "Inf")
})

test_that("get_weighted_average works and matches get_cog()", {
  skip_on_cran()
  m <- index_fit()
  nd <- index_grid()
  nd$test_vector <- nd$depth

  lifecycle::expect_deprecated(
    predictions <- predict(m, newdata = nd, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  wa <- get_weighted_average(predictions, vector = nd$test_vector, bias_correct = FALSE)
  wa_direct <- get_weighted_average(m, newdata = nd, vector = nd$test_vector,
    bias_correct = FALSE)
  expect_warning(
    wa_positional <- get_weighted_average(predictions, nd$test_vector),
    "positional"
  )
  expect_s3_class(wa, "data.frame")
  expect_equal(wa_direct, wa)
  expect_equal(wa_positional, wa)
  expect_true(all(c("est", "se", "year") %in% names(wa)))
  expect_equal(nrow(wa), length(unique(pcod_2011$year)))

  # area weighting:
  set.seed(1)
  nd$area <- runif(nrow(nd), 0.9, 1.1)
  lifecycle::expect_deprecated(
    predictions_area <- predict(m, newdata = nd, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  wa_area <- get_weighted_average(predictions_area, vector = nd$test_vector,
    area = nd$area, bias_correct = FALSE)
  wa_area_direct <- get_weighted_average(m, newdata = nd,
    vector = nd$test_vector, area = "area", bias_correct = FALSE)
  expect_s3_class(wa_area, "data.frame")
  expect_equal(wa_area_direct, wa_area)

  expect_error(get_weighted_average(predictions, vector = c(1, 2, 3)), regexp = "length")
  expect_error(get_weighted_average(predictions, vector = NULL), regexp = "vector")

  # a weighted average of Y is the Y center of gravity:
  cog <- get_cog(predictions, bias_correct = FALSE, format = "wide")
  wa_lat <- get_weighted_average(predictions, vector = nd$Y, bias_correct = FALSE)
  expect_equal(wa_lat$est, cog$est_y, tolerance = 1e-6)
  expect_equal(wa_lat$se, cog$se_y, tolerance = 1e-6)
  expect_equal(wa_lat$lwr, cog$lwr_y, tolerance = 1e-6)
  expect_equal(wa_lat$upr, cog$upr_y, tolerance = 1e-6)
})

test_that("get_index() etc. errors if the data have been subset after predicting", {
  skip_on_cran()
  m <- index_fit()
  nd <- index_grid()
  lifecycle::expect_deprecated(
    p <- predict(m, newdata = nd, return_tmb_object = TRUE),
    "return_tmb_object"
  )
  p$data <- p$data[p$data$Y > 5700, , drop = FALSE]
  expect_error(get_index(p), regexp = "data")
})

test_that("index functions work directly on do_index = TRUE fits", {
  skip_on_cran()

  pcod_spde <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    spatiotemporal = "off", # speed
    time = "year", mesh = pcod_spde,
    family = tweedie(link = "log"),
    do_index = TRUE,
    predict_args = list(newdata = nd),
    index_args = list(area = 1)
  )

  # precomputed fast path still used when nothing is overridden:
  ind <- get_index(m, bias_correct = FALSE)
  ind_nd <- get_index(m, newdata = nd, bias_correct = FALSE)
  expect_equal(ind$est, ind_nd$est, tolerance = 1e-6)
  expect_equal(ind$se, ind_nd$se, tolerance = 1e-6)

  # an explicit `area` must be honoured, not silently ignored:
  ind4 <- get_index(m, area = 4, bias_correct = FALSE)
  expect_equal(ind4$est, 4 * ind$est, tolerance = 1e-6)
  ind4_nd <- get_index(m, newdata = nd, area = 4, bias_correct = FALSE)
  expect_equal(ind4$est, ind4_nd$est, tolerance = 1e-6)

  # derived quantities on the bare fit match the newdata path:
  eao <- get_eao(m, bias_correct = FALSE)
  eao_nd <- get_eao(m, newdata = nd, bias_correct = FALSE)
  expect_equal(eao$est, eao_nd$est, tolerance = 1e-6)

  wa <- get_weighted_average(m, vector = nd$depth, bias_correct = FALSE)
  wa_nd <- get_weighted_average(m, newdata = nd, vector = nd$depth,
    bias_correct = FALSE)
  expect_equal(wa$est, wa_nd$est, tolerance = 1e-6)

  cog <- get_cog(m, bias_correct = FALSE)
  cog_nd <- get_cog(m, newdata = nd, bias_correct = FALSE)
  expect_equal(cog$est, cog_nd$est, tolerance = 1e-6)
  cog_wide <- get_cog(m, bias_correct = FALSE, format = "wide")
  expect_s3_class(cog_wide, "data.frame")
  expect_true(all(c("est_x", "est_y", "year") %in% names(cog_wide)))
})
