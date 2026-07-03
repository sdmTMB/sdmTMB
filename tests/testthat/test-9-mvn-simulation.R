test_that("rmvnorm sim prediction works with no random effects", {
  skip_on_cran()
  m <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 0 + as.factor(year),
    mesh = pcod_mesh_2011, family = tweedie(link = "log"),
    spatial = "off", spatiotemporal = "off"
  )
  set.seed(1)
  nd <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  p <- predict(m, newdata = nd, nsim = 30L)
  p1 <- predict(m, newdata = nd)
  expect_identical(class(p)[[1]], "matrix")
  expect_identical(ncol(p), 30L)
  .mean <- apply(p, 1, mean)
  expect_gt(sd(p[1, , drop = TRUE]), 0.1)
  expect_gt(cor(.mean, p1$est), 0.99)

  # still works with type = "response"
  # p2 <- predict(m, newdata = nd[nd$year >= 2011, ], nsim = 30L, type = "response")
  # .mean2 <- apply(p2, 1, mean)
  # expect_gt(cor(.mean2, exp(p1$est)), 0.99)
})

test_that("rmvnorm sim prediction works", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    mesh = mesh, family = tweedie(link = "log"))
  set.seed(1)
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  p <- predict(m, newdata = nd, nsim = 15L)
  p1 <- predict(m, newdata = nd)
  expect_identical(class(p)[[1]], "matrix")
  expect_identical(ncol(p), 15L)

  # expect_equal(round(p[1:2, 1:10], 5),
  #   structure(c(1.71569, 1.83575, -1.33492, -1.27293, 0.38908, 0.70163,
  #     1.45686, 1.60475, 2.30503, 2.1425, 1.34876, 1.33194, 4.38547,
  #     4.14214, 1.29596, 1.0981, 0.50995, 0.37118, 1.85081, 1.80129), .Dim = c(2L,
  #       10L), .Dimnames = list(c("0", "0"), NULL)))

  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.99)
})

test_that("get_index_sims works", {
  skip_on_cran()
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie(link = "log"),
    time = "year", spatiotemporal = "off"
  )
  qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  set.seed(1029)
  p <- predict(m, newdata = qcs_grid_2011, nsim = 200L)
  expect_equal(ncol(p), 200L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))

  # library(dplyr)
  # a <- reshape2::melt(p)
  # a <- group_by(a, Var1, Var2) |> summarize(est = sum(exp(value)))
  # a <- group_by(a, Var1) |> summarise(est = median(est))

  x <- get_index_sims(p)
  expect_equal(nrow(x), length(unique(qcs_grid_2011$year)))
  expect_true(sum(is.na(x$se)) == 0L)

  p_regular <- predict(m, newdata = qcs_grid_2011, return_tmb_object = TRUE)
  x_regular <- get_index(p_regular)
  # expect_equal(round(x_regular$est/x$est, 5),
  #   c(0.91494, 0.92115, 0.92244, 0.9049))

  x_sims <- get_index_sims(p, return_sims = TRUE)
  expect_equal(nrow(x_sims), nrow(x) * ncol(p))

  # all areas doubled:
  xa2 <- get_index_sims(p, area = rep(2, nrow(p)))
  expect_equal(xa2$est / x$est, rep(2, 4))

  # 1 year different:
  areas <- rep(1, nrow(p))
  areas[qcs_grid_2011$year == 2011] <- 3.14
  xpi <- get_index_sims(p, area = areas)
  expect_equal(xpi$est / x$est, c(3.14, 1, 1, 1))

  # 5 cells different:
  areas[1:5] <- 0.01
  xpi2 <- get_index_sims(p, area = areas)
  expect_lt(xpi2$est[1], xpi$est[1])
  expect_equal(xpi2$est[-1], xpi$est[-1])

  # check that index still works with type = response
  expect_match(attr(p,"link"), "log")

  # p_response <- predict(m, newdata = qcs_grid_2011, nsim = 200L, type = "response")
  # expect_match(attr(p_response,"link"), "response")
  # expect_warning(get_index_sims(p_response))
  # suppressWarnings(
  #   x_response <- get_index_sims(p_response, agg_function = function(x) sum(x))
  #   )
  # x_sims2 <- get_index_sims(p)
  # # x_response$est
  # # x_sims2$est
  # expect_equal(mean(x_response$est/x_sims2$est), 1, tolerance = 0.01)
  #
  # # # check that area still works with type = response
  # suppressWarnings(
  # xpi2_response <- get_index_sims(p_response,
  #                                 area = areas,
  #                                 area_function = function(x, area) x * area,
  #                                 agg_function = function(x) sum(x ))
  # )
  # expect_equal(mean(xpi2$est[-1]/xpi2_response$est[-1]), 1, tolerance = 0.01)
  # expect_equal(xpi2$est[1]/xpi2_response$est[1], 1, tolerance = 0.05)

  # arg checking:
  expect_error(get_index_sims(3))
  expect_error(get_index_sims(level = 1.1))
  expect_error(get_index_sims(level = -0.1))
  expect_error(get_index_sims(area = c(1, 2)))
  expect_error(get_index_sims(area = rep(NA, nrow(p))))
  expect_error(get_index_sims(est_function = 3))
  expect_error(get_index_sims(agg_function = 3))

  # check that doesn't fail when attribute is missing, but does give warning
  attr(p,"link") <- NULL
  expect_warning(get_index_sims(p))
})

test_that("rmvnorm sim prediction works", {
  skip_on_cran()
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    mesh = mesh, family = tweedie(link = "log"))
  set.seed(1)
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  p <- predict(m, newdata = nd, nsim = 15L)
  p1 <- predict(m, newdata = nd)
  expect_identical(class(p)[[1]], "matrix")
  expect_identical(ncol(p), 15L)

  # expect_equal(round(p[1:2, 1:10], 5),
  #   structure(c(1.71569, 1.83575, -1.33492, -1.27293, 0.38908, 0.70163,
  #     1.45686, 1.60475, 2.30503, 2.1425, 1.34876, 1.33194, 4.38547,
  #     4.14214, 1.29596, 1.0981, 0.50995, 0.37118, 1.85081, 1.80129), .Dim = c(2L,
  #       10L), .Dimnames = list(c("0", "0"), NULL)))

  .mean <- apply(p, 1, mean)
  .sd <- apply(p, 1, sd)
  expect_gt(cor(.mean, p1$est), 0.99)
})

test_that("predict link attribute and get_index_sims work with delta", {
  skip_on_cran()
  m <- sdmTMB(density ~ 0 + as.factor(year),
              data = pcod_2011, mesh = pcod_mesh_2011, family = delta_gamma(),
              time = "year", spatiotemporal = "off"
  )
  qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  set.seed(1029)

  p <- predict(m, newdata = qcs_grid_2011, nsim = 50L, model = 1)
  expect_equal(ncol(p), 50L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))
  expect_match(attr(p,"link"), "logit")
  range(p)
  expect_lt(min(p),0)

  p <- predict(m, newdata = qcs_grid_2011, nsim = 50L, model = 2)
  expect_equal(ncol(p), 50L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))
  expect_match(attr(p,"link"), "log")
  expect_no_match(attr(p,"link"), "logit")
  range(p)
  expect_gt(min(p),0)

  # p1 <- predict(m, newdata = qcs_grid_2011, nsim = 50L, model = 1, type = "response")
  # expect_equal(ncol(p1), 50L)
  # expect_equal(nrow(p1), nrow(qcs_grid_2011))
  # expect_no_match(attr(p1,"link"), "log")
  # expect_no_match(attr(p1,"link"), "logit")
  # expect_match(attr(p1,"link"), "response")
  #
  # p2 <- predict(m, newdata = qcs_grid_2011, nsim = 50L, model = 2, type = "response")
  # expect_equal(ncol(p2), 50L)
  # expect_equal(nrow(p2), nrow(qcs_grid_2011))
  # expect_no_match(attr(p2,"link"), "log")
  # expect_no_match(attr(p2,"link"), "logit")
  # expect_match(attr(p2,"link"), "response")
  #
  # p3 <- predict(m, newdata = qcs_grid_2011, nsim = 50L, model = 1, type = "response")
  # expect_equal(ncol(p3), 50L)
  # expect_equal(nrow(p3), nrow(qcs_grid_2011))
  # expect_no_match(attr(p3,"link"), "log")
  # expect_no_match(attr(p3,"link"), "logit")
  # expect_match(attr(p3,"link"), "response")
  #
  # # check the predictions in response space are indeed in the correct ranges
  # expect_lt(max(p1),1)
  # expect_gt(min(p1),0)
  # expect_gt(min(p2),1)
  # expect_lt(min(p3),1)
  # expect_gt(min(p3),0)
  # expect_lt(max(p3),max(p2))

  p <- predict(m, newdata = qcs_grid_2011, nsim = 50L)
  expect_equal(ncol(p), 50L)
  expect_equal(nrow(p), nrow(qcs_grid_2011))
  expect_match(attr(p,"link"), "log")
  expect_no_match(attr(p,"link"), "logit")
  range(p)
  expect_lt(min(p),0)

  x <- get_index_sims(p)
  expect_equal(nrow(x), length(unique(qcs_grid_2011$year)))
  expect_true(sum(is.na(x$se)) == 0L)

  p_regular <- predict(m, newdata = qcs_grid_2011, return_tmb_object = TRUE)
  x_regular <- get_index(p_regular, bias_correct = T)

  x_sims <- get_index_sims(p, return_sims = TRUE)
  expect_equal(nrow(x_sims), nrow(x) * ncol(p))
  expect_equal(mean(x$est/x_regular$est), 1, tolerance = 0.05)
})

test_that("rmvnorm sim prediction works with various sims_vars", {
  skip_on_cran()

  # https://github.com/sdmTMB/sdmTMB/issues/107
  d <- subset(pcod, year == 2003)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 12)
  m1 <- sdmTMB(
    density ~ 1, data = d,
    mesh = pcod_spde, family = tweedie(link = "log"),
    spatial_varying = ~ 0 + depth_mean,
  )
  p1b <- predict(m1, nsim = 10, sims_var = 'zeta_s')
  expect_identical(dim(p1b), c(nrow(d), 10L))

  m2 <- sdmTMB(
    density ~ 1, data = d, control = sdmTMBcontrol(multiphase = FALSE),
    mesh = pcod_spde, family = tweedie(link = "log"),
    spatial_varying = ~ 0 + depth_mean + depth_sd,
  )
  p2b <- predict(m2, nsim = 10, sims_var = 'zeta_s')
  expect_identical(dim(p2b[[1]]), c(nrow(d), 10L))
  expect_identical(dim(p2b[[2]]), c(nrow(d), 10L))
  expect_identical(length(p2b), 2L)

  p2b <- predict(m2, nsim = 10, sims_var = 'est_rf')
  expect_identical(dim(p1b), c(nrow(d), 10L))

  p2b <- predict(m2, nsim = 10, sims_var = 'epsilon_st')
  expect_identical(dim(p1b), c(nrow(d), 10L))

  # Delta models:

  # need more data to converge:
  d <- pcod
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 10)
  m3 <- sdmTMB(
    density ~ 1, data = d,
    mesh = pcod_spde, family = delta_gamma(),
    spatial_varying = ~ 0 + depth_mean + depth_sd,
  )
  p3b <- predict(m3, nsim = 3, sims_var = 'zeta_s', model = 1)
  expect_identical(dim(p3b[[1]]), c(nrow(d), 3L))
  expect_identical(dim(p3b[[2]]), c(nrow(d), 3L))
  expect_identical(length(p3b), 2L)

  expect_warning(p3b <- predict(m3, nsim = 3, sims_var = 'zeta_s', model = NA))
  expect_identical(dim(p3b[[1]]), c(nrow(d), 3L))
  expect_identical(dim(p3b[[2]]), c(nrow(d), 3L))
  expect_identical(length(p3b), 2L)

  p3b <- predict(m3, nsim = 3, sims_var = 'zeta_s', model = 2)
  expect_identical(dim(p3b[[1]]), c(nrow(d), 3L))
  expect_identical(dim(p3b[[2]]), c(nrow(d), 3L))
  expect_identical(length(p3b), 2L)

  p3c <- predict(m3, nsim = 3, sims_var = 'omega_s', model = 1)
  expect_identical(dim(p3c), c(nrow(d), 3L))

  p3c <- predict(m3, nsim = 3, sims_var = 'epsilon_st', model = 1)
  expect_identical(dim(p3c), c(nrow(d), 3L))
})

test_that("nsim with s() and no other random effects works", {
  # https://github.com/sdmTMB/sdmTMB/issues/233
  # non-spatial model with smooth
  fit <- sdmTMB(
    density ~ s(depth),
    spatial = "off",
    data = pcod_2011,
    family = tweedie(link = "log")
  )
  p <- predict(fit, nsim = 3L)
  expect_true(ncol(p) == 3L)
})

test_that("gather/spread sims work", {
  skip_on_cran()

  m <- sdmTMB(density ~ 0 + depth_scaled + depth_scaled2,
    data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie(),
    spatiotemporal = "off")
  x <- spread_sims(m, nsim = 10)
  expect_true(nrow(x) == 10L)
  expect_s3_class(x, "data.frame")

  x <- gather_sims(m, nsim = 10)
  expect_true(ncol(x) == 3L)
  expect_s3_class(x, "data.frame")
})
