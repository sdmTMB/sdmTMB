test_that("Offset works", {
  skip_on_cran()

  pcod$offset <- rnorm(nrow(pcod))
  fit2 <- sdmTMB(density ~ 1,
    offset = pcod$offset,
    data = pcod, spatial = "off",
    family = tweedie()
  )
  expect_error(fit2 <- sdmTMB(density ~ 1,
    offset = year,
    data = pcod, spatial = "off",
    family = tweedie()
  ), regexp = "year")
})


test_that("Offset matches glmmTMB", {
  skip_on_cran()

  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  pcod_pos <- subset(pcod, density > 0)

  fit1 <- glmmTMB::glmmTMB(density ~ 1,
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  fit1_off <- glmmTMB::glmmTMB(density ~ 1 + offset(offset),
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  pcod_pos$offset2 <- log(1)
  fit1_off0 <- glmmTMB::glmmTMB(density ~ 1 + offset(offset2),
    data = pcod_pos,
    family = Gamma(link = "log")
  )

  fit2 <- sdmTMB(density ~ 1,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  fit2_off <- sdmTMB(density ~ 1,
    offset = pcod_pos$offset,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  fit2_off0 <- sdmTMB(density ~ 1,
    offset = pcod_pos$offset2,
    data = pcod_pos, spatial = "off",
    family = Gamma(link = "log")
  )

  b1 <- summary(fit1)$coefficients$cond[1]
  b1_offset <- summary(fit1_off)$coefficients$cond[1]
  b1_offset0 <- summary(fit1_off0)$coefficients$cond[1]

  b2 <- tidy(fit2)$estimate[1]
  b2_offset <- tidy(fit2_off)$estimate[1]
  b2_offset0 <- tidy(fit2_off0)$estimate[1]

  # test that glmmTMB and sdmTMB agree
  expect_equal(b2, b1, tolerance = 1e-4)
  expect_equal(b2_offset, b1_offset, tolerance = 1e-4)
  expect_equal(b2_offset0, b1_offset0, tolerance = 1e-4)

  # the offset of 0 is same as no offset
  expect_equal(b2_offset0, b2, tolerance = 1e-8)

  # the offset is doing something
  expect_false(((b2_offset - b2) == 0))
})

test_that("Offset works with extra_time", {
  skip_on_cran()
  set.seed(1)
  pcod$offset <- rnorm(nrow(pcod))
  mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), n_knots = 80)
  fit <- sdmTMB(density ~ 1,
    offset = pcod$offset,
    mesh = mesh,
    time = "year",
    extra_time = c(2006L, 2008L, 2010L, 2012L, 2014L, 2016L),
    data = pcod, spatial = "off",
    spatiotemporal = "ar1",
    family = tweedie()
  )
  expect_true(inherits(fit, "sdmTMB"))
  b <- tidy(fit, "ran_pars")
  expect_equal(round(b$estimate[b$term == "rho"], 2), 0.91)
})

test_that("Offset prediction matches glm()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  dat <- pcod[pcod$density > 0,]
  dat$.offset <- rnorm(nrow(dat))
  fit <- sdmTMB(
    present ~ 1,
    offset = dat$.offset,
    data = dat, spatial = "off",
    family = Gamma(link = "log")
  )
  fit_glm <- glm(
    present ~ 1 + offset(.offset),
    data = dat,
    family = Gamma(link = "log")
  )
  fit_glmmTMB <- glmmTMB::glmmTMB(
    present ~ 1 + offset(.offset),
    data = dat,
    family = Gamma(link = "log")
  )

  p <- predict(fit)
  p_glm <- predict(fit_glm)
  p_glmmTMB <- predict(fit_glmmTMB)
  expect_equal(p$est, unname(p_glm))
  expect_equal(p$est, p_glmmTMB)

  p_glmmTMB <- predict(fit_glmmTMB, newdata = dat)
  expect_equal(p$est, unname(p_glm))
  expect_equal(p$est, p_glmmTMB)

  set.seed(1)
  p <- predict(fit, nsim = 1000, offset = dat$.offset)
  mu <- apply(p, 1, mean)
  plot(mu, p_glm)
  expect_equal(unname(mu), unname(p_glm), tolerance = 0.01)

  # sdmTMB ignores offset here (but not glm() or glmmTMB()!)
  # p <- predict(fit, newdata = dat)
  # p_glmmTMB <- predict(fit_glmmTMB, newdata = dat)
  # expect_equal(p$est, unname(p_glmmTMB))
})

test_that("offset gets passed through cross validation as expected #372", {
  skip_on_cran()
  dat <- subset(dogfish, catch_weight > 0)
  expect_error(
    x <- sdmTMB_cv(catch_weight ~ 1,
      data = dat,
      family = Gamma("log"), offset = log(dat$area_swept), spatial = "off"
    ),
    "offset"
  )
  set.seed(1)
  x <- sdmTMB_cv(catch_weight ~ 1,
    data = dat, family = Gamma("log"),
    offset = "area_swept", spatial = "off",
    mesh = make_mesh(dat, c("X", "Y"), cutoff = 10), k_folds = 2
  )
  y <- x$data[, c("catch_weight", "cv_predicted")]
  # plot(y$catch_weight, y$cv_predicted)
  # if offset is applied, will have unique values because an intercept-only model:
  expect_true(length(unique(y$cv_predicted)) >= length(unique(dat$area_swept)))
})

test_that("predicting on newdata with a non-null offset in fit but a null offset in predict informs the user appropriately", {
  skip_on_cran()
  dat <- subset(dogfish, catch_weight > 0)
  fit <- sdmTMB(
    catch_weight ~ 1,
    data = dat,
    family = Gamma("log"),
    offset = "area_swept",
    spatial = "off"
  )
  rlang::local_options(rlib_message_verbosity = "verbose") # ignore "once"
  expect_no_message(pred <- predict(fit))
  expect_no_message(pred <- predict(fit, offset = rep(0, nrow(dat))))
  expect_no_message(pred0 <- predict(fit, newdata = qcs_grid, offset = 0))
  expect_message({pred <- predict(fit, newdata = qcs_grid)}, regexp = "offset")
  expect_equal(pred$est, pred0$est)

  # a scalar offset is recycled
  pred1 <- predict(fit, newdata = qcs_grid, offset = log(100))
  expect_equal(pred1$est, pred0$est + log(100))
  expect_error(predict(fit, newdata = qcs_grid, offset = c(0, 1)), regexp = "length one")

  # get_index() defaults to an offset of 0 without a message
  nd <- replicate_df(qcs_grid, "year", unique(dat$year))
  fit_t <- sdmTMB(
    catch_weight ~ 0 + as.factor(year),
    data = dat,
    family = Gamma("log"),
    offset = "area_swept",
    time = "year",
    spatial = "off",
    spatiotemporal = "off"
  )
  expect_no_message(
    ind <- get_index(fit_t, newdata = nd, bias_correct = FALSE),
    message = "offset"
  )
  ind0 <- get_index(fit_t, newdata = nd, offset = 0, bias_correct = FALSE)
  expect_equal(ind, ind0)
  ind1 <- get_index(fit_t, newdata = nd, offset = log(2), bias_correct = FALSE)
  expect_equal(ind1$est, 2 * ind0$est)
})

# #
# # offset/prediction setting checks:
#
# pos <- dogfish[dogfish$catch_weight > 0,]
#
# m1 <- sdmTMB(catch_weight ~ 1, family = Gamma("log"), data = pos, offset = log(pos$area_swept), spatial = "off")
# # m2 <- glmmTMB::glmmTMB(catch_weight ~ 1, family = Gamma("log"), data = pos, offset = log(pos$area_swept))
# # m3 <- glm(catch_weight ~ 1, family = Gamma("log"), data = pos, offset = log(pos$area_swept))
#
# head(predict(m1, newdata = pos, offset = rep(0, nrow(pos)))$est) # right
# head(predict(m1, offset = rep(0, nrow(pos)))$est) # right (was wrong)
#
# head(predict(m1, offset = log(pos$area_swept))$est) # right
# head(predict(m1, newdata = pos, offset = log(pos$area_swept))$est) # right
#
# head(predict(m1)$est) # right
# head(predict(m1, newdata = pos)$est) # right
#
# head(predict(m1, newdata = pos, offset = rep(0, nrow(pos)), nsim = 2)) # right
# head(predict(m1, offset = rep(0, nrow(pos)), nsim = 2)) # right (was wrong)
#
# head(predict(m1, offset = log(pos$area_swept), nsim = 2)) # right
# head(predict(m1, newdata = pos, offset = log(pos$area_swept), nsim = 2)) # right
#
# head(predict(m1, nsim = 2)) # right
# head(predict(m1, newdata = pos, nsim = 2)) # right
#
#
# # m2 <- glmmTMB::glmmTMB(catch_weight ~ 1, family = Gamma("log"), data = pos)
# # head(predict(m2))
# # head(predict(m2, newdata = pos))
# # predict(m2, newdata = pos[1:3,,drop=FALSE])
# #
# # e1 <- predict(m1, newdata = pos, offset = rep(0, nrow(pos)))$est
# # plot(exp(e1), pos$catch_weight)
# # mean(exp(e1))
# # mean(pos$catch_weight)
# #
# #
# # e1 <- predict(m1, newdata = pos, offset = log(pos$area_swept))$est
# # plot(exp(e1), pos$catch_weight)
# # mean(exp(e1))
# # mean(pos$catch_weight)
# #
# # head(predict(m1, newdata = pos, offset = rep(0, nrow(pos)))$est)
# #
# # head(predict(m2))
# # head(predict(m2, newdata = pos))
# #
# # head(predict(m2))
# # head(predict(m2, newdata = pos))
# # head(predict(m3))
# #
# # expect_equal(p7[,1,drop=TRUE], p6[,1,drop=TRUE])
# #
# # set.seed(1)
# # suppressWarnings(p8 <- predict(m, newdata = pcod, nsim = 2L))
# # suppressWarnings(p9 <- predict(m, newdata = pcod, offset = rep(0, nrow(pcod)), nsim = 2L))
# #
# # pos <- subset(dogfish, catch_weight > 0)
# #
# # mm <- glm(catch_weight ~ 1, family = Gamma("log"), data = pos, offset = log(pos$area_swept))
# # mm_s <- sdmTMB(catch_weight ~ 1, family = Gamma("log"), data = pos, offset = log(pos$area_swept), spatial = "off")
# # pp1 <- predict(mm)
# # pp1_s <- predict(mm_s)
# #
# #
# #
# # head(p8)
# # head(p9)
