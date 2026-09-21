test_that("RE group factor levels are properly checked.", {
  expect_error(check_valid_factor_levels(c(1, 2, 3), "test"))
  expect_error(check_valid_factor_levels(c("A", "B")))
  x <- factor(c("a", "b", "c"))
  expect_true(check_valid_factor_levels(x))
  x <- x[-1]
  expect_error(check_valid_factor_levels(x, "test"))
})

test_that("Model with random intercepts fits appropriately.", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  s <- sdmTMB_simulate(
    ~1,
    data = loc,
    mesh = spde,
    range = 1.4,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 1,
    B = 0
  )

  g <- rep(gl(30, 10), 999)
  set.seed(134)
  RE_vals <- rnorm(30, 0, 0.4)
  h <- rep(gl(40, 10), 999)
  set.seed(1283)
  RE_vals2 <- rnorm(40, 0, 0.2)
  s$g <- g[seq_len(nrow(s))]
  s$h <- h[seq_len(nrow(s))]
  s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

  # ignore RE:
  m1 <- sdmTMB(data = s, formula = observed ~ 1, mesh = spde)
  tidy(m1, "fixed", conf.int = TRUE)
  .t1 <- tidy(m1, "ran_pars", conf.int = TRUE)

  # with RE:
  m <- sdmTMB(
    data = s, time = NULL,
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde
  )
  tidy(m, "fixed", conf.int = TRUE)
  .t <- tidy(m, "ran_pars", conf.int = TRUE)
  print(m)

  expect_gt(.t1$estimate[.t1$term == "phi"], .t$estimate[.t$term == "phi"])
  expect_gt(.t1$estimate[.t1$term == "sigma_O"], .t$estimate[.t$term == "sigma_O"])
  expect_equal(nrow(.t1), 3)
  expect_equal(nrow(.t), 5)

  b <- as.list(m$sd_report, "Estimate")
  .cor <- cor(c(RE_vals, RE_vals2), b$re_b_pars[, 1])
  expect_equal(round(.cor, 5), 0.8313)
  expect_equal(round(b$re_b_pars[seq_len(5)], 5),
    c(-0.28645, 0.68619, 0.10028, -0.31436, -0.61168),
    tolerance = 1e-5
  )

  p <- predict(m)
  p.nd <- predict(m, newdata = s)
  # newdata is not the same as fitted data:
  p.nd2 <- predict(m, newdata = s[1:3, , drop = FALSE])

  expect_equal(p.nd2$est[1:3], p$est[1:3], tolerance = 1e-4)
  expect_equal(p.nd2$est_non_rf[1:3], p$est_non_rf[1:3], tolerance = 1e-4)
  expect_equal(p.nd2$est[1:3], p.nd$est[1:3], tolerance = 1e-9)
  expect_equal(p$est, p.nd$est, tolerance = 1e-4)
  expect_equal(p$est_rf, p.nd$est_rf, tolerance = 1e-4)
  expect_equal(p$est_non_rf, p.nd$est_non_rf, tolerance = 1e-4)

  # prediction with missing level in `newdata` works:
  s_drop <- s[s$g != 1, , drop = FALSE]
  p.nd <- predict(m, newdata = s_drop)
  p <- p[s$g != 1, , drop = FALSE]
  expect_equal(p$est, p.nd$est, tolerance = 1e-4)

  # prediction without random intercepts included:
  p.nd.null <- predict(m, newdata = s, re_form_iid = NULL)
  p.nd.na <- predict(m, newdata = s, re_form_iid = NA)
  p.nd.0 <- predict(m, newdata = s, re_form_iid = ~0)
  expect_identical(p.nd.na, p.nd.0)
  expect_false(identical(p.nd.null$est, p.nd.0$est))

  # random ints match glmmTMB exactly:
  m <- sdmTMB(
    data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde, spatial = "off",
    control = sdmTMBcontrol(nlminb_loops = 3, newton_loops = 5)
  )
  .t <- tidy(m, "ran_vcov")
  m.glmmTMB <- glmmTMB::glmmTMB(
    data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h)
  )
  .v <- glmmTMB::VarCorr(m.glmmTMB)
  expect_equal(as.numeric(.t[[1]])[1],
    sqrt(as.numeric(.v$cond$g)),
    tolerance = 1e-5
  )
  expect_equal(as.numeric(.t[[1]])[2],
    sqrt(as.numeric(.v$cond$h)),
    tolerance = 1e-5
  )

  a <- as.list(m.glmmTMB$sdr, "Std. Error")
  b <- as.list(m$sd_report, "Std. Error")
  expect_equal(b$re_b_pars[,1,drop=TRUE], a$b, tolerance = 1e-3)

  # with CIs on "vcov" works:
  .t <- tidy(m, effects="ran_vcov", conf.int = FALSE)
  expect_equal(length(.t$est), 2)
  expect_equal(names(.t$est), c("Model 1 Group g", "Model 1 Group h"))
  .t <- tidy(m, effects="ran_vcov", conf.int = TRUE)
  expect_equal(names(.t), c("est","lo","hi"))
  expect_equal(length(.t$lo), 2)
  expect_equal(names(.t$lo), c("Model 1 Group g", "Model 1 Group h"))
  expect_equal(length(.t$hi), 2)
  expect_equal(names(.t$hi), c("Model 1 Group g", "Model 1 Group h"))
  expect_equal(as.numeric(unlist(.t$hi)), c(0.5808372, 0.2717274), tolerance = 0.001)
  expect_equal(as.numeric(unlist(.t$lo)), c(0.3367532, 0.1414036), tolerance = 0.001)
  expect_equal(as.numeric(unlist(.t$est)), c(0.4422655, 0.1960184), tolerance = 0.001)


  sdmTMB_re <- as.list(m$sd_report, "Estimate")
  glmmTMB_re <- glmmTMB::ranef(m.glmmTMB)$cond
  expect_equal(c(glmmTMB_re$g$`(Intercept)`, glmmTMB_re$h$`(Intercept)`),
    sdmTMB_re$re_b_pars[, 1],
    tolerance = 1e-5
  )

  # predicting with new levels uses population-level predictions:
  m <- sdmTMB(data = s, formula = observed ~ 1 + (1 | g), spatial = "off")
  nd <- data.frame(g = factor(c(800)))
  expect_warning(
    p_warn <- predict(m, newdata = nd),
    regexp = "Found new levels"
  )
  p_silent <- predict(m, newdata = nd, allow_new_levels = TRUE)
  expect_equal(p_warn$est, p_silent$est)
  expect_equal(p_silent$est, as.numeric(coef(m)))
})

test_that("Random intercepts and cross validation play nicely", {
  skip_on_cran()
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")
  s <- sdmTMB_simulate(
    ~1,
    data = loc,
    mesh = spde,
    range = 1.4,
    phi = 0.1,
    sigma_O = 0,
    seed = 1,
    B = 0
  )
  g <- rep(gl(3, 10), 999)
  RE_vals <- rnorm(3, 0, 0.4)
  s$g <- g[seq_len(nrow(s))]
  s$observed <- s$observed + RE_vals[s$g]
  # one level in its own fold:
  fold_ids <- as.integer(s$g %in% c(1, 2)) + 1
  out <- sdmTMB_cv(
    observed ~ 1 + (1 | g),
    fold_ids = fold_ids, k_folds = 2L, spatial = "off", data = s, mesh = spde,
    parallel = FALSE
  )
  expect_equal(round(out$sum_loglik, 3), -51.36)
  # Because the function fits with all the data but sets the missing fold to
  # have likelihood weights of 0, the fitted model is aware of all levels and
  # the missing levels just get left at a value of 0 because the data never
  # inform the model, which is exactly what you want.
})

test_that("Tidy returns random intercepts appropriately.", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  x <- stats::runif(500, -1, 1)
  y <- stats::runif(500, -1, 1)
  loc <- data.frame(x = x, y = y)
  spde <- make_mesh(loc, c("x", "y"), n_knots = 50, type = "kmeans")

  s <- sdmTMB_simulate(
    ~1,
    data = loc,
    mesh = spde,
    range = 1.4,
    phi = 0.1,
    sigma_O = 0.2,
    seed = 1,
    B = 0
  )

  g <- rep(gl(30, 10), 999)
  set.seed(134)
  RE_vals <- rnorm(30, 0, 0.4)
  h <- rep(gl(40, 10), 999)
  set.seed(1283)
  RE_vals2 <- rnorm(40, 0, 0.2)
  s$g <- g[seq_len(nrow(s))]
  s$h <- h[seq_len(nrow(s))]
  s$observed <- s$observed + RE_vals[s$g] + RE_vals2[s$h]

  # with RE; check against glmmTMB
  m <- sdmTMB(
    data = s, time = NULL,
    formula = observed ~ 1 + (1 | g) + (1 | h),
    mesh = spde,
    spatial = "off"
  )
  m2 <- glmmTMB::glmmTMB(
    data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h)
  )

  ranpars <- tidy(m, "ran_vcov", conf.int = TRUE)
  s2 <- as.list(m2$sdr, "Estimate")
  expect_equal(as.numeric(ranpars$est), exp(s2$theta), tolerance = 0.001)
  s2se <- as.list(m2$sdr, "Std. Error")
  upr <- exp(s2$theta + 2 * s2se$theta)
  lwr <- exp(s2$theta - 2 * s2se$theta)
  expect_equal(as.numeric(ranpars$lo), lwr, tolerance = 0.05)
  expect_equal(as.numeric(ranpars$hi), upr, tolerance = 0.05)

  ranint <- tidy(m, "ran_vals", conf.int = TRUE)

  expect_equal(ranef(m2)$cond$g[[1]],
    ranint$estimate[1:30],
    tolerance = 0.01
  )

  # also check that ranef returns the same thing with same names
  expect_equal(names(ranef(m2)$cond), lapply(ranef(m), names)[[1]])

  # and check that they return the same values
  expect_equal(ranef(m2)$cond$g[[1]], ranef(m)[[1]]$g[,1], tolerance = 1e-5)
})


test_that("Random intercept classes in predict() are checked appropriately", {
  skip_on_cran()
  set.seed(1)

  pcod$year_f <- as.factor(pcod$year)
  pcod_yrf_as_num <- pcod_yrf_as_chr <- pcod
  pcod_yrf_as_num$year_f <- as.numeric(pcod$year)
  pcod_yrf_as_chr$year_f <- as.character(pcod$year)

  m_yrf_re <- sdmTMB(
    data = pcod,
    formula = density ~ poly(log(depth), 2) + (1 | year_f),
    family = tweedie(link = "log"),
    spatial = "off"
  )

  expect_error(
    p8 <- predict(m_yrf_re, newdata = pcod_yrf_as_num),
    regexp = "newdata"
  )

  expect_error(
    p9 <- predict(m_yrf_re, newdata = pcod_yrf_as_chr),
    regexp = "newdata"
  )

  expect_error(
    p10 <- predict(m_yrf_re, newdata = pcod_yrf_as_num, re_form = NA),
    regexp = "newdata"
  )

  p11 <- predict(m_yrf_re, newdata = pcod_yrf_as_num,  # This should work
                 re_form = NA, re_form_iid = NA)

  expect_s3_class(p11, "tbl_df")
})

test_that("Random intercepts are incorporated into est_non_rf1 and est_non_rf2 correctly #342", {
  ## test based on example by @tom-peatman #342
  skip_on_cran()
  skip_on_ci()
  pcod_2011$year <- factor(pcod_2011$year)
  ## fit model with random intercepts for year in both model components
  fit <- sdmTMB(
    density ~ (1 | year),
    spatial = "off",
    data = pcod_2011, mesh = pcod_mesh_2011,
    family = delta_gamma(link1 = "logit", link2 = "log")
  )
  ## fitted random intercepts
  t1 <- tidy(fit, effects = "ran_vals")
  expect_equal(t1$estimate[t1$model==1], c(-0.24066, 0.33083, 0.20459, -0.29288), tolerance = 0.001)
  expect_equal(t1$estimate[t1$model==2], c(0.17724, -0.14142, 0.11862, -0.16844), tolerance = 0.001)

  ## data frame to get predictions for each year (at a specific location and depth)
  new_dat <- expand.grid(
    X = mean(pcod_2011$X), Y = mean(pcod_2011$Y),
    year = unique(pcod_2011$year)
  )

  p <- predict(fit, newdata = new_dat)
  expect_equal(p$est_non_rf1, c(-0.40011, 0.17137, 0.04514, -0.45234), tolerance = 0.001)
  expect_equal(p$est_non_rf2, c(4.6455, 4.32684, 4.58688, 4.29982), tolerance = 0.001)
})
