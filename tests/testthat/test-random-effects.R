# library(glmmTMB)
# library(sdmTMB)
# library(lme4)
# test_that("RE group factor levels are properly checked.", {
#   expect_error(check_valid_factor_levels(c(1, 2, 3), "test"))
#   expect_error(check_valid_factor_levels(c("A", "B")))
#   x <- factor(c("a", "b", "c"))
#   expect_true(check_valid_factor_levels(x))
#   x <- x[-1]
#   expect_error(check_valid_factor_levels(x, "test"))
# })

test_that("predict() handles NA values in response column of newdata", {
  skip_on_cran()

  set.seed(42)
  n <- 120
  d <- data.frame(
    cpue = ifelse(stats::runif(n) < 0.4, 0, exp(stats::rnorm(n))),
    x = stats::rnorm(n),
    g = factor(sample(letters[1:6], n, replace = TRUE))
  )
  fit <- sdmTMB(
    cpue ~ x + (1 | g),
    data = d,
    spatial = "off",
    family = delta_lognormal(),
    silent = TRUE
  )

  nd <- d
  nd$cpue[1] <- NA_real_
  p_with_response <- predict(fit, newdata = nd)
  p_without_response <- predict(fit, newdata = nd[, setdiff(names(nd), "cpue"), drop = FALSE])

  expect_equal(p_with_response$est1, p_without_response$est1, tolerance = 1e-10)
  expect_equal(p_with_response$est2, p_without_response$est2, tolerance = 1e-10)
})

test_that("Model with random intercepts fits appropriately.", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("lme4")
  # Classic lme4 random effects model
  data("sleepstudy", package="lme4")
  lmer_fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
  glmmTMB_fit <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
  # Same model in sdmTMB
  sdmTMB_fit <- sdmTMB(Reaction ~ Days + (Days | Subject), sleepstudy, spatial="off")

  # with smoothers! (was broken)
  sdmTMB_fit_smooth <- sdmTMB(Reaction ~ s(Days) + (Days | Subject), sleepstudy, spatial="off")

  # demo MVN draws:
  if (FALSE) {
    object <- sdmTMB_fit
    n_sims <- 200
    tmb_sd <- object$sd_report
    set.seed(1)
    samps <- sdmTMB:::rmvnorm_prec(object$tmb_obj$env$last.par.best, tmb_sd, n_sims)
    pars <- c(tmb_sd$par.fixed, tmb_sd$par.random)
    pn <- names(pars)
    samps[1:20,1:3]
    b_j <- which(pn == "b_j")[2]
    re_b_pars <- which(pn == "re_b_pars") # intercepts followed by random slopes here...
    re <- samps[re_b_pars, ]
    intercepts <- re[1:18,]
    slopes <- re[19:36,]
    b <- samps[b_j,,drop=FALSE]
    out <- matrix(nrow = nrow(slopes), ncol = ncol(slopes))
    for (i in 1:ncol(slopes)) {
      out[,i] <- b[,i] + slopes[,i]
    }
    med <- apply(out, 1, quantile, probs = 0.50)
    lwr <- apply(out, 1, quantile, probs = 0.05)
    upr <- apply(out, 1, quantile, probs = 0.95)
    plot(med, 1:18, xlim = range(c(lwr, upr)))
    segments(lwr, 1:18, upr, 1:18)
  }

  ps1 <- predict(sdmTMB_fit_smooth)
  ps2 <- predict(sdmTMB_fit_smooth, newdata = sleepstudy)
  expect_equal(ps1$est, ps2$est, tolerance = 1e-4)

  p0 <- predict(sdmTMB_fit, newdata = NULL)
  p1 <- predict(lmer_fit, newdata = sleepstudy)
  p2 <- predict(sdmTMB_fit, newdata = sleepstudy)
  expect_equal(as.numeric(p1), p2$est, tolerance = 1e-4)
  expect_equal(p0$est, p2$est, tolerance = 1e-4)

  # fixed effect only prediction:
  p0 <- predict(glmmTMB_fit, re.form = NA)
  p1 <- predict(sdmTMB_fit, re_form_iid = NA)
  expect_equal(p0, p1$est, tolerance = 1e-4)

  # missing factor level:
  ndtest <- sleepstudy
  ndtest <- ndtest[ndtest$Subject != "308", ]
  p1 <- predict(lmer_fit, newdata = ndtest)
  p2 <- predict(sdmTMB_fit, newdata = ndtest)
  expect_equal(as.numeric(p1), p2$est, tolerance = 1e-4)

  # levels themselves missing:
  ndtest <- sleepstudy
  ndtest <- ndtest[ndtest$Subject == "308", ]
  ndtest$Subject <- factor(as.character(ndtest$Subject))
  p1 <- predict(lmer_fit, newdata = ndtest)
  p2 <- predict(sdmTMB_fit, newdata = ndtest)
  expect_equal(as.numeric(p1), p2$est, tolerance = 1e-4)

  # Check fixed effects are identical
  expect_equal(lme4::fixef(lmer_fit)[1], coef(sdmTMB_fit)[1])
  expect_equal(lme4::fixef(lmer_fit)[2], coef(sdmTMB_fit)[2])
  # Check likelihood / AIC identical
  expect_equal(AIC(lmer_fit), AIC(sdmTMB_fit))
  # Check variances and covariance of REs is equal
  REs <- sdmTMB_fit$sd_report$value[grep("cov_pars", names(sdmTMB_fit$sd_report$value))]
  expect_equal(as.numeric(attr(summary(lmer_fit)$varcor[[1]], "stddev")),
               as.numeric(exp(REs[c(1,3)])), tolerance = 1.0e-4)
  expect_equal(as.numeric(attr(summary(lmer_fit)$varcor[[1]], "correlation")[1,2]),
               as.numeric(REs[2]), tolerance = 1.0e-2)

  # Check that ranef() returns the same thing
  expect_equal(mean(diag(cor(ranef(sdmTMB_fit)[[1]]$Subject, ranef(lmer_fit)$Subject))), 1)

  # verify model with no random int works
  sdmTMB_fit <- sdmTMB(Reaction ~ Days + (-1 + Days | Subject), sleepstudy, spatial="off")
  lmer_fit <- lme4::lmer(Reaction ~ Days + (-1 + Days | Subject), sleepstudy, REML = FALSE)
  expect_equal(fixef(lmer_fit)[1], coef(sdmTMB_fit)[1])
  expect_equal(fixef(lmer_fit)[2], coef(sdmTMB_fit)[2])
  REs <- sdmTMB_fit$sd_report$value[grep("cov_pars", names(sdmTMB_fit$sd_report$value))]
  expect_equal(as.numeric(attr(summary(lmer_fit)$varcor[[1]], "stddev"))[1],
               as.numeric(exp(REs[c(1)])), tolerance = 1.0e-4)

  # add a new level and verify multiple groups works
  data("sleepstudy", package="lme4")
  sleepstudy$age <- as.factor(rep(letters[1:5],36))
  set.seed(1)
  devs <- rnorm(5,0,1)
  sleepstudy$Reaction <- sleepstudy$Reaction + devs[rep(1:5,36)] + rnorm(nrow(sleepstudy),0,0.03)

  sdmTMB_fit <- sdmTMB(Reaction ~ Days + (1 + Days | Subject) + (1 | age), sleepstudy, spatial="off")
  # glmmtmb_fit <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject) + (1|age), sleepstudy, REML = FALSE)

  x <- capture.output(sdmTMB_fit)
  expect_true(any(grepl("Corr", x)))
  expect_true(any(grepl("Std.Dev.", x)))
  expect_true(any(grepl("Variance", x)))
  expect_true(any(grepl("565.71", x)))
  expect_true(any(grepl("0.08", x)))

  # library(sdmTMB)
  # data("sleepstudy", package="lme4")
  # sleepstudy$age <- as.factor(rep(letters[1:5],36))
  # set.seed(1)
  # devs <- rnorm(5,0,1)
  # sleepstudy$Reaction <- sleepstudy$Reaction + devs[rep(1:5,36)] + rnorm(nrow(sleepstudy),0,0.03)
  #
  # # random slopes and intercepts by subject
  # sdmTMB_fit <- sdmTMB(Reaction ~ Days + (Days | Subject) + (1 | age), sleepstudy, spatial="off")
  # sdmTMB_fit
  #
  # # random intercepts only by subject
  # sdmTMB_fit <- sdmTMB(Reaction ~ Days + (1 | Subject) + (1 | age), sleepstudy, spatial="off")
  # sdmTMB_fit
  #
  # # random slopes only by subject
  # sdmTMB_fit <- sdmTMB(Reaction ~ Days + (0 + Days | Subject) + (1 | age), sleepstudy, spatial="off")
  # sdmTMB_fit

  # sdmTMB_fit <- sdmTMB(Reaction ~ Days + (1 | Subject) + (1 | age), sleepstudy, spatial="off")
  # glmmtmb_fit <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | age), sleepstudy)
  # summary(glmmtmb_fit)
  # sdmTMB_fit


  sdmTMB_fit <- sdmTMB(Reaction ~ Days + (1 + Days | Subject) + (1 | age), sleepstudy, spatial="off")


  expect_equal(glmmTMB::fixef(glmmTMB_fit)$cond[1], coef(sdmTMB_fit)[1], tolerance = 1e-2)
  expect_equal(glmmTMB::fixef(glmmTMB_fit)$cond[2], coef(sdmTMB_fit)[2], tolerance = 1e-2)
  REs <- sdmTMB_fit$sd_report$value[grep("cov_pars", names(sdmTMB_fit$sd_report$value))]
  expect_equal(as.numeric(attr(summary(glmmTMB_fit)$varcor$cond$Subject, 'stddev')),
               as.numeric(exp(REs[c(1,3)])), tolerance = 1.0e-3)

  # Add in spatial field
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
  #tidy(m1, "fixed", conf.int = TRUE)
  .t1 <- tidy(m1, "ran_pars", conf.int = TRUE)

  # with RE:
  m <- sdmTMB(
     data = s, time = NULL,
     formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde
  )
  #tidy(m, "fixed", conf.int = TRUE)
  .t <- tidy(m, "ran_pars", conf.int = TRUE)
  #print(m)

  expect_gt(.t1$estimate[.t1$term == "phi"], .t$estimate[.t$term == "phi"])
  expect_gt(.t1$estimate[.t1$term == "sigma_O"], .t$estimate[.t$term == "sigma_O"])
  expect_equal(nrow(.t1), 3)
  expect_equal(nrow(.t), 5)

  b <- as.list(m$sd_report, "Estimate")
  .cor <- cor(c(RE_vals, RE_vals2), b$re_b_pars)
  expect_equal(round(c(.cor), 5), 0.8313)
  expect_equal(round(b$re_b_pars[seq_len(5)], 5),
    c(-0.28645, 0.68619, 0.10028, -0.31436, -0.61168),
    tolerance = 1e-5
  )
#
#   p <- predict(m)
#   p.nd <- predict(m, newdata = s)
#   # newdata is not the same as fitted data:
#   p.nd2 <- predict(m, newdata = s[1:3, , drop = FALSE])

#   expect_equal(p.nd2$est[1:3], p$est[1:3], tolerance = 1e-4)
#   expect_equal(p.nd2$est_non_rf[1:3], p$est_non_rf[1:3], tolerance = 1e-4)
#   expect_equal(p.nd2$est[1:3], p.nd$est[1:3], tolerance = 1e-9)
#   expect_equal(p$est, p.nd$est, tolerance = 1e-4)
#   expect_equal(p$est_rf, p.nd$est_rf, tolerance = 1e-4)
#   expect_equal(p$est_non_rf, p.nd$est_non_rf, tolerance = 1e-4)
#
#   # prediction with missing level in `newdata` works:
#   s_drop <- s[s$g != 1, , drop = FALSE]
#   p.nd <- predict(m, newdata = s_drop)
#   p <- p[s$g != 1, , drop = FALSE]
#   expect_equal(p$est, p.nd$est, tolerance = 1e-4)
#
#   # prediction without random intercepts included:
#   p.nd.null <- predict(m, newdata = s, re_form_iid = NULL)
#   p.nd.na <- predict(m, newdata = s, re_form_iid = NA)
#   p.nd.0 <- predict(m, newdata = s, re_form_iid = ~0)
#   expect_identical(p.nd.na, p.nd.0)
#   expect_false(identical(p.nd.null$est, p.nd.0$est))
#
  # random ints match glmmTMB exactly:
  m <- sdmTMB(
    data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h), mesh = spde, spatial = "off"
  )
  .t <- tidy(m, "ran_pars")
  m.glmmTMB <- glmmTMB::glmmTMB(
    data = s,
    formula = observed ~ 1 + (1 | g) + (1 | h)
  )
  # .v <- glmmTMB::VarCorr(m.glmmTMB)
  # expect_equal(.t$estimate[.t$term == "sigma_G"][1],
  #   sqrt(as.numeric(.v$cond$g)),
  #   tolerance = 1e-5
  # )
  # expect_equal(.t$estimate[.t$term == "sigma_G"][2],
  #   sqrt(as.numeric(.v$cond$h)),
  #   tolerance = 1e-5
  # )

  sdmTMB_re <- as.list(m$sd_report, "Estimate")
  glmmTMB_re <- glmmTMB::ranef(m.glmmTMB)$cond
  expect_equal(c(glmmTMB_re$g$`(Intercept)`, glmmTMB_re$h$`(Intercept)`),
    sdmTMB_re$re_b_pars[,1],
    tolerance = 1e-5
  )
#
  # predicting with new levels throws error for now:
  m <- sdmTMB(data = s, formula = observed ~ 1 + (1 | g), spatial = "off")
  nd <- data.frame(g = factor(c(1, 2, 3, 800)), observed=1)
  expect_error(predict(m, newdata = nd), regexp = "Extra")

  # predicting with missing factors works with the right re_form_iid
  m <- sdmTMB(data = s, formula = observed ~ 1 + (1 | g), spatial = "off")
  nd <- s[, !names(s) %in% "g", drop = FALSE]
  p1 <- predict(m, newdata = nd, re_form_iid = ~0)
  p2 <- predict(m, newdata = nd, re_form_iid = NA)
  # and this fails when not correct
  expect_error(predict(m, newdata = nd), regexp = "variable lengths differ")
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
  # ranpars <- tidy(m, "ran_pars", conf.int = TRUE)
  # s2 <- as.list(m2$sdr, "Estimate")
  # expect_equal(ranpars$estimate[-1], exp(s2$theta), tolerance = 0.001)
  # s2se <- as.list(m2$sdr, "Std. Error")
  # upr <- exp(s2$theta + 2 * s2se$theta)
  # lwr <- exp(s2$theta - 2 * s2se$theta)
  # expect_equal(ranpars$conf.low[-1], lwr, tolerance = 0.01)
  # expect_equal(ranpars$conf.high[-1], upr, tolerance = 0.01)

  ranint <- tidy(m, "ran_vals", conf.int = TRUE)

  expect_equal(ranef(m2)$cond$g[[1]],
    ranint$estimate[1:30],
    tolerance = 0.01
  )
#
#   # also check that ranef returns the same thing with same names
#   expect_equal(names(ranef(m2)$cond), names(ranef(m)$cond))
#
#   # and check that they return the same values
#   expect_equal(ranef(m2)$cond$g[[1]], ranef(m)$cond$g[[1]], tolerance = 1e-5)
})
#
#
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

test_that("Delta model works with random effects", {
  skip_on_cran()
  set.seed(1)

  data(pcod)
  pcod$year_f <- as.factor(pcod$year)

  # with single formula, the random effects should get carried through to all pieces
  m_yrf_re <- sdmTMB(
    data = pcod,
    formula = density ~ (1 | year_f),
    family = delta_gamma(),
    spatial = "off"
  )
  expect_equal(nrow(tidy(m_yrf_re, "ran_vals")), length(unique(pcod$year))*2)


  # test 2 different RE intercepts
 m_yrf_re_1 <- sdmTMB(
    data = pcod,
    formula = list(density ~ (1 | year_f), density ~ (1|year_f)),
    family = delta_gamma(),
    spatial = "off"
  )

  # 2 diff intercepts, same number of levels
  intcpts <- rnorm(9)
  pcod$vessel <- sample(1:9, size = nrow(pcod), replace=T)
  pcod$density[which(pcod$present==1)] <- exp(log(pcod$density[which(pcod$present==1)]) + intcpts[pcod$vessel[which(pcod$present==1)]])
  pcod$vessel <- as.factor(pcod$vessel)
  pcod$density[which(pcod$present==1)] <- exp(log(pcod$density[which(pcod$present==1)]) + intcpts[pcod$vessel[which(pcod$present==1)]])

  m_yrf_re2 <- sdmTMB(
    data = pcod,
    formula = list(density ~ (1 | year_f), density ~ (1|vessel)),
    family = delta_gamma(),
    spatial = "off"
  )
  glmm_pres <- glmmTMB::glmmTMB(
    data = pcod,
    formula = present ~ (1 | year),
    family = binomial()
  )
  log_vars <- m_yrf_re2$sd_report$value[grep("re_cov_pars", names(m_yrf_re2$sd_report$value))]
  expect_equal(as.numeric(attr(summary(glmm_pres)$varcor[[1]]$year, "stddev")), as.numeric(exp(log_vars[1])),
               tolerance = 1e-5)

  # 2 diff intercepts, different number of levels
  intcpts <- rnorm(10)
  pcod$vessel <- sample(1:10, size = nrow(pcod), replace=T)
  pcod$density[which(pcod$present==1)] <- exp(log(pcod$density[which(pcod$present==1)]) + intcpts[pcod$vessel[which(pcod$present==1)]])
  pcod$vessel <- as.factor(pcod$vessel)
  pcod$density[which(pcod$present==1)] <- exp(log(pcod$density[which(pcod$present==1)]) + intcpts[pcod$vessel[which(pcod$present==1)]])

  m_yrf_re3 <- sdmTMB(
    data = pcod,
    formula = list(density ~ (-1+depth | year_f), density ~ (1|vessel)),
    family = delta_gamma(),
    spatial = "off"
  )

  # test 2 different numbers of random ints
  m_yrf_re4 <- sdmTMB(
    data = pcod,
    formula = list(density ~ (1 | year_f) + (1|vessel), density ~ (1|vessel)),
    family = delta_gamma(),
    spatial = "off"
  )


  # test 2 different numbers of random ints, different number of levels
  # m_yrf_re5 <- sdmTMB(
  #   data = pcod,
  #   formula = list(density ~ (depth | year_f) + (1|vessel), density ~ (1|vessel)),
  #   family = delta_gamma(),
  #   spatial = "off"
  # )

  # # Test same model with characters
  # pcod$year_chr <- paste(pcod$year)
  # m_yrf_re6 <- sdmTMB(
  #   data = pcod,
  #   formula = list(density ~ (depth | year_chr) + (1|vessel), density ~ (1|vessel)),
  #   family = delta_gamma(),
  #   spatial = "off"
  # )
  #
  #
  # # Test same model with integers
  # m_yrf_re7 <- sdmTMB(
  #   data = pcod,
  #   formula = list(density ~ (depth | year) + (1|vessel), density ~ (1|vessel)),
  #   family = delta_gamma(),
  #   spatial = "off"
  # )
})

test_that("issue breakpt() version of formula doesn't break random effect prediction #423", {
  d <- pcod
  d$year_f <- as.factor(pcod$year)
  m <- sdmTMB(
    data = d,
    formula = density ~ 0 + breakpt(depth_scaled) + (1 | year_f),
    spatial = "off",
    priors = sdmTMBpriors(
      threshold_breakpt_slope = normal(0, 1),
      threshold_breakpt_cut = normal(0, 1)
    ),
    family = tweedie(link = "log")
  )
  nd <- data.frame(
    depth_scaled = seq(min(pcod$depth_scaled) + 0.5,
      max(pcod$depth_scaled) - 0.2,
      length.out = 100
    ),
    year = 2015L
  )
  # Works because I turn off random intercepts
  p <- predict(m,
    newdata = nd, se_fit = TRUE, re_form = NA,
    re_form_iid = NA, xy_cols = c("X", "Y")
  )
  # Throws breakpoint-related error but it's actually because I
  # don't have the random factor in newdata but trying to predict with it
  expect_error({
  p <- predict(m,
    newdata = nd, se_fit = TRUE, re_form = NA,
    xy_cols = c("X", "Y")
  )}, regexp = "year_f")
})

