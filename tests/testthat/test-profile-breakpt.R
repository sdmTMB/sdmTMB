sim_breakpt <- function(n = 300, cut = 0.3, sd = 0.4, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(x = x, y = 1 + pmin(x, cut) + rnorm(n, sd = sd))
}

test_that("profile_breakpt() recovers the cutpoint and fits the final model", {
  skip_on_cran()
  d <- sim_breakpt()
  p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off")
  e <- p$estimate
  expect_s3_class(p, "breakpt_profile")
  expect_true(e$lwr < 0.3 && e$upr > 0.3)
  expect_true(e$lwr < e$cut && e$cut < e$upr)
  expect_false(e$open_interval)
  expect_false(e$few_below)
  expect_true(e$refine_complete)
  expect_equal(e$n_failed, 0)
  expect_equal(min(p$profile$dnll, na.rm = TRUE), 0)
  expect_false(is.unsorted(p$profile$cut))

  # the final fit holds the cutpoint at the best value and matches the profile
  b <- p$fit$model$par[names(p$fit$model$par) == "b_threshold"]
  expect_length(b, 1L) # only the slope is estimated
  expect_equal(p$fit$tmb_obj$env$parList()$b_threshold[2, 1], e$cut)
  expect_equal(e$final_nll, min(p$profile$nll, na.rm = TRUE), tolerance = 1e-4)
  expect_equal(e$aic, 2 * e$nll + 2 * e$k)
  expect_equal(e$k, length(p$fit$model$par) + 1L)
  # AIC() and BIC() count the cutpoint; AIC() on the final fit doesn't
  expect_equal(AIC(p), e$aic)
  expect_equal(AIC(p), AIC(p$fit) + 2, tolerance = 1e-6)
  expect_equal(BIC(p), -2 * logLik(p)[[1]] + log(nrow(d)) * e$k)
  # the final fit's call names the caller's objects, so print() and update() work
  expect_identical(p$fit$call$data, quote(d))
  expect_null(p$fit$call$final)
  expect_output(print(p$fit), "Data: d")
  # the stored call keeps the caller's control, not the starts and map that fix the cutpoint
  expect_null(p$fit$call$control)
  # the estimate reports the final fit's slope
  expect_equal(e$slope, tidy(p$fit)$estimate[[2L]], tolerance = 1e-8)
  # the profile, the estimate, and the final fit agree on the minimum
  expect_equal(e$nll, min(p$profile$nll, na.rm = TRUE))
  expect_lte(e$nll, e$final_nll)

  # the profile minimum is no worse than fitting the cutpoint directly from the truth
  m <- suppressMessages(suppressWarnings(sdmTMB(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(start = list(b_threshold = matrix(c(1, 0.3), 2L, 1L))))))
  expect_lt(e$nll, m$model$objective + 0.01) # within the grid's resolution
  expect_output(print(p), "Cutpoint")
})

test_that("profile_breakpt() warns when the cutpoint is poorly identified", {
  skip_on_cran()
  set.seed(2)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + d$x + rnorm(200, sd = 0.5) # linear: no bend within the data
  w <- character()
  p <- withCallingHandlers(
    profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", final = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(p$estimate$open_interval)
  expect_true(any(grepl("end of the grid", w)))
  expect_null(p$fit)
  expect_output(print(p), "reaches the end of the grid")

  # an open interval: the shaded band runs to the end of the grid
  skip_if_not_installed("ggplot2")
  g <- plot(p)
  expect_s3_class(g, "ggplot")
  band <- ggplot2::layer_data(g, 1)
  expect_equal(band$xmax, max(p$profile$cut))
  expect_equal(band$xmin, p$estimate$lwr)
})

test_that("profile_breakpt() counts parameters like logLik() with `control$profile`", {
  skip_on_cran()
  d <- sim_breakpt(n = 150)
  # `profile = TRUE` moves b_j into the inner problem, but logLik() still counts it
  p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(profile = TRUE))
  expect_equal(p$estimate$k, attr(logLik(p$fit), "df") + 1L)
  expect_equal(AIC(p), AIC(p$fit) + 2, tolerance = 1e-6)
  expect_identical(p$fit$call$control, quote(sdmTMBcontrol(profile = TRUE)))
  p0 <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(profile = TRUE), final = FALSE)
  expect_equal(p0$estimate$k, p$estimate$k)
})

test_that("profile_breakpt() applies `control` bounds to the profile fits", {
  skip_on_cran()
  d <- sim_breakpt(n = 150) # true slope 1
  # a bound on b_threshold bounds only the slope; the cutpoint is fixed in each fit
  p <- suppressWarnings(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(upper = list(b_threshold = 0.5)), final = FALSE))
  expect_true(all(p$profile$slope <= 0.5 + 1e-8, na.rm = TRUE))
  expect_equal(p$estimate$slope, 0.5, tolerance = 1e-4)
})

test_that("profile_breakpt() flags few observations above the cutpoint", {
  skip_on_cran()
  set.seed(3)
  x <- rnorm(300)
  cut <- quantile(x, 0.975, names = FALSE)
  d <- data.frame(x = x, y = 1 + 2 * pmin(x, cut) + rnorm(300, sd = 0.2))
  w <- character()
  p <- withCallingHandlers(
    profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", final = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(p$estimate$few_above)
  expect_false(p$estimate$few_below)
  expect_true(any(grepl("fall above", w)))
  expect_output(print(p), "few observations fall above")
})

test_that("profile_breakpt() checks its inputs", {
  skip_on_cran()
  d <- sim_breakpt(n = 100)
  expect_error(profile_breakpt(y ~ x, data = d, spatial = "off"), "breakpt")
  expect_error(profile_breakpt(y ~ logistic(x), data = d, spatial = "off"), "breakpt")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    control = sdmTMBcontrol(map = list(b_threshold = factor(c(1, NA))))), "map")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", tol = 0), "positive")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", max_fits = Inf), "finite")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", cuts = c(0, 1)), "3")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off", level = 95), "level")
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    priors = sdmTMBpriors(threshold_breakpt_cut = normal(0, 1))), "prior")
  d$y <- ifelse(d$y > 1, d$y, 0)
  expect_error(profile_breakpt(y ~ breakpt(x), data = d, spatial = "off",
    family = delta_gamma()), "Delta")
})

test_that("sdmTMB() points to profile_breakpt() when estimating a cutpoint directly", {
  skip_on_cran()
  d <- sim_breakpt(n = 100)
  # every fit, not once per session
  for (i in 1:2) {
    expect_message(suppressWarnings(sdmTMB(y ~ breakpt(x), data = d, spatial = "off")),
      "profile_breakpt")
  }
  # not when the cutpoint is fixed, as profile_breakpt() does
  msgs <- character()
  withCallingHandlers(
    sdmTMB(y ~ breakpt(x), data = d, spatial = "off",
      control = sdmTMBcontrol(start = list(b_threshold = matrix(c(1, 0.3), 2L, 1L)),
        map = list(b_threshold = factor(c(1, NA))))),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("profile_breakpt", msgs)))

  # delta models still get the warning, but not the advice to use profile_breakpt()
  d$y <- ifelse(d$y > 1, d$y, 0)
  msgs <- character()
  m <- withCallingHandlers(
    suppressWarnings(sdmTMB(y ~ breakpt(x), data = d, spatial = "off", family = delta_gamma())),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  msgs <- paste(msgs, collapse = " ")
  expect_match(msgs, "isn't smooth in the cutpoint")
  expect_match(msgs, "doesn't support delta models")
  expect_no_match(msgs, "Consider `profile_breakpt\\(\\)`")
  expect_false(sanity(m, silent = TRUE)$breakpt_ok)
  expect_message(sanity(m), "doesn't support delta models")
})

test_that("sdmTMB_cv() points to profile_breakpt() once, not for every fold", {
  skip_on_cran()
  d <- sim_breakpt(n = 100)
  msgs <- character()
  withCallingHandlers(
    suppressWarnings(sdmTMB_cv(y ~ breakpt(x), data = d, spatial = "off", k_folds = 3)),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  expect_equal(sum(grepl("profile_breakpt", msgs)), 1L)
})

test_that("sanity() flags a directly estimated cutpoint but not a profiled one", {
  skip_on_cran()
  d <- sim_breakpt()
  m <- suppressMessages(suppressWarnings(sdmTMB(y ~ breakpt(x), data = d, spatial = "off")))
  expect_message(s <- sanity(m), "profile_breakpt")
  expect_false(s$breakpt_ok)
  expect_false(s$all_ok)

  p <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off")
  s <- sanity(p$fit, silent = TRUE)
  expect_true(s$breakpt_ok)

  m <- sdmTMB(y ~ x, data = d, spatial = "off")
  expect_true(sanity(m, silent = TRUE)$breakpt_ok)
})

test_that("profile_breakpt() evaluates its `...` arguments once", {
  skip_on_cran()
  d <- sim_breakpt(n = 150)
  n_eval <- 0
  get_data <- function() {
    n_eval <<- n_eval + 1
    d
  }
  p <- profile_breakpt(y ~ breakpt(x), data = get_data(), spatial = "off")
  expect_s3_class(p$fit, "sdmTMB")
  expect_equal(n_eval, 1)
})

test_that("profile_breakpt() works with breakpt(positive = TRUE)", {
  skip_on_cran()
  d <- sim_breakpt()
  p0 <- profile_breakpt(y ~ breakpt(x), data = d, spatial = "off")
  p <- profile_breakpt(y ~ breakpt(x, positive = TRUE), data = d, spatial = "off")
  e <- p$estimate
  # the slope is reported as itself, not log(slope)
  expect_gt(e$slope, 0)
  expect_true(all(p$profile$slope > 0, na.rm = TRUE))
  expect_equal(e$slope, tidy(p$fit)$estimate[[2L]], tolerance = 1e-8)
  expect_equal(exp(p$fit$model$par[["b_threshold"]]),
    tidy(p$fit)$estimate[[2L]], tolerance = 1e-8)
  # with a clearly positive slope, it matches the unconstrained profile
  expect_equal(e$cut, p0$estimate$cut, tolerance = 1e-3)
  expect_equal(e$slope, p0$estimate$slope, tolerance = 1e-3)
  expect_equal(e$nll, p0$estimate$nll, tolerance = 1e-4)
  expect_equal(c(e$lwr, e$upr), c(p0$estimate$lwr, p0$estimate$upr), tolerance = 1e-2)
  expect_output(print(p), sprintf("Slope %.3f", e$slope))
})

test_that("profile_breakpt(positive = TRUE) doesn't carry a zero slope across the profile", {
  skip_on_cran()
  # rises steeply below the cutpoint and declines above it, so with the cutpoint at the
  # median the best positive slope is 0; fits warm-started from there used to stay at 0
  set.seed(1)
  x <- rnorm(400)
  d <- data.frame(x = x, y = 3 * pmin(x, -1.3) - pmax(x, -1.3) + rnorm(400, sd = 0.3))
  # few observations fall below the best cutpoint, which the profile warns about
  p <- suppressWarnings(profile_breakpt(y ~ breakpt(x, positive = TRUE), data = d,
    spatial = "off", final = FALSE))
  # a fit with the cutpoint fixed at the truth, started from a moderate slope
  m <- suppressMessages(sdmTMB(y ~ breakpt(x, positive = TRUE), data = d, spatial = "off",
    control = sdmTMBcontrol(start = list(b_threshold = matrix(c(0, -1.3), 2L, 1L)),
      map = list(b_threshold = factor(c(1, NA))))))
  expect_lt(p$estimate$nll, m$model$objective + 1e-4)
  expect_gt(p$estimate$slope, 1)
  # the decline above the cutpoint isn't in the model, so the best cutpoint is below the truth
  expect_lt(p$estimate$cut, -1)

  # the same profile with x in other units
  d2 <- transform(d, x = 100 * x)
  p2 <- suppressWarnings(profile_breakpt(y ~ breakpt(x, positive = TRUE), data = d2,
    spatial = "off", final = FALSE))
  expect_equal(p2$estimate$nll, p$estimate$nll, tolerance = 1e-4)
  expect_equal(p2$estimate$cut, 100 * p$estimate$cut, tolerance = 1e-3)
  expect_equal(p2$estimate$slope, p$estimate$slope / 100, tolerance = 1e-3)
})
