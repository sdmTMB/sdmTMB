test_that("fits don't capture the calling frame through formula environments", {
  skip_on_cran()
  d <- pcod_2011
  d$depth_s <- as.numeric(scale(d$depth))
  fit_in_fn <- function() {
    big <- numeric(2e6) # should not be serialized with the fit
    deg <- 2 # referenced by the formula; must still resolve
    sdmTMB(density ~ s(depth_s), data = d,
      dispformula = ~ poly(depth_s, deg), spatial = "off", family = tweedie())
  }
  fit <- fit_in_fn()
  e <- environment(fit$formula[[1]])
  expect_identical(ls(e), "deg")
  expect_identical(environment(fit$terms[[1]]), e)
  expect_identical(environment(fit$smoothers$formula_no_sm), e)
  expect_identical(environment(fit$dispformula), e)
  expect_lt(as.numeric(object.size(serialize(fit, NULL))), 1e7)

  f <- tempfile(fileext = ".rds")
  saveRDS(fit, f)
  p <- predict(fit, newdata = d)
  fit2 <- readRDS(f)
  expect_equal(predict(fit2, newdata = d)$est, p$est)
})

test_that("symbols from an enclosing function frame are kept", {
  skip_on_cran()
  d <- pcod_2011
  d$depth_s <- as.numeric(scale(d$depth))
  outer_fn <- function() {
    big <- numeric(2e6) # should not be serialized with the fit
    deg <- 2 # defined one frame above the formula's environment
    inner_fn <- function() {
      sdmTMB(density ~ poly(depth_s, deg), data = d,
        spatial = "off", family = tweedie())
    }
    inner_fn()
  }
  fit <- outer_fn()
  expect_identical(ls(environment(fit$formula[[1]])), "deg")
  expect_lt(as.numeric(object.size(serialize(fit, NULL))), 1e7)

  f <- tempfile(fileext = ".rds")
  saveRDS(fit, f)
  p <- predict(fit, newdata = d)
  fit2 <- readRDS(f)
  expect_equal(predict(fit2, newdata = d)$est, p$est)
})

test_that("formulas from the global environment are untouched", {
  skip_on_cran()
  fit <- eval(quote(sdmTMB(density ~ depth_scaled, data = pcod_2011,
    spatial = "off", family = tweedie())), globalenv())
  expect_identical(environment(fit$formula[[1]]), globalenv())
})

test_that("emmeans works on a fit made inside a function", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  fit_in_fn <- function() {
    dat <- pcod_2011
    dat$fyear <- as.factor(dat$year)
    sdmTMB(density ~ fyear, data = dat, spatial = "off", family = tweedie())
  }
  fit <- fit_in_fn()
  e <- as.data.frame(emmeans::emmeans(fit, ~fyear))
  expect_equal(nrow(e), 4L)
})
