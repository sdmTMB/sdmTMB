test_that("coef and vcov and confint work", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  x <- coef(fit)
  expect_equal(round(unname(x), 3), c(5.347, -0.010))
  expect_equal(names(x), c("(Intercept)", "depth"))
  expect_equal(round(sigma(fit), 3), 16.887)

  x <- vcov(fit)
  expect_equal(nrow(x), 2L)
  expect_equal(ncol(x), 2L)
  expect_equal(colnames(x)[1], "(Intercept)")
  expect_equal(rownames(x)[1], "(Intercept)")

  x <- vcov(fit, complete = TRUE)
  expect_equal(nrow(x), 4L)
  expect_equal(ncol(x), 4L)

  x <- confint(fit)
  expect_equal(nrow(x), 2L)
  expect_equal(ncol(x), 3L)
  expect_true(grepl("2\\.5", colnames(x))[1])
  expect_true(grepl("97\\.5", colnames(x))[2])
  expect_true(grepl("Estimate", colnames(x))[3])
})

test_that("coef works with delta models and informs as needed", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = delta_gamma()
  )
  expect_message(x <- coef(fit), regexp = "model")
  expect_message(x <- coef(fit, model = 1), regexp = "model")
  expect_message(x <- coef(fit, model = 2), regexp = "model")

  expect_equal(round(sigma(fit), 3), 1.261)
})

test_that("set_delta_model is respected by model methods", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = delta_gamma()
  )
  fit1 <- set_delta_model(fit, model = 1)
  fit2 <- set_delta_model(fit, model = 2)

  expect_equal(suppressMessages(coef(fit1)), suppressMessages(coef(fit, model = 1)))
  expect_equal(suppressMessages(coef(fit2)), suppressMessages(coef(fit, model = 2)))
  expect_equal(vcov(fit1), vcov(fit, model = 1))
  expect_equal(vcov(fit2), vcov(fit, model = 2))
  expect_equal(fixef(fit1), fixef(fit, model = 1))
  expect_equal(fixef(fit2), fixef(fit, model = 2))

  fit$formula[[2]] <- density ~ year
  fit$terms[[2]] <- stats::terms(fit$formula[[2]])
  fit1 <- set_delta_model(fit, model = 1)
  fit2 <- set_delta_model(fit, model = 2)
  expect_equal(formula(fit1), fit$formula[[1]])
  expect_equal(formula(fit2), fit$formula[[2]])
  expect_equal(terms(fit1), fit$terms[[1]])
  expect_equal(terms(fit2), fit$terms[[2]])
  expect_equal(model.matrix(fit1), model.matrix(fit, model = 1))
  expect_equal(model.matrix(fit2), model.matrix(fit, model = 2))

  newdata <- fit$data[seq_len(5), , drop = FALSE]
  newdata$depth <- seq(min(fit$data$depth), max(fit$data$depth), length.out = nrow(newdata))
  expect_equal(nrow(model.matrix(fit, data = newdata)), nrow(newdata))
  expect_equal(
    model.matrix(fit, data = newdata),
    stats::model.matrix(
      fit$terms[[1]], data = newdata,
      contrasts.arg = fit$contrasts[[1]]
    )
  )

  expect_equal(vcov(fit1, model = 2), vcov(fit, model = 2))
})

test_that("various methods work", {
  skip_on_cran()
  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  f <- fitted(fit)
  expect_equal(nrow(pcod_2011), length(f))

  fit <- sdmTMB(
    density ~ depth,
    data = pcod_2011, spatial = "off",
    family = delta_gamma()
  )
  f <- fitted(fit)
  expect_equal(nrow(pcod_2011), length(f))

  a <- AIC(fit)
  expect_equal(round(a, 3), 6062.726)

  f <- fixef(fit)
  expect_length(f, 2L)

  f <- family(fit)
  expect_identical(f$family, c("binomial", "Gamma"))

  x <- terms(fit)
  expect_identical(as.character(x), c("~", "density", "depth"))

  pcod_2011$fyear <- as.factor(pcod_2011$year)
  fit <- sdmTMB(
    density ~ (1 | fyear),
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  x <- ranef(fit)
  expect_equal(nrow(x[[1]]$fyear), 4L)
})

test_that("more sigma methods work", {
  skip_on_cran()

  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, spatial = "off",
    family = tweedie(link = "log")
  )
  expect_equal(round(sigma(fit), 3), 17.592)

  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, spatial = "off",
    family = delta_lognormal()
  )
  expect_equal(round(sigma(fit), 3), 1.414)

  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, spatial = "off",
    family = poisson(link = "log")
  )
  expect_equal(sigma(fit), 1)
})
