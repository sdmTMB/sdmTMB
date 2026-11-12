test_that("Standard errors on overall predictions from delta models work", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 12)
  fit <- sdmTMB(
    density ~ depth_scaled + I(depth_scaled^2),
    data = pcod,
    mesh = mesh,
    family = delta_gamma(),
  )
  fit

  nd <- data.frame(depth_scaled = seq(min(pcod$depth_scaled), max(pcod$depth_scaled), length.out = 50))

  # link
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "link")
  head(p)

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "link")
  head(p)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "link")
  head(p)


  # response
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "response")
  head(p)


  # re_form = NULL, se_fit = FALSE
  nd$X <- mean(pcod$X)
  nd$Y <- mean(pcod$Y)
  p <- predict(fit, model = NA, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 1, newdata = nd, type = "response")
  head(p)

  p <- predict(fit, model = 2, newdata = nd, type = "response")
  head(p)


  # link with se_fit = TRUE, re_form = NA
  p <- predict(fit, model = NA, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  expect_true("est" %in% names(p))
  expect_true("est_se" %in% names(p))
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(as.numeric(p$est), as.numeric(log(plogis(p$est1) * exp(p$est2))))

  p <- predict(fit, model = 1, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(p$est, p$est1)

  p <- predict(fit, model = 2, re_form = NA, newdata = nd, type = "link", se_fit = TRUE)
  head(p)
  plot(p$depth_scaled, p$est)
  lines(p$depth_scaled, p$est - 2 * p$est_se)
  lines(p$depth_scaled, p$est + 2 * p$est_se)
  expect_equal(p$est, p$est2)

  visreg_delta(fit, xvar = "depth_scaled", nn = 10, model = 1)
  visreg_delta(fit, xvar = "depth_scaled", nn = 10, model = 2)

  # standard errors are only available on the link scale:
  expect_error(
    predict(fit, model = NA, re_form = NA, newdata = nd, type = "response", se_fit = TRUE),
    regexp = "link scale"
  )

})

test_that("Draws honour re_form = NA", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  nd <- pcod_2011[1:20, ]
  check_draws <- function(fit, ...) {
    lp <- fit$tmb_obj$env$last.par.best
    samples <- cbind(lp, lp)
    for (re_form in list(NULL, NA)) {
      for (type in c("link", "response")) {
        p <- predict(fit, newdata = nd, re_form = re_form, type = type, ...)
        d <- predict(fit, newdata = nd, re_form = re_form, type = type,
          mcmc_samples = samples, ...)
        expect_equal(unname(d[, 1]), p$est, tolerance = 1e-8)
        expect_equal(d[, 1], d[, 2])
      }
    }
  }

  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
    family = tweedie(), time = "year")
  check_draws(fit)

  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
    family = delta_gamma(), time = "year")
  check_draws(fit)
  check_draws(fit, model = 1)
  check_draws(fit, model = 2)
})
