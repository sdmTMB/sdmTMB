test_that("predict() routes omitted newdata and resolves offsets by row source", {
  skip_on_cran()

  d <- pcod_2011
  d$log_effort <- log(runif(nrow(d), 0.5, 2))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 30)
  fit <- sdmTMB(
    density ~ depth_scaled,
    data = d, mesh = mesh, offset = "log_effort",
    family = tweedie(), spatial = "on"
  )
  fitted_offset <- fit$tmb_data$offset_i

  # fast path: fitted report, fitted offset
  p <- predict(fit)
  r <- fit$tmb_obj$report(fit$tmb_obj$env$last.par.best)
  expect_equal(p$est, r$eta_i[, 1])
  expect_equal(predict(fit, return_tmb_report = TRUE), r)

  # projection on the fitted data carries over the fitted offset
  p_resp <- predict(fit, type = "response")
  expect_equal(p_resp$est, exp(p$est), tolerance = 1e-6)
  r_proj <- predict(fit, re_form_iid = NA, return_tmb_report = TRUE)
  expect_true("proj_eta" %in% names(r_proj))
  expect_equal(r_proj$proj_eta[, 1], p$est, tolerance = 1e-6)

  # return_tmb_data prepares fitted rows even when newdata is omitted
  td <- predict(fit, return_tmb_data = TRUE)
  expect_type(td, "list")
  expect_equal(td$proj_offset_i, fitted_offset)
  expect_equal(nrow(td$proj_X_ij[[1]]), nrow(d))
  expect_equal(predict(fit, offset = rep(1, nrow(d)),
    return_tmb_data = TRUE)$proj_offset_i, rep(1, nrow(d)))

  # supplied newdata defaults to zero offset (with a message) unless overridden
  expect_message(td_nd <- predict(fit, newdata = d, return_tmb_data = TRUE),
    "offset")
  expect_equal(td_nd$proj_offset_i, rep(0, nrow(d)))
  td_off <- predict(fit, newdata = d, offset = d$log_effort, return_tmb_data = TRUE)
  expect_equal(td_off$proj_offset_i, d$log_effort)
  p_nd <- predict(fit, newdata = d, offset = d$log_effort)
  expect_equal(p_nd$est, p$est, tolerance = 1e-6)

  # data return takes precedence over report and simulation output
  td_sim <- predict(fit, nsim = 2, return_tmb_data = TRUE,
    return_tmb_report = TRUE)
  expect_true("proj_X_ij" %in% names(td_sim))
  r_sim <- predict(fit, nsim = 2, return_tmb_report = TRUE)
  expect_length(r_sim, 2L)
  expect_true("proj_eta" %in% names(r_sim[[1]]))
})

test_that("predict() data preparation does not need a fitted objective", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 30)
  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
    family = tweedie(), spatial = "on")
  nd <- pcod_2011[1:10, ]
  expected <- predict(fit, newdata = nd, return_tmb_data = TRUE)

  # any use of the objective (including reinitialization) now errors:
  no_obj <- fit
  no_obj$tmb_obj <- list(env = list(
    ADFun = list(ptr = new("externalptr")),
    beSilent = function() stop("fitted objective used")
  ))
  no_obj$model <- NULL
  expect_equal(predict(no_obj, newdata = nd, return_tmb_data = TRUE), expected)
  expect_equal(predict(no_obj, return_tmb_data = TRUE)$proj_offset_i,
    fit$tmb_data$offset_i)

  # fits saved with `parlist` predict without the fitted objective:
  expect_equal(predict(no_obj, newdata = nd)$est, predict(fit, newdata = nd)$est)
  expect_error(predict(no_obj), "fitted objective used")
})

test_that("predict() standard errors are selected from sdreport() untransformed", {
  skip_on_cran()

  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 30)
  fit <- sdmTMB(density ~ depth_scaled, data = pcod_2011, mesh = mesh,
    family = delta_gamma(), spatial = "off")
  nd <- pcod_2011[1:10, ]
  x <- predict_sdmTMB(fit, newdata = nd, re_form = NA, se_fit = TRUE,
    return_tmb_object = TRUE)
  sr <- sdreport_sdmTMB(x$obj, bias.correct = FALSE)
  se <- as.list(sr, "Std. Error", report = TRUE)
  expect_equal(x$data$est_se, as.numeric(se$proj_fe_combined))
  for (m in 1:2) {
    p <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, model = m)
    expect_equal(p$est_se, se$proj_fe[, m])
  }
})
