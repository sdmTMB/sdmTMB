test_that(".apply_design() reuses fitted levels, contrasts, and bases", {
  d <- data.frame(
    x = seq(-2, 2, length.out = 12),
    f = factor(rep(c("a", "b", "c"), 4)),
    ch = rep(c("u", "v"), 6)
  )
  des <- .make_design(~ f + ch + poly(x, 2), d, "`test`")
  expect_identical(des$spec$columns, colnames(des$X))

  # subset and reordered rows match the fitted rows
  idx <- c(7, 2, 11, 2)
  expect_equal(.apply_design(des$spec, d[idx, ]), des$X[idx, ],
    ignore_attr = TRUE)

  # a single row uses the fitted poly() basis
  expect_equal(.apply_design(des$spec, d[5, ]), des$X[5, , drop = FALSE],
    ignore_attr = TRUE)

  # absent factor levels keep all fitted columns
  nd <- d[d$f == "b", ]
  X_nd <- .apply_design(des$spec, nd)
  expect_identical(colnames(X_nd), des$spec$columns)
  expect_equal(X_nd, des$X[d$f == "b", ], ignore_attr = TRUE)

  # global contrasts don't change the saved encoding
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op), add = TRUE)
  expect_equal(.apply_design(des$spec, d[idx, ]), des$X[idx, ],
    ignore_attr = TRUE)
  options(op)

  # new levels and missing values are errors
  nd <- d[1:2, ]
  nd$f <- factor(c("a", "z"))
  expect_error(.apply_design(des$spec, nd), "`test`")
  nd <- d[1:2, ]
  nd$x[2] <- NA
  expect_error(.apply_design(des$spec, nd), "Missing values")
  expect_error(.apply_design(des$spec, d[, c("f", "ch")]), "`test`")
})

test_that("Auxiliary prediction designs reuse the fitted encoding", {
  skip_on_cran()
  d <- pcod_2011
  d$fyear <- factor(d$year)
  set.seed(1)
  d$reg <- factor(sample(c("n", "m", "s"), nrow(d), replace = TRUE))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 30)
  op <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(op), add = TRUE)
  # The designs are under test, not convergence; the random `reg` effects
  # are not identifiable, so suppress the resulting fit warnings.
  fit <- suppressWarnings(sdmTMB(
    log(density + 1) ~ 1,
    data = d, mesh = mesh, time = "year",
    spatiotemporal = "off",
    spatial_varying = ~ 1 + reg,
    time_varying = ~ 0 + poly(depth_scaled, 2),
    dispformula = ~ reg + poly(depth_scaled, 2),
    silent = TRUE
  ))
  # SVC intercept is aliased to omega_s and dropped:
  expect_identical(fit$design_specs$spatial_varying$columns, c("regn", "regs"))
  expect_identical(colnames(fit$tmb_data$z_i), c("regn", "regs"))

  check_rows <- function(fit, idx) {
    nd <- d[idx, ]
    expect_equal(predict_svc_matrix(fit, nd), fit$tmb_data$z_i[idx, , drop = FALSE],
      ignore_attr = TRUE)
    expect_equal(predict_tv_matrix(fit, nd), fit$tmb_data$X_rw_ik[idx, , drop = FALSE],
      ignore_attr = TRUE)
    expect_equal(predict_disp_matrix(fit, nd), fit$tmb_data$Xdisp_ij[idx, , drop = FALSE],
      ignore_attr = TRUE)
    expect_identical(colnames(predict_svc_matrix(fit, nd)), c("regn", "regs"))
  }
  idx <- c(50, 3, 200, 3)
  check_rows(fit, idx)
  check_rows(fit, 17) # one row with the fitted poly() basis
  check_rows(fit, which(d$reg == "s")) # absent factor levels

  p_fit <- predict(fit)
  p_nd <- predict(fit, newdata = d[idx, ])
  expect_equal(p_nd$est, p_fit$est[idx], tolerance = 1e-8)

  # switching global contrasts after fitting doesn't change predictions
  options(contrasts = c("contr.sum", "contr.poly"))
  check_rows(fit, idx)
  expect_equal(predict(fit, newdata = d[idx, ])$est, p_nd$est, tolerance = 1e-8)
  options(op)

  nd <- d[1:3, ]
  nd$reg <- factor(c("n", "m", "new"))
  expect_error(predict(fit, newdata = nd), "new levels")

  # the RTMB backend projects the same saved designs
  rt <- fit
  rt$backend <- "rtmb"
  expect_equal(predict(rt, newdata = d[idx, ])$est, p_nd$est, tolerance = 1e-8)
  expect_rtmb_fit_data_matches(fit, as.data.frame(d[idx, ]), info = "design specs")

  # fits saved before `design_specs` reconstruct the designs
  legacy <- fit
  legacy$design_specs <- NULL
  check_rows(legacy, idx)
  expect_equal(predict(legacy, newdata = d[idx, ])$est, p_nd$est, tolerance = 1e-8)
})
