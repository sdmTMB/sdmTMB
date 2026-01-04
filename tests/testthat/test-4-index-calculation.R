test_that("get_index(), get_index_sims(), and get_cog() work", {
  local_edition(3)
  skip_on_cran()
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    time = "year", mesh = pcod_spde, family = tweedie(link = "log")
  )
  # expect_snapshot(m)
  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)

  p <- predict(m, newdata = nd, return_tmb_object = FALSE)
  expect_error(get_index(p), regexp = "return_tmb_object")

  ind <- get_index(predictions, bias_correct = FALSE)
  expect_equal(class(ind), "data.frame")

  ind_corrected <- get_index(predictions, bias_correct = TRUE)
  cached <- c(259275.869, 385596.993, 428192.403, 116533.436, 189156.419,
    333697.54, 336522.903, 400682.008, 191430.747)
  expect_lt(sum(abs(ind_corrected$est - cached) / cached), 1e-03)
  expect_gt(mean(ind_corrected$est - ind$est), 0)

  cached <- c(209276.6742, 316827.7388, 355482.8274, 91541.2387, 150187.3753,
    274230.3944, 274541.7833, 326960.9376, 151502.4236)
  expect_lt(sum(abs(ind$est - cached) / cached), 1e-03)

  cached <- c(151512.8161, 244649.3652, 274277.5841, 66923.457, 109088.6819,
    212808.0659, 204421.8674, 241466.6371, 110227.6953)
  expect_lt(sum(abs(ind$lwr - cached) / cached), 1e-03)

  set.seed(1)
  pred_sim <- predict(m, nsim = 2000L)
  ind_sim <- get_index_sims(pred_sim)
  expect_gt(cor(ind_sim$est, ind$est), 0.9)
  expect_gt(cor(ind_sim$lwr, ind$lwr), 0.9)
  expect_gt(cor(ind_sim$upr, ind$upr), 0.9)
  # sims mimics bias corrected index, which would be higher:
  expect_gt(mean(ind$est), mean(ind_sim$est))

  cog <- get_cog(predictions, bias_correct = FALSE)
  expect_equal(class(cog), "data.frame")
  cached <- c(462.0865, 481.3353, 471.6714, 481.9013, 485.464, 469.7996,
    475.9824, 457.4335, 463.3973, 5756.9909, 5729.2348, 5760.9857,
    5735.0945, 5726.1765, 5745.2497, 5744.2829, 5757.1354, 5755.7778
  )
  expect_lt(sum(abs(cog$est - cached) / cached), 1e-04)

  cog_wide <- get_cog(predictions, bias_correct = FALSE, format="wide")
  expect_equal(class(cog_wide), "data.frame")
  expect_equal(names(cog_wide), c("year", "est_x", "lwr_x","upr_x", "se_x","est_y","lwr_y","upr_y","se_y", "type"))
  expect_equal(cog$est[which(cog$coord=="X")], cog_wide$est_x)
})

test_that("get_cog() with bias correction works", {
  skip_on_cran()
  skip_on_ci() # a little slow
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    time = "year", mesh = pcod_spde, family = tweedie(link = "log")
  )

  nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
  predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
  system.time({
    cog_bc <- get_cog(predictions, bias_correct = TRUE)
  })
  cached <- c(463.4799, 482.7737, 471.7401, 481.6446, 485.0965, 470.2728,
    475.3523, 460.7445, 465.8458, 5755.4058, 5727.3907, 5758.6983,
    5733.3987, 5726.0713, 5744.0075, 5745.9221, 5753.5249, 5752.7965
  )
  expect_lt(sum(abs(cog_bc$est - cached) / cached), 1e-04)
})

test_that("index errors are returned as needed", {
  skip_on_cran()

  g <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))

  expect_error(
  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011,
    spatial = "off", spatiotemporal = "off",
    family = tweedie(link = "log"),
    time = "year",
    predict_args = list(newdata = g),
    index_args = list(area = 1)
  ), regexp = "do_index" # missing!
  )

  fit <- sdmTMB(
    density ~ 1,
    data = pcod_2011, spatial = "off", spatiotemporal = "off",
    family = tweedie(link = "log"),
    time = "year"
  )
  p1 <- predict(fit, newdata = NULL, return_tmb_object = TRUE)
  expect_error(get_index(p1), "newdata") # missing!

  p2 <- predict(fit, newdata = g, return_tmb_object = TRUE)
  suppressMessages(
    i <- get_index(p2)
  )
  expect_s3_class(i, "data.frame")
})
