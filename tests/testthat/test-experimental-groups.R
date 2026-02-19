test_that("experimental groups is parsed and plumbed", {
  local_edition(2)

  d <- subset(pcod, year %in% c(2003, 2005, 2007))
  d <- d[seq_len(min(300L, nrow(d))), ]
  d$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 30)

  fit <- sdmTMB(
    density ~ 1,
    data = d,
    time = "year",
    mesh = mesh,
    family = tweedie(link = "log"),
    spatial = "off",
    spatiotemporal = "rw",
    do_fit = FALSE,
    experimental = list(groups = "grp")
  )
  expect_identical(fit$groups, "grp")
  expect_identical(fit$group_levels, levels(d$grp))
  expect_equal(fit$tmb_data$n_groups, nlevels(d$grp))
  expect_identical(fit$tmb_data$group_i, as.integer(d$grp) - 1L)
  expect_true("upsilon_stg" %in% fit$tmb_random)
  expect_false("epsilon_st" %in% fit$tmb_random)
})

test_that("experimental groups validates input", {
  local_edition(2)
  d <- pcod_2011
  d$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = binomial(),
      do_fit = FALSE,
      experimental = list(groups = 1)
    ),
    regexp = "single column name"
  )

  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = binomial(),
      do_fit = FALSE,
      experimental = list(groups = "missing")
    ),
    regexp = "must match a column"
  )

  d2 <- d
  d2$grp_chr <- as.character(d2$grp)
  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d2,
      mesh = pcod_mesh_2011,
      family = binomial(),
      do_fit = FALSE,
      experimental = list(groups = "grp_chr")
    ),
    regexp = "must be a factor"
  )

  d3 <- d
  d3$grp[1] <- NA
  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d3,
      mesh = pcod_mesh_2011,
      family = binomial(),
      do_fit = FALSE,
      experimental = list(groups = "grp")
    ),
    regexp = "cannot contain `NA` values"
  )

  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = binomial(),
      do_fit = FALSE,
      experimental = list(groups = "grp"),
      spatiotemporal = "off"
    ),
    regexp = "requires `spatiotemporal = 'rw'`"
  )

  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = delta_gamma(),
      time = "year",
      do_fit = FALSE,
      experimental = list(groups = "grp"),
      spatiotemporal = "rw"
    ),
    regexp = "not supported for delta models"
  )

  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = binomial(),
      time = "year",
      do_fit = FALSE,
      do_index = TRUE,
      experimental = list(groups = "grp"),
      spatiotemporal = "rw"
    ),
    regexp = "not currently compatible with `do_index = TRUE`"
  )

  d$time_std <- as.numeric(scale(d$year))
  expect_error(
    sdmTMB(
      present ~ depth_scaled,
      data = d,
      mesh = pcod_mesh_2011,
      family = binomial(),
      time = "year",
      do_fit = FALSE,
      experimental = list(groups = "grp", epsilon_model = "trend", epsilon_predictor = "time_std"),
      spatiotemporal = "rw"
    ),
    regexp = "not currently compatible with `experimental\\$epsilon_model`"
  )
})

test_that("grouped RW C++ path fits as experimental", {
  local_edition(2)
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  d <- data.frame(
    X = runif(120),
    Y = runif(120),
    year = rep(1:4, each = 30)
  )
  d$grp <- factor(rep(c("a", "b"), each = nrow(d) / 2))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.15)
  s <- sdmTMB_simulate(
    formula = ~1,
    data = d,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    range = 0.4,
    sigma_E = 0.2,
    phi = 0.1,
    sigma_O = 0,
    spatiotemporal = "rw",
    seed = 1,
    B = 0
  )
  s$grp <- d$grp

  fit <- sdmTMB(
    observed ~ 1,
    data = s,
    time = "year",
    mesh = mesh,
    family = gaussian(),
    spatiotemporal = "rw",
    extra_time = 5L,
    experimental = list(groups = "grp")
  )
  expect_s3_class(fit, "sdmTMB")
  expect_true("upsilon_stg" %in% fit$tmb_random)
  est <- as.list(fit$sd_report, "Estimate", report = TRUE)
  expect_true("sigma_U" %in% names(est))

  p_null <- predict(fit)
  p_data <- predict(fit, newdata = s)
  expect_equal(p_null$est, p_data$est, tolerance = 1e-6)
  expect_equal(p_null$est_rf, p_data$est_rf, tolerance = 1e-6)
  expect_equal(p_null$epsilon_st, p_data$epsilon_st, tolerance = 1e-6)

  p_tmb <- predict(fit, newdata = s, return_tmb_object = TRUE)
  expect_identical(p_tmb$pred_tmb_data$proj_group_i, as.integer(s$grp) - 1L)

  nd_missing_group <- s
  nd_missing_group$grp <- NULL
  expect_error(
    predict(fit, newdata = nd_missing_group),
    regexp = "Missing column `grp`"
  )

  nd_extra_level <- s
  nd_extra_level$grp <- factor(as.character(nd_extra_level$grp), levels = c(levels(s$grp), "c"))
  expect_error(
    predict(fit, newdata = nd_extra_level),
    regexp = "not identical between fitted and prediction data"
  )

  nd_reordered_levels <- s
  nd_reordered_levels$grp <- factor(as.character(nd_reordered_levels$grp), levels = rev(levels(s$grp)))
  expect_error(
    predict(fit, newdata = nd_reordered_levels),
    regexp = "not identical between fitted and prediction data"
  )

  nd_subset <- subset(s, grp == "a")
  nd_subset$grp <- factor(nd_subset$grp, levels = levels(s$grp))
  p_subset <- predict(fit, newdata = nd_subset)
  expect_equal(p_subset$est, p_data$est[s$grp == "a"], tolerance = 1e-6)

  nd_extra_time <- replicate_df(s, "year", 1:5)
  nd_extra_time$grp <- factor(nd_extra_time$grp, levels = levels(s$grp))
  p_extra_time <- predict(fit, newdata = nd_extra_time)
  expect_equal(nrow(p_extra_time), nrow(nd_extra_time))
})

test_that("make_groups validates prediction levels against fit levels", {
  x <- factor(c("a", "b", "a"), levels = c("a", "b"))
  expect_identical(make_groups(x, prev_levels = c("a", "b")), c(0L, 1L, 0L))

  expect_error(
    make_groups(factor(c("a", "b"), levels = c("a", "b", "unused"))),
    regexp = "Extra factor levels found"
  )

  expect_error(
    make_groups(c("a", "b")),
    regexp = "not a factor"
  )

  expect_error(
    make_groups(factor(c("a", NA), levels = c("a", "b"))),
    regexp = "cannot contain `NA` values"
  )

  x_bad <- factor(c("a", "b"), levels = c("b", "a"))
  expect_error(
    make_groups(x_bad, prev_levels = c("a", "b")),
    regexp = "not identical between fitted and prediction data"
  )
})
