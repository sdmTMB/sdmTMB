# Tests for make_zero_one_map()

test_that("make_zero_one_map() detects all-zero and all-one levels (delta)", {
  d <- pcod
  d$density[d$year == 2003] <- 0   # force all-zero
  d$density[d$year == 2017] <- 1.5 # force all-positive

  m <- make_zero_one_map(
    formula = density ~ 0 + as.factor(year),
    data = d, group = "year", family = delta_gamma()
  )
  expect_true(2003 %in% m$all_zero_levels)
  expect_true(2017 %in% m$all_one_levels)

  yrs <- sort(unique(d$year))
  i_zero <- which(yrs == 2003)
  i_one  <- which(yrs == 2017)

  expect_true(is.na(as.integer(m$map$b_j)[i_zero]))
  expect_true(is.na(as.integer(m$map$b_j)[i_one]))
  expect_equal(m$start$b_j[i_zero], -20)
  expect_equal(m$start$b_j[i_one], 20)

  expect_true(is.na(as.integer(m$map$b_j2)[i_zero]))
  expect_false(is.na(as.integer(m$map$b_j2)[i_one]))
  expect_equal(m$start$b_j2[i_zero], -20)
  expect_equal(m$start$b_j2[i_one], 0)
})

test_that("make_zero_one_map() works for non-delta (Tweedie) families", {
  d <- pcod
  d$density[d$year == 2003] <- 0
  d$density[d$year == 2017] <- 1.5

  m <- make_zero_one_map(
    formula = density ~ 0 + as.factor(year),
    data = d, group = "year", family = tweedie()
  )
  expect_null(m$map$b_j2)
  expect_null(m$start$b_j2)
  expect_true(2003 %in% m$all_zero_levels)
  expect_length(m$all_one_levels, 0)
  yrs <- sort(unique(d$year))
  i_one <- which(yrs == 2017)
  expect_false(is.na(as.integer(m$map$b_j)[i_one]))
  expect_equal(m$start$b_j[i_one], 0)
})

test_that("level ordering matches factor level order even when data are shuffled", {
  d <- pcod
  d$density[d$year == 2003] <- 0
  d <- d[order(-d$year), ]

  m <- make_zero_one_map(
    formula = density ~ 0 + as.factor(year),
    data = d, group = "year", family = tweedie()
  )
  yrs <- sort(unique(d$year))
  i_zero <- which(yrs == 2003)
  expect_equal(m$start$b_j[i_zero], -20)
  expect_true(is.na(as.integer(m$map$b_j)[i_zero]))
  expect_equal(sum(is.na(as.integer(m$map$b_j))), 1L)
})

test_that("works on a non-time grouping factor (region)", {
  d <- pcod
  d$region <- ifelse(d$Y > mean(d$Y), "N", "S")
  # Force the "S" region to be all zeros
  d$density[d$region == "S"] <- 0

  m <- make_zero_one_map(
    formula = density ~ 0 + as.factor(region),
    data = d, group = "region", family = tweedie()
  )
  expect_true("S" %in% m$all_zero_levels)
  lvls <- sort(unique(d$region))
  i_s <- which(lvls == "S")
  expect_equal(m$start$b_j[i_s], -20)
  expect_true(is.na(as.integer(m$map$b_j)[i_s]))
})

test_that("works on a region-year interaction column", {
  d <- pcod
  d$region <- ifelse(d$Y > mean(d$Y), "N", "S")
  d$ry <- interaction(d$region, d$year, drop = TRUE)
  # Force one region-year cell to all zero
  d$density[d$region == "S" & d$year == 2007] <- 0

  m <- make_zero_one_map(
    formula = density ~ 0 + ry,
    data = d, group = "ry", family = tweedie()
  )
  expect_true("S.2007" %in% as.character(m$all_zero_levels))
})

test_that("NA-filtered rows match the fitting model frame", {
  d <- data.frame(
    y = c(0, 1, 1),
    year = c(2001, 2001, 2002),
    x = c(1, NA, 2)
  )

  m <- make_zero_one_map(
    formula = y ~ 0 + as.factor(year) + x,
    data = d, group = "year", family = tweedie()
  )

  expect_identical(m$all_zero_levels, 2001)
  expect_true(is.na(as.integer(m$map$b_j)[1]))
  expect_equal(m$start$b_j[1], -20)
})

test_that("unused factor levels are ignored", {
  d <- data.frame(
    y = c(0, 1),
    year = factor(c("2001", "2001"), levels = c("2001", "2002", "2003"))
  )

  m <- make_zero_one_map(
    formula = y ~ 0 + as.factor(year),
    data = d, group = "year", family = tweedie()
  )

  expect_length(m$all_zero_levels, 0L)
  expect_length(m$all_one_levels, 0L)
  expect_false(any(is.na(as.integer(m$map$b_j))))
})

test_that("weighted binomial proportions detect all-zero and all-one levels", {
  d <- data.frame(
    prop = c(1, 1, 0),
    year = c(2001, 2001, 2002)
  )
  w <- c(3, 4, 5)

  m <- make_zero_one_map(
    formula = prop ~ 0 + as.factor(year),
    data = d, group = "year", family = binomial(), weights = w
  )

  expect_true(2001 %in% m$all_one_levels)
  expect_true(2002 %in% m$all_zero_levels)
  expect_true(all(is.na(as.integer(m$map$b_j))))
  expect_equal(m$start$b_j, c(20, -20))
})

test_that("Poisson-link delta warns/informs for all-positive levels", {
  d <- pcod
  d$density[d$year == 2017] <- 1.5

  expect_message(
    make_zero_one_map(
      formula = density ~ 0 + as.factor(year),
      data = d, group = "year",
      family = delta_gamma(type = "poisson-link")
    ),
    regexp = "Poisson-link"
  )
})

test_that("input validation: missing columns and bad inputs", {
  d <- pcod
  expect_error(
    make_zero_one_map(density ~ 0 + as.factor(year), d, "nope", tweedie()),
    "not found"
  )
  expect_error(
    make_zero_one_map(missing_resp ~ 0 + as.factor(year), d, "year", tweedie()),
    "not found"
  )
  expect_error(
    make_zero_one_map(
      density ~ 0 + as.factor(year), d, "year", binomial(),
      weights = rep(1, nrow(d) - 1L)
    )
  )
  expect_error(
    make_zero_one_map(density ~ 0 + as.factor(year), d, "year", tweedie(), value = -1)
  )
  expect_error(
    make_zero_one_map(
      list(density ~ 0 + as.factor(year),
           density ~ 0 + as.factor(year),
           density ~ 0 + as.factor(year)),
      d, "year", delta_gamma()
    )
  )
})

test_that("design without matching factor column errors informatively", {
  d <- pcod
  d$density[d$year == 2003] <- 0
  expect_error(
    make_zero_one_map(
      formula = density ~ 1 + depth,
      data = d, group = "year", family = tweedie()
    ),
    "design matrix"
  )
})

test_that("returned map/start can be used in sdmTMBcontrol() and sdmTMB() fits", {
  skip_on_cran()
  d <- pcod
  d$density[d$year == 2003] <- 0

  m <- make_zero_one_map(
    formula = density ~ 0 + as.factor(year),
    data = d, group = "year", family = tweedie()
  )

  fit <- sdmTMB(
    density ~ 0 + as.factor(year),
    time = "year",
    data = d, family = tweedie(), spatial = "off", spatiotemporal = "off",
    control = sdmTMBcontrol(map = m$map, start = m$start)
  )
  expect_s3_class(fit, "sdmTMB")
  yrs <- sort(unique(d$year))
  i_zero <- which(yrs == 2003)
  expect_equal(unname(fit$tmb_obj$env$parList()$b_j[i_zero]), -20)

  # index?
  p <- predict(
    fit,
    newdata = replicate_df(qcs_grid, "year", unique(pcod$year)),
    return_tmb_object = TRUE
  )
  ind <- get_index(p, bias_correct = FALSE)
  expect_true(ind$est[ind$year == 2003] < 0.0001)
  expect_true(ind$se[ind$year == 2003] < 0.000001)
})
