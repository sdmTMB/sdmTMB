# One non-spatial fit and one set of simulated predictions shared by most
# tests.
range_edge_grid <- function() replicate_df(qcs_grid_small, "year", unique(pcod$year))
range_edge_fit <- fit_once(function() {
  sdmTMB(
    density ~ 0 + as.factor(year),
    data = pcod, family = tweedie(link = "log"),
    time = "year", spatiotemporal = "off", spatial = "off"
  )
})
range_edge_sims <- fit_once(function() {
  set.seed(123)
  predict(range_edge_fit(), newdata = range_edge_grid(), nsim = 50)
})

test_that("get_range_edge() basic functionality works", {
  skip_on_cran()

  # Create prediction grid
  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Calculate range edges
  edges <- get_range_edge(p, axis = nd$Y)

  # Test output structure
  expect_s3_class(edges, "data.frame")
  expect_named(edges, c("year", "quantile", "est", "lwr", "upr", "se"))
  expect_equal(nrow(edges), length(unique(pcod$year)) * 2) # 2 quantiles by default
  expect_equal(unique(edges$quantile), c(0.025, 0.975))

  # Test that estimates are within reasonable range
  expect_true(all(edges$est >= min(nd$Y, na.rm = TRUE)))
  expect_true(all(edges$est <= max(nd$Y, na.rm = TRUE)))

  # Test that confidence intervals make sense
  expect_true(all(edges$lwr <= edges$est))
  expect_true(all(edges$upr >= edges$est))
  expect_true(all(edges$se >= 0))
})

test_that("get_range_edge() works with custom quantiles", {
  skip_on_cran()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Test with custom quantiles
  edges <- get_range_edge(p, axis = nd$Y, quantiles = c(0.05, 0.5, 0.95))

  expect_equal(unique(edges$quantile), c(0.05, 0.5, 0.95))
  expect_equal(nrow(edges), length(unique(pcod$year)) * 3)

  # Test that median is between lower and upper quantiles
  for (yr in unique(edges$year)) {
    yr_data <- edges[edges$year == yr, ]
    expect_true(yr_data$est[yr_data$quantile == 0.5] >= yr_data$est[yr_data$quantile == 0.05])
    expect_true(yr_data$est[yr_data$quantile == 0.5] <= yr_data$est[yr_data$quantile == 0.95])
  }
})

test_that("get_range_edge() works with return_sims = TRUE", {
  skip_on_cran()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Get simulation draws
  edges_sims <- get_range_edge(p, axis = nd$Y, return_sims = TRUE)

  # Test output structure
  expect_s3_class(edges_sims, "data.frame")
  expect_named(edges_sims, c("year", "quantile", ".value", ".iteration"))
  expect_equal(nrow(edges_sims), length(unique(pcod$year)) * 2 * ncol(p)) # 2 quantiles
  expect_equal(unique(edges_sims$.iteration), seq_len(ncol(p)))

  # Compare with summary output
  edges <- get_range_edge(p, axis = nd$Y, return_sims = FALSE)

  # Check that summary statistics are approximately correct
  for (yr in unique(edges$year)) {
    for (q in unique(edges$quantile)) {
      sim_vals <- edges_sims$.value[edges_sims$year == yr & edges_sims$quantile == q]
      summary_row <- edges[edges$year == yr & edges$quantile == q, ]

      # Median should be close to est
      expect_equal(median(sim_vals, na.rm = TRUE), summary_row$est, tolerance = 0.01)
      # SD should be close to se
      expect_equal(sd(sim_vals, na.rm = TRUE), summary_row$se, tolerance = 0.01)
    }
  }
})

test_that("get_range_edge() input validation works", {
  skip_on_cran()

  m <- range_edge_fit()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Test error for predictions without simulations
  p_no_sim <- predict(m, newdata = nd)
  expect_error(get_range_edge(p_no_sim, axis = nd$Y),
    regexp = "matrix output.*nsim > 0"
  )

  # Test error for wrong axis length
  expect_error(get_range_edge(p, axis = nd$Y[1:10]),
    regexp = "same length"
  )

  # Test error for non-numeric axis
  expect_error(get_range_edge(p, axis = rep("a", nrow(nd))),
    regexp = "numeric vector"
  )

  # Test error for invalid quantiles
  expect_error(get_range_edge(p, axis = nd$Y, quantiles = c(0, 1)),
    regexp = "between 0 and 1"
  )
  expect_error(get_range_edge(p, axis = nd$Y, quantiles = c(-0.1, 0.5)),
    regexp = "between 0 and 1"
  )
  expect_error(get_range_edge(p, axis = nd$Y, quantiles = "0.5"),
    regexp = "numeric"
  )

  # Test error for invalid level
  expect_error(get_range_edge(p, axis = nd$Y, level = 1.5),
    regexp = "level"
  )
  expect_error(get_range_edge(p, axis = nd$Y, level = 0),
    regexp = "level"
  )

  # Test error for invalid return_sims
  expect_error(get_range_edge(p, axis = nd$Y, return_sims = "TRUE"),
    regexp = "logical"
  )
})

test_that("get_range_edge() works with a spatial model, confidence levels, and edge cases", {
  skip_on_cran()

  m <- sdmTMB(
    density ~ depth_scaled,
    mesh = pcod_mesh_2011,
    data = pcod_2011, family = tweedie(link = "log"),
    time = "year", spatiotemporal = "off", spatial = "on"
  )
  # full grid: range edges on a coarse grid snap to a few Y values
  nd <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
  set.seed(202)
  p <- predict(m, newdata = nd, nsim = 50)

  # 70% CI should be narrower than 95% CI
  edges_70 <- get_range_edge(p, axis = nd$Y, level = 0.70)
  edges_95 <- get_range_edge(p, axis = nd$Y, level = 0.95)
  expect_true(mean(edges_70$upr - edges_70$lwr) < mean(edges_95$upr - edges_95$lwr))

  # Test with extreme quantiles
  edges <- get_range_edge(p, axis = nd$Y, quantiles = c(0.001, 0.999))
  expect_s3_class(edges, "data.frame")
  expect_true(all(!is.na(edges$est)))

  # Test with single quantile
  edges <- get_range_edge(p, axis = nd$Y, quantiles = 0.5)
  expect_equal(nrow(edges), length(unique(pcod_2011$year)))
  expect_equal(unique(edges$quantile), 0.5)
})

test_that("get_range_edge() works with different link functions", {
  skip_on_cran()

  # Test with binomial/logit link
  pcod$present <- as.numeric(pcod$density > 0)
  m_binomial <- sdmTMB(
    present ~ 1,
    data = pcod, family = binomial(),
    time = "year", spatiotemporal = "off", spatial = "off"
  )

  nd <- range_edge_grid()
  set.seed(404)
  p <- predict(m_binomial, newdata = nd, nsim = 50)

  edges <- get_range_edge(p, axis = nd$Y)

  # Should still produce valid output
  expect_s3_class(edges, "data.frame")
  expect_true(all(!is.na(edges$est)))
  expect_true(all(edges$lwr <= edges$est))
  expect_true(all(edges$upr >= edges$est))
})

test_that("get_range_edge() axis ordering is handled correctly", {
  skip_on_cran()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Test with Y axis (original)
  edges_y <- get_range_edge(p, axis = nd$Y)

  # Test with negative Y axis (reversed)
  edges_neg_y <- get_range_edge(p, axis = -nd$Y)

  # Range edges should be opposite in sign
  expect_equal(edges_y$est[edges_y$quantile == 0.025],
    -edges_neg_y$est[edges_neg_y$quantile == 0.975],
    tolerance = 0.1
  )

  # Test with X axis instead
  edges_x <- get_range_edge(p, axis = nd$X)

  # Should produce different results than Y axis
  expect_false(all(abs(edges_x$est - edges_y$est) < 1))
})

test_that("get_range_edge() warning for missing link attribute", {
  skip_on_cran()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  # Remove link attribute to trigger warning
  attr(p, "link") <- NULL

  expect_warning(get_range_edge(p, axis = nd$Y),
    regexp = "No link attribute"
  )
})

test_that("get_range_edge() handles consistent time column naming", {
  skip_on_cran()

  nd <- range_edge_grid()
  p <- range_edge_sims()

  edges <- get_range_edge(p, axis = nd$Y)

  # First column should be named "year" (the time column from the model)
  expect_equal(names(edges)[1], "year")
  expect_equal(unique(edges$year), sort(unique(pcod$year)))

  # Test with return_sims = TRUE
  edges_sims <- get_range_edge(p, axis = nd$Y, return_sims = TRUE)
  expect_equal(names(edges_sims)[1], "year")
})
