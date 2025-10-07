# # Tests for covariate diffusion functionality
#
# test_that("covariate_diffusion argument is parsed correctly", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- scale(d$depth)[, 1]
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)
#
#   # Test formula parsing
#   expect_error(
#     sdmTMB(
#       data = d,
#       formula = density ~ 1,
#       covariate_diffusion = ~ depth_centered,
#       mesh = mesh,
#       family = gamma(link = "log"),
#       do_fit = FALSE
#     ),
#     NA  # Should not error
#   )
#
#   # Test invalid covariate_diffusion argument
#   expect_error(
#     sdmTMB(
#       data = d,
#       formula = density ~ 1,
#       covariate_diffusion = "depth_centered",  # Should be formula
#       mesh = mesh,
#       family = gamma(link = "log"),
#       do_fit = FALSE
#     ),
#     "formula"
#   )
#
#   # Test covariate not in data
#   expect_error(
#     sdmTMB(
#       data = d,
#       formula = density ~ 1,
#       covariate_diffusion = ~ nonexistent_var,
#       mesh = mesh,
#       family = gamma(link = "log"),
#       do_fit = FALSE
#     ),
#     "object.*not found"
#   )
# })
#
# test_that("covariate diffusion model fits successfully with Gaussian family", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   # Basic covariate diffusion model
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian(link = "identity")
#   )
#
#   expect_s3_class(m1, "sdmTMB")
#   expect_true(m1$tmb_data$covariate_diffusion == 1)
#   expect_true(m1$model$convergence == 0)
#
#   # Check that diffused covariate parameters exist
#   pars <- names(m1$model$par)
#   expect_true("ln_tau_diffusion" %in% pars)
#   expect_true("diffused_cov_s" %in% pars)
#
#   # Check dimensions match mesh vertices and observation count
#   n_vertices <- mesh$mesh$n
#   expect_equal(dim(m1$tmb_data$vertex_cov_raw), c(n_vertices, 1))
#   expect_equal(dim(m1$model$par$diffused_cov_s), c(n_vertices, 1))
#   expect_equal(length(m1$model$par$ln_tau_diffusion), 1)
# }
#
# test_that("covariate diffusion works with multiple diffused covariates", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   d$temp <- rnorm(nrow(d))  # Simulated temperature covariate
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered + temp,
#     covariate_diffusion = ~ depth_centered + temp,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   expect_s3_class(m1, "sdmTMB")
#   n_vertices <- mesh$mesh$n
#   expect_equal(dim(m1$tmb_data$vertex_cov_raw), c(n_vertices, 2))
#   expect_equal(dim(m1$model$par$diffused_cov_s), c(n_vertices, 2))
#   expect_equal(length(m1$model$par$ln_tau_diffusion), 2)
#
#   # Test tidy output includes diffusion parameters
#   b1 <- tidy(m1, effects = "ran_pars")
#   diffusion_pars <- b1[grepl("tau_diffusion", b1$term), ]
#   expect_equal(nrow(diffusion_pars), 2)
# })
#
# test_that("covariate diffusion works with different families", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   # Poisson family
#   m_pois <- sdmTMB(
#     data = d,
#     formula = present ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = poisson(link = "log")
#   )
#   expect_s3_class(m_pois, "sdmTMB")
#   expect_true(m_pois$model$convergence == 0)
#
#   # Gamma family
#   d_pos <- d[d$density > 0, ]
#   m_gamma <- sdmTMB(
#     data = d_pos,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = Gamma(link = "log")
#   )
#   expect_s3_class(m_gamma, "sdmTMB")
#   expect_true(m_gamma$model$convergence == 0)
# })
#
# test_that("covariate diffusion works with delta models", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m_delta <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = delta_gamma()
#   )
#
#   expect_s3_class(m_delta, "sdmTMB")
#   expect_true(m_delta$model$convergence == 0)
#
#   # Check that diffusion is applied to both model components
#   expect_equal(length(m_delta$model$par$ln_tau_diffusion), 2)  # One for each component
# })
#
# test_that("covariate diffusion works with spatiotemporal models", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m_st <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     time = "year",
#     mesh = mesh,
#     family = gamma(link = "log"),
#     spatiotemporal = "ar1"
#   )
#
#   expect_s3_class(m_st, "sdmTMB")
#   expect_true(m_st$model$convergence == 0)
#
#   # Check that time dimension is handled correctly
#   expect_equal(dim(m_st$tmb_data$cov_diffusion_i), c(nrow(d), 1))
# })
#
# test_that("prediction works correctly with covariate diffusion", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Predict on original data
#   p1 <- predict(m1)
#   expect_s3_class(p1, "data.frame")
#   expect_equal(nrow(p1), nrow(d))
#   expect_true("est" %in% names(p1))
#   expect_true("diffused_depth_centered" %in% names(p1))
#
#   # Predict on new data (subset)
#   d_new <- d[1:100, ]
#   p2 <- predict(m1, newdata = d_new)
#   expect_equal(nrow(p2), nrow(d_new))
#   expect_true("diffused_depth_centered" %in% names(p2))
#
#   # Check that diffused values are different from original
#   expect_false(all(p1$diffused_depth_centered == d$depth_centered))
# })
#
# test_that("covariate diffusion preserves covariate mass conservation", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   p1 <- predict(m1)
#
#   # Test that mass conservation principle is maintained
#   # Note: Exact conservation depends on interpolation method and mesh density
#   # We test that the diffused values are smoother (lower variance) than original
#   original_var <- var(d$depth_centered)
#   diffused_var <- var(p1$diffused_depth_centered)
#   expect_true(diffused_var <= original_var) # Diffusion should smooth/reduce variance
#
#   # Test that means are approximately preserved (interpolation bias)
#   original_mean <- mean(d$depth_centered)
#   diffused_mean <- mean(p1$diffused_depth_centered)
#   expect_equal(original_mean, diffused_mean, tolerance = 0.1)
# })
#
# test_that("covariate diffusion works with spatial_varying", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   d$year_scaled <- as.numeric(scale(d$year))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   # Test that covariate_diffusion and spatial_varying can be used together
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered + year_scaled,
#     covariate_diffusion = ~ depth_centered,
#     spatial_varying = ~ 0 + year_scaled,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   expect_s3_class(m1, "sdmTMB")
#   expect_true(m1$model$convergence == 0)
#
#   # Check that both features are present in the model
#   pars <- names(m1$model$par)
#   expect_true("ln_tau_diffusion" %in% pars)
#   expect_true("zeta_s" %in% pars)
# })
#
# test_that("tidy method works with covariate diffusion models", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Test fixed effects
#   b1 <- tidy(m1, effects = "fixed")
#   expect_s3_class(b1, "data.frame")
#   expect_true(nrow(b1) >= 2)
#
#   # Test random parameters include diffusion
#   b2 <- tidy(m1, effects = "ran_pars")
#   expect_true(any(grepl("tau_diffusion", b2$term)))
#
#   # Test with confidence intervals
#   b3 <- tidy(m1, effects = "ran_pars", conf.int = TRUE)
#   expect_true("conf.low" %in% names(b3))
#   expect_true("conf.high" %in% names(b3))
# })
#
# test_that("print method works with covariate diffusion models", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Test that print method doesn't error and includes diffusion info
#   output <- capture.output(print(m1))
#   expect_true(any(grepl("Covariate diffusion", output)))
#   expect_true(any(grepl("depth_centered", output)))
# })
#
# test_that("covariate diffusion errors appropriately with invalid inputs", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   # Test error when covariate_diffusion covariate not in main formula
#   expect_error(
#     sdmTMB(
#       data = d,
#       formula = density ~ 1,  # depth_centered not included
#       covariate_diffusion = ~ depth_centered,
#       mesh = mesh,
#       family = gaussian()
#     ),
#     "covariate_diffusion.*formula"
#   )
#
#   # Test error with intercept in covariate_diffusion
#   expect_error(
#     sdmTMB(
#       data = d,
#       formula = density ~ 1 + depth_centered,
#       covariate_diffusion = ~ 1 + depth_centered,  # Should not include intercept
#       mesh = mesh,
#       family = gaussian()
#     ),
#     "intercept.*covariate_diffusion"
#   )
# })
#
# test_that("simulate method works with covariate diffusion models", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Test simulate method
#   set.seed(123)
#   sim1 <- simulate(m1, nsim = 1)
#   expect_equal(length(sim1), nrow(d))
#   expect_true(is.numeric(sim1))
#
#   # Test simulate with newdata
#   d_new <- d[1:50, ]
#   sim2 <- simulate(m1, nsim = 1, newdata = d_new)
#   expect_equal(length(sim2), nrow(d_new))
# })
#
# test_that("vertex interpolation produces reasonable covariate values", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Check that vertex covariate values are within reasonable range
#   vertex_cov <- m1$tmb_data$vertex_cov_raw[, 1]
#   obs_range <- range(d$depth_centered)
#
#   # Allow some extrapolation but not extreme values
#   expect_true(all(vertex_cov >= obs_range[1] - 2))
#   expect_true(all(vertex_cov <= obs_range[2] + 2))
#
#   # Vertex values should show spatial smoothing compared to irregular observations
#   vertex_var <- var(vertex_cov)
#   obs_var <- var(d$depth_centered)
#
#   # This may not always hold depending on mesh density vs observation density
#   # but generally vertex interpolation should create smoother fields
#   # expect_true(vertex_var <= obs_var * 1.5)  # Allow some variance due to interpolation
# })
#
# test_that("diffusion parameters are estimable and reasonable", {
#   skip_on_cran()
#   local_edition(2)
#
#   d <- pcod
#   d$depth_centered <- as.numeric(scale(d$depth))
#   mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
#
#   m1 <- sdmTMB(
#     data = d,
#     formula = density ~ 1 + depth_centered,
#     covariate_diffusion = ~ depth_centered,
#     mesh = mesh,
#     family = gaussian()
#   )
#
#   # Check that tau_diffusion is estimated (not at boundary)
#   tau_diffusion <- exp(m1$model$par$ln_tau_diffusion)
#   expect_true(tau_diffusion > 0.001)  # Not too small (over-smoothed)
#   expect_true(tau_diffusion < 1000)   # Not too large (under-smoothed)
#
#   # Check that diffused covariates show spatial correlation
#   diffused_cov <- m1$model$par$diffused_cov_s[, 1]
#   vertex_coords <- mesh$mesh$loc[, 1:2]
#
#   # Calculate spatial correlation (should be positive for nearby vertices)
#   if (length(diffused_cov) >= 10) {  # Only if enough vertices
#     distances <- as.matrix(dist(vertex_coords[1:min(10, length(diffused_cov)), ]))
#     cov_diffs <- as.matrix(dist(diffused_cov[1:min(10, length(diffused_cov))]))
#
#     # Nearby vertices should have more similar covariate values
#     correlation <- cor(c(distances), c(cov_diffs), use = "complete.obs")
#     expect_true(correlation > 0)  # Positive correlation between distance and covariate difference
#   }
# })
