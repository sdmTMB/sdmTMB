#!/usr/bin/env Rscript

# Test script for covariate diffusion functionality
library(sdmTMB)

# Load data
d <- pcod

# Create mesh
mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)

cat("Testing covariate diffusion functionality...\n")

# First test: Check that we get proper error when mesh doesn't have vertex_covariates
cat("\n1. Testing error when mesh lacks vertex_covariates...\n")
tryCatch({
  m_error <- sdmTMB(
    data = d,
    formula = density ~ 1 + depth_scaled,
    covariate_diffusion = ~ depth_scaled,
    mesh = mesh,
    family = gaussian(),
    do_fit = FALSE
  )
  cat("ERROR: Should have failed but didn't!\n")
}, error = function(e) {
  cat("PASS: Got expected error:", e$message, "\n")
})

# Add vertex covariates to mesh using tinyVAST
cat("\n2. Adding vertex covariates to mesh using tinyVAST::add_mesh_covariates()...\n")
if (!requireNamespace("tinyVAST", quietly = TRUE)) {
  # If tinyVAST is not available, create mock vertex_covariates
  cat("WARNING: tinyVAST not available, creating mock vertex_covariates\n")
  set.seed(123)
  n_vertices <- mesh$mesh$n
  mesh$vertex_covariates <- data.frame(
    depth_scaled = rnorm(n_vertices, mean = mean(d$depth_scaled), sd = sd(d$depth_scaled))
  )
} else {
  mesh <- tinyVAST::add_mesh_covariates(
    mesh = mesh,
    data = sdmTMB::qcs_grid,
    covariates = "depth_scaled",
    coords = c("X", "Y")
  )
}
cat("Mesh now has vertex_covariates with", ncol(mesh$vertex_covariates), "variables\n")
cat("Variable names:", names(mesh$vertex_covariates), "\n")

# Test basic covariate diffusion model setup
cat("\n3. Testing covariate diffusion TMB setup (do_fit = FALSE)...\n")
tryCatch({
  m1 <- sdmTMB(
    data = d,
    formula = density ~ 1 + depth_scaled,
    covariate_diffusion = ~ depth_scaled,
    mesh = mesh,
    family = gaussian(),
    do_fit = FALSE
  )
  cat("SUCCESS: Covariate diffusion TMB setup completed!\n")

  # Inspect TMB data structure
  cat("\nTMB data elements:\n")
  cat("covariate_diffusion flag:", m1$tmb_data$covariate_diffusion, "\n")
  cat("n_cov_diffusion:", m1$tmb_data$n_cov_diffusion, "\n")
  cat("vertex_cov_raw dimensions:", dim(m1$tmb_data$vertex_cov_raw), "\n")
  cat("diffusion_col_indices:", m1$tmb_data$diffusion_col_indices, "\n")
  cat("invM0 dimensions:", dim(m1$tmb_data$invM0), "\n")
  cat("invsqrtM0 dimensions:", dim(m1$tmb_data$invsqrtM0), "\n")
  cat("M1 dimensions:", dim(m1$tmb_data$M1), "\n")

  # Inspect TMB parameters
  cat("\nTMB parameters:\n")
  cat("ln_tau_diffusion length:", length(m1$tmb_params$ln_tau_diffusion), "\n")
  cat("diffused_cov_s dimensions:", dim(m1$tmb_params$diffused_cov_s), "\n")

  # Check parameter mapping
  cat("\nParameter mapping:\n")
  cat("ln_tau_diffusion mapped:", !is.null(m1$tmb_map$ln_tau_diffusion), "\n")
  cat("diffused_cov_s mapped:", !is.null(m1$tmb_map$diffused_cov_s), "\n")

  # Now try fitting
  cat("\n3b. Now attempting to fit the model...\n")
  m1_fit <- sdmTMB(
    data = d,
    formula = density ~ 1 + depth_scaled,
    covariate_diffusion = ~ depth_scaled,
    mesh = mesh,
    family = gaussian()
  )
  cat("SUCCESS: Basic covariate diffusion model fitted!\n")
  cat("Convergence:", m1_fit$model$convergence, "\n")

  # Check that diffusion parameters exist
  pars <- names(m1_fit$model$par)
  has_diffusion_tau <- "ln_tau_diffusion" %in% pars
  has_diffused_cov <- "diffused_cov_s" %in% pars

  cat("Has ln_tau_diffusion parameter:", has_diffusion_tau, "\n")
  cat("Has diffused_cov_s parameter:", has_diffused_cov, "\n")

  if (has_diffusion_tau) {
    tau_val <- exp(m1_fit$model$par$ln_tau_diffusion)
    cat("tau_diffusion value:", tau_val, "\n")
  }

  # Test prediction
  cat("\n4. Testing prediction with covariate diffusion...\n")
  p1 <- predict(m1_fit)
  cat("Prediction successful, got", nrow(p1), "predictions\n")

  # Check if diffused covariate values are in prediction
  has_diffused_col <- "diffused_depth_scaled" %in% names(p1)
  cat("Has diffused_depth_scaled column:", has_diffused_col, "\n")

  if (has_diffused_col) {
    cat("Original depth_scaled range:", range(d$depth_scaled), "\n")
    cat("Diffused depth_scaled range:", range(p1$diffused_depth_scaled), "\n")

    # Check that diffusion smooths the values (reduces variance)
    original_var <- var(d$depth_scaled)
    diffused_var <- var(p1$diffused_depth_scaled)
    cat("Original variance:", original_var, "Diffused variance:", diffused_var, "\n")
    cat("Diffusion reduces variance:", diffused_var < original_var, "\n")
  }

}, error = function(e) {
  cat("ERROR in basic model fitting:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
})

cat("\n5. Testing error when covariate not in main formula...\n")
tryCatch({
  m_error2 <- sdmTMB(
    data = d,
    formula = density ~ 1,  # depth_scaled not included
    covariate_diffusion = ~ depth_scaled,
    mesh = mesh,
    family = gaussian(),
    do_fit = FALSE
  )
  cat("ERROR: Should have failed but didn't!\n")
}, error = function(e) {
  cat("PASS: Got expected error:", e$message, "\n")
})

cat("\nCovariate diffusion testing completed!\n")
