multi_family_anisotropy_fixture <- function() {
  set.seed(21)
  dat <- pcod_2011
  dat$dist <- rep(c("gauss", "delta", "delta"), length.out = nrow(dat))
  gauss_rows <- dat$dist == "gauss"
  delta_rows <- !gauss_rows

  dat$y <- numeric(nrow(dat))
  dat$y[gauss_rows] <- 2 + stats::rnorm(sum(gauss_rows), sd = 0.25)
  present <- stats::rbinom(sum(delta_rows), size = 1, prob = 0.7)
  y_pos <- stats::rgamma(sum(delta_rows), shape = 6, scale = 0.3)
  dat$y[delta_rows] <- ifelse(present == 1, y_pos, 0)

  suppressWarnings(sdmTMB(
    y ~ 1,
    data = dat,
    mesh = pcod_mesh_2011,
    family = list(gauss = gaussian(), delta = delta_gamma()),
    distribution_column = "dist",
    spatial = "on",
    spatiotemporal = "off",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 0L, getsd = TRUE)
  ))
}

test_that("Print anisotropy prints correctly", {
  skip_on_cran()
  skip_on_ci() # slow

  # No anisotropy
  fit1 <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    anisotropy = FALSE
  )

  expect_output(print(fit1), regexp = "range: 33.23")
  # expect_null(plot_anisotropy(fit1))
  # expect_null(plot_anisotropy2(fit1))

  # -------------------
  # Anisotropy with spatial only
  fit_sp_only <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = tweedie(),
    spatial = "on",
    spatiotemporal = "off",
    anisotropy = TRUE
  )

  plot_anisotropy(fit_sp_only)
  plot_anisotropy2(fit_sp_only)
  expect_output(print(fit_sp_only), regexp = "\\(spatial\\): 6.1 to 86.0 at 126")

  # Anisotropy with only spatiotemporal random field
  test_mesh <- make_mesh(data = pcod, xy_cols = c("X", "Y"), cutoff = 20)
  fit_st_only <- sdmTMB(
    data = pcod,
    formula = density ~ 1,
    mesh = test_mesh,
    family = tweedie(),
    spatial = "off",
    spatiotemporal = "iid",
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1)
  )

  expect_output(print(fit_st_only), regexp = "\\(spatiotemporal\\): 16.5 to 29.1 at 54")
# -------------------

  # Anisotropy when not shared across random fields
  fit2 <- sdmTMB(
    data = pcod,
    formula = density ~ 1,
    mesh = test_mesh,
    family = tweedie(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(nlminb_loops = 2, newton_loops = 0)
  )

  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatial\\): ")
  expect_output(cat(print_anisotropy(fit2)), regexp = "\\(spatiotemporal\\): ")

  # Anisotropy when shared across random fields in delta model
  fit_dg_shared <- sdmTMB(
    data = pcod_2011,
    formula = density ~ 1,
    mesh = pcod_mesh_2011,
    family = delta_gamma(),
    share_range = TRUE,
    time = "year",
    anisotropy = TRUE,
    control = sdmTMBcontrol(newton_loops = 1L)
  )

  expect_output(cat(print_anisotropy(fit_dg_shared, m = 1)), regexp = "\\(spatial\\): 36")
  expect_output(cat(print_anisotropy(fit_dg_shared, m = 2)), regexp = "\\(spatial\\): 2")

  # Anisotropy when not shared across random fields in delta model
  fit_dg_not_shared <- sdmTMB(
    data = pcod,
    formula = density ~ 1,
    mesh = test_mesh,
    family = delta_gamma(),
    share_range = FALSE,
    time = "year",
    anisotropy = TRUE
  )

  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatial\\): 19")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 1)), regexp = "\\(spatiotemporal\\): 12")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatial\\): 0")
  expect_output(cat(print_anisotropy(fit_dg_not_shared, m = 2)), regexp = "\\(spatiotemporal\\): 9")
})

test_that("multi-family anisotropy reports both linear predictors", {
  skip_on_cran()
  skip_on_ci() # slow

  fit_multi <- multi_family_anisotropy_fixture()
  aniso_df <- plot_anisotropy(fit_multi, return_data = TRUE)

  expect_equal(sort(unique(aniso_df$model_num)), c(1L, 2L))
  expect_true(all(aniso_df$random_field == "spatial"))
  expect_output(cat(print_anisotropy(fit_multi, m = 1)), regexp = "\\(spatial\\): ")
  expect_output(cat(print_anisotropy(fit_multi, m = 2)), regexp = "\\(spatial\\): ")
})
