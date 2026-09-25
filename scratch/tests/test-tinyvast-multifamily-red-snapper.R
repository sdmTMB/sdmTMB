if (!requireNamespace("tinyVAST", quietly = TRUE)) {
  testthat::skip("tinyVAST not installed")
}

library(testthat)
library(sdmTMB)
library(tinyVAST)
library(fmesher)

expect_same_loglik <- function(fit_sd, fit_tv, tolerance = 1e-5) {
  expect_equal(
    as.numeric(logLik(fit_sd)),
    as.numeric(logLik(fit_tv)),
    tolerance = tolerance
  )
}

test_that("tinyVAST/sdmTMB red snapper multi-family regression check", {
  skip_on_cran()

  data(red_snapper)
  data(red_snapper_shapefile)

  family_list <- list(
    Encounter = binomial(link = "cloglog"),
    Count = poisson(link = "log"),
    Biomass_KG = tweedie(link = "log")
  )

  red_snapper$Data_type <- relevel(red_snapper$Data_type, ref = "Biomass_KG")
  red_snapper$var <- "logdens"

  form_sd <- Response_variable ~ Data_type + factor(Year)
  form_tv <- Response_variable ~ Data_type + factor(Year) + offset(log(AreaSwept_km2))

  mesh_tv <- fmesher::fm_mesh_2d(red_snapper[, c("Lon", "Lat")], cutoff = 1.2)
  mesh_sd <- make_mesh(red_snapper, c("Lon", "Lat"), mesh = mesh_tv)

  fit_sd <- sdmTMB(
    form_sd,
    data = red_snapper,
    mesh = mesh_sd,
    time = "Year",
    spatial = "on",
    spatiotemporal = "off",
    family = family_list,
    distribution_column = "Data_type",
    offset = log(red_snapper$AreaSwept_km2),
    weights = rep(1, nrow(red_snapper)),
    control = sdmTMBcontrol(multiphase = FALSE, newton_loops = 1L)
  )

  fit_tv <- tinyVAST(
    data = red_snapper,
    formula = form_tv,
    space_term = "logdens <-> logdens, sd_space",
    # spacetime_term = "logdens <-> logdens, 0, sd_spacetime",
    spacetime_term = NULL,
    space_columns = c("Lon", "Lat"),
    spatial_domain = mesh_tv,
    time_column = "Year",
    distribution_column = "Data_type",
    family = family_list,
    variable_column = "var",
    delta_options = list(formula = ~ Data_type + factor(Year) + offset(log(AreaSwept_km2))),
    control = tinyVASTcontrol(newton_loops = 1L)
  )

  expect_same_loglik(fit_sd, fit_tv, tolerance = 1e-4)

  sf_grid <- sf::st_make_grid(red_snapper_shapefile, cellsize = c(2, 2))
  sf_grid <- suppressWarnings(sf::st_intersection(sf_grid, red_snapper_shapefile))
  sf_grid <- sf::st_make_valid(sf_grid)

  grid_coords <- sf::st_coordinates(sf::st_centroid(sf_grid))
  areas_km2 <- as.numeric(sf::st_area(sf_grid)) / 1e6

  years <- sort(unique(red_snapper$Year))
  data_types <- levels(red_snapper$Data_type)

  for (dt in data_types) {
    est_sd <- numeric(length(years))
    se_sd <- numeric(length(years))
    est_tv <- numeric(length(years))
    se_tv <- numeric(length(years))

    for (ii in seq_along(years)) {
      yr <- years[ii]
      nd <- data.frame(
        Lat = grid_coords[, "Y"],
        Lon = grid_coords[, "X"],
        Year = yr,
        Data_type = dt,
        AreaSwept_km2 = mean(red_snapper$AreaSwept_km2),
        var = "logdens"
      )

      pred_sd <- predict(
        fit_sd,
        newdata = nd,
        return_tmb_object = TRUE,
        offset = rep(log(mean(red_snapper$AreaSwept_km2)), nrow(nd))
      )
      idx_sd <- get_index(pred_sd, area = areas_km2, bias_correct = TRUE)
      idx_tv <- integrate_output(fit_tv, newdata = nd, area = areas_km2, apply.epsilon = TRUE)

      est_sd[ii] <- idx_sd$est
      se_sd[ii] <- idx_sd$se_natural
      est_tv[ii] <- unname(idx_tv[["Est. (bias.correct)"]])
      se_tv[ii] <- unname(idx_tv[["Std. Error"]])
    }

    expect_equal(est_sd, est_tv, tolerance = 1e-3)
    good <- is.finite(se_sd) & is.finite(se_tv)
    if (any(good)) {
      expect_equal(se_sd[good], se_tv[good], tolerance = 5e-3)
    }
  }
})
