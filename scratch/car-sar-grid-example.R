# CAR/SAR example: dogfish data on a regular areal grid.

pkgload::load_all(".", quiet = TRUE)

library(sf)
library(ggplot2)
theme_set(theme_light())
library(dplyr)

dogfish_points <- sf::st_as_sf(dogfish, coords = c("X", "Y"), crs = NA)

dogfish_boundary <- sf::st_union(dogfish_points) |>
  sf::st_convex_hull() |>
  sf::st_as_sf()

dogfish_grid_obj <- make_areal_grid(
  dogfish,
  xy_cols = c("X", "Y"),
  spatial_domain = dogfish_boundary,
  n = c(25L, 20L),
  space_column = "cell_id"
)
plot(dogfish_grid_obj$grid[[1]])

raw_points <- dogfish_grid_obj$data

# Aggregate the raw tow-level data to one row per cell and year.
cell_data <- raw_points |>
  group_by(cell_id, year) |>
  summarise(
    catch_weight = sum(catch_weight, na.rm = TRUE),
    area_swept = sum(area_swept, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    .groups = "drop"
  )

cell_data$log_depth <- log(cell_data$depth)

fit_dogfish <- sdmTMB(
  catch_weight ~ log_depth,
  data = cell_data,
  mesh = dogfish_grid_obj$domain,
  spatial_model = "car",
  time = "year",
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  offset = log(cell_data$area_swept),
  silent = TRUE
)
sanity(fit_dogfish)
fit_dogfish

# Predict on the full grid by expanding to every cell-year combination.
# For cells with no observations, use the cell-average depth and a common
# effort offset so the map covers the entire spatial domain.
cell_covariates <- raw_points |>
  group_by(cell_id) |>
  summarise(
    depth = mean(depth, na.rm = TRUE),
    .groups = "drop"
  )

pred_data <- tidyr::crossing(
  cell_id = dogfish_grid_obj$grid$cell_id,
  year = sort(unique(raw_points$year))
) |>
  left_join(cell_covariates, by = "cell_id") |>
  mutate(
    depth = ifelse(is.na(depth), mean(raw_points$depth, na.rm = TRUE), depth),
    log_depth = log(depth)
  )

pred_dogfish <- predict(
  fit_dogfish,
  newdata = pred_data,
  type = "response",
  offset = rep(0, nrow(pred_data))
)

pred_dogfish_grid <- left_join(
  dogfish_grid_obj$grid,
  mutate(pred_data,
    est = pred_dogfish$est,
    omega_s = pred_dogfish$omega_s,
    epsilon_st = pred_dogfish$epsilon_st
  )
)

ggplot(pred_dogfish_grid) +
  geom_sf(aes(fill = omega_s), colour = "grey40") +
  scale_fill_gradient2() +
  labs(fill = "Spatial random effects")

ggplot(pred_dogfish_grid) +
  geom_sf(aes(fill = epsilon_st), colour = "grey40") +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  labs(fill = "Spatiotemporal random effects")

ggplot(pred_dogfish_grid) +
  geom_sf(aes(fill = est), colour = "grey40") +
  facet_wrap(~year) +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "Predicted kg catch weight\nper 1km^2 fished")
