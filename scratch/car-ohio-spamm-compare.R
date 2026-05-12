## Compare sdmTMB CAR against spaMM on the Ohio county mortality data.
##
## This keeps the Ohio example structure, but switches the areal model to CAR
## and compares the resulting fixed effects and fitted values against a spaMM
## model with the same adjacency matrix.
##
## Note: spaMM does not have the same spatiotemporal IID term used in the
## original Ohio example, so this comparison uses a purely spatial areal model
## with year as a fixed effect.

pkgload::load_all(".", quiet = TRUE)

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop("Package `spaMM` is required for this comparison script.", call. = FALSE)
}

library(dplyr)
library(ggplot2)
library(sf)
theme_set(theme_light())

# Package data: county-year observations plus county polygons.
dat <- ohio_df

# Build the areal domain from the county polygons.
# Use rook adjacency so the spatial graph matches shared-border contiguity.
domain <- make_areal_domain(
  data = dat,
  spatial_domain = ohio_sf,
  space_column = "county",
  id_column = "county",
  adjacency = "rook"
)

# Align county levels with the adjacency matrix for spaMM.
dat_spamm <- dat
dat_spamm$county <- factor(dat_spamm$county, levels = domain$unit_names)

# Fit the CAR model in sdmTMB.
fit_car <- sdmTMB(
  cases ~ 0 + as.factor(year) + pct_male,
  data = dat,
  mesh = domain,
  spatial_model = "car",
  family = poisson(link = "log"),
  spatial = "on",
  offset = log(dat$pop),
  silent = TRUE
)

# Fit the equivalent spaMM model using the same adjacency matrix.
fit_spamm <- spaMM::fitme(
  cases ~ 0 + as.factor(year) + pct_male + offset(log(pop)) + adjacency(1 | county),
  data = dat_spamm,
  family = poisson(),
  adjMatrix = domain$W_raw
)

print(summary(fit_spamm))
print(fit_car)

# Compare fixed effects on the linear predictor scale.
car_fixef <- coef(fit_car)
spamm_fixef <- spaMM::fixef(fit_spamm)

common_terms <- intersect(names(car_fixef), names(spamm_fixef))
coef_compare <- data.frame(
  term = common_terms,
  sdmTMB = unname(car_fixef[common_terms]),
  spaMM = unname(spamm_fixef[common_terms])
) |>
  mutate(
    diff = sdmTMB - spaMM,
    abs_diff = abs(diff)
  ) |>
  arrange(term)

print(coef_compare)

# Compare fitted values on the response scale.
pred_car <- predict(
  fit_car,
  newdata = dat,
  type = "response",
  offset = log(dat$pop)
)

pred_spamm <- predict(
  fit_spamm,
  newdata = dat_spamm,
  type = "response"
)

pred_compare <- dat |>
  mutate(
    county = as.character(county),
    observed = cases,
    pred_sdmTMB = pred_car$est,
    pred_spaMM = as.numeric(pred_spamm),
    rate_obs = cases / pop,
    rate_sdmTMB = pred_car$est / pop,
    rate_spaMM = as.numeric(pred_spamm) / pop
  )

cat("\nPrediction comparison\n")
cat(sprintf("Count correlation: %.4f\n", cor(pred_compare$pred_sdmTMB, pred_compare$pred_spaMM)))
cat(sprintf("Count RMSE: %.4f\n", sqrt(mean((pred_compare$pred_sdmTMB - pred_compare$pred_spaMM)^2))))
cat(sprintf("Rate correlation: %.4f\n", cor(pred_compare$rate_sdmTMB, pred_compare$rate_spaMM)))

pred_map <- left_join(
  ohio_sf,
  pred_compare,
  by = "county"
)

if (requireNamespace("tidyr", quietly = TRUE)) {
  pred_long <- tidyr::pivot_longer(
    pred_map,
    cols = c(rate_sdmTMB, rate_spaMM),
    names_to = "model",
    values_to = "rate"
  )

  ggplot(pred_long) +
    geom_sf(aes(fill = rate), colour = "grey40", linewidth = 0.15) +
    facet_wrap(~model) +
    scale_fill_viridis_c() +
    labs(title = "Predicted county mortality rate", fill = "rate")
}

ggplot(pred_compare, aes(pred_sdmTMB, pred_spaMM)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  coord_equal() +
  labs(
    title = "sdmTMB CAR vs spaMM",
    x = "sdmTMB predicted cases",
    y = "spaMM predicted cases"
  )

