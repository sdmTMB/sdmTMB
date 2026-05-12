## Compare sdmTMB SAR against tinyVAST on the Ohio county mortality data.
##
## This keeps the Ohio areal example structure, but compares the SAR version of
## the model between the two packages. The main check here is whether fitted
## responses agree closely when both models use the same county-year data,
## population offset, and rook-adjacency areal domain.

pkgload::load_all(".", quiet = TRUE)
library(dplyr)
library(ggplot2)
library(sf)
library(tinyVAST)
library(igraph)

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

# tinyVAST uses an igraph spatial domain for areal SAR models. We pass the
# same county neighborhood structure used by sdmTMB.
spatial_graph <- graph_from_adjacency_matrix(
  as.matrix(domain$W_raw),
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)
V(spatial_graph)$name <- domain$unit_names

# tinyVAST needs explicit RAM terms to keep the spatial and spatiotemporal
# latent fields in the model when there is only one response variable.
space_term <- "
  response <-> response, space_sd
"
spacetime_term <- "
  response <-> response, 0, spacetime_sd
"

# Fit the SAR model in tinyVAST.
fit_tinyvast <- tinyVAST(
  formula = cases ~ 0 + as.factor(year) + pct_male,
  data = dat,
  space_term = space_term,
  # spacetime_term = spacetime_term,
  family = poisson(link = "log"),
  spatial_domain = spatial_graph,
  space_columns = "county",
  # time_column = "year",
  control = tinyVASTcontrol(newton_loops = 3L)
)

# Fit the matching SAR model in sdmTMB.
fit_sdmTMB <- sdmTMB(
  cases ~ 0 + as.factor(year) + pct_male,
  data = dat,
  mesh = domain,
  spatial_model = "sar",
  # time = "year",
  family = poisson(link = "log"),
  spatial = "on",
  # spatiotemporal = "iid",
  # offset = log(dat$pop),
  control = sdmTMBcontrol(
    newton_loops = 3L,
    sar_weight_style = "raw"
  )
)

print(fit_tinyvast)
print(fit_sdmTMB)
logLik(fit_sdmTMB)
logLik(fit_tinyvast)

# Compare fixed-effect parameter estimates on the linear predictor scale.
tinyvast_fixed <- summary(fit_tinyvast, what = "fixed") |>
  tibble::rownames_to_column("term") |>
  transmute(term, tinyVAST = Estimate)

sdmTMB_fixed <- tidy(fit_sdmTMB) |>
  transmute(term, sdmTMB = estimate)

fixed_compare <- full_join(tinyvast_fixed, sdmTMB_fixed, by = "term") |>
  mutate(diff = tinyVAST - sdmTMB) |>
  arrange(term)

cat("\nFixed-effect comparison\n")
print(fixed_compare)

# Compare the package-reported latent-field summaries.
# These are not identical internal parameterizations, so read them as model-
# specific SAR summaries rather than exact one-to-one transformed parameters.
tinyvast_latent <- tibble::tibble(
  term = c("sigma_O", "sigma_E"),
  tinyVAST_term = c("space_term", "spacetime_term"),
  tinyVAST = c(
    summary(fit_tinyvast, what = "space_term")$Estimate[1],
    summary(fit_tinyvast, what = "spacetime_term")$Estimate[1]
  )
)

sdmTMB_latent <- tidy(fit_sdmTMB, "ran_pars") |>
  filter(term %in% c("sigma_O", "sigma_E")) |>
  select(term, sdmTMB = estimate, std.error, conf.low, conf.high)

latent_compare <- left_join(tinyvast_latent, sdmTMB_latent, by = "term") |>
  arrange(term)

rho_compare <- tidy(fit_sdmTMB, "ran_pars") |>
  filter(term == "rho_sar") |>
  select(term, estimate, std.error, conf.low, conf.high)

cat("\nLatent-field summaries\n")
print(latent_compare)
print(rho_compare)

# Compare fitted responses on the observed data.
pred_tinyvast <- predict(fit_tinyvast, newdata = dat)
pred_sdmTMB_df <- predict(
  fit_sdmTMB,
  newdata = dat,
  type = "response"
  # offset = log(dat$pop)
)

pred_compare <- dat |>
  mutate(
    pred_tinyVAST = pred_tinyvast,
    pred_sdmTMB = pred_sdmTMB_df$est,
    diff = pred_tinyVAST - pred_sdmTMB,
    abs_diff = abs(diff),
    rate_obs = cases / pop,
    rate_tinyVAST = pred_tinyVAST / pop,
    rate_sdmTMB = pred_sdmTMB / pop
  )

cat("\nPrediction comparison\n")
cat(sprintf("Count correlation: %.6f\n", cor(pred_compare$pred_tinyVAST, pred_compare$pred_sdmTMB)))
cat(sprintf("Count RMSE: %.6f\n", sqrt(mean(pred_compare$diff^2))))
cat(sprintf("Count max abs diff: %.6f\n", max(pred_compare$abs_diff)))
cat(sprintf("Rate correlation: %.6f\n", cor(pred_compare$rate_tinyVAST, pred_compare$rate_sdmTMB)))

print(
  pred_compare |>
    select(county, year, cases, pop, pred_tinyVAST, pred_sdmTMB, diff) |>
    arrange(desc(abs(diff))) |>
    head()
)

pred_map <- left_join(
  ohio_sf,
  pred_compare,
  by = "county"
)

ggplot(pred_compare, aes(pred_tinyVAST, pred_sdmTMB)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  coord_equal() +
  labs(
    title = "tinyVAST SAR vs sdmTMB SAR",
    x = "tinyVAST predicted cases",
    y = "sdmTMB predicted cases"
  )

if (requireNamespace("tidyr", quietly = TRUE)) {
  pred_long <- tidyr::pivot_longer(
    pred_map,
    cols = c(rate_tinyVAST, rate_sdmTMB),
    names_to = "model",
    values_to = "rate"
  )

  ggplot(pred_long) +
    geom_sf(aes(fill = rate), colour = "grey40", linewidth = 0.15) +
    facet_wrap(~model) +
    scale_fill_viridis_c() +
    labs(title = "Predicted county mortality rate", fill = "rate")
}
