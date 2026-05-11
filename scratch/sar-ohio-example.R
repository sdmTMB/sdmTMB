## SAR example: Ohio county lung cancer mortality.
##
## This is an areal, spatiotemporal Poisson example:
## - counties are already stored as polygons in `geodaData::ohio_lung`
## - the year-specific counts are reshaped to long format
## - the model uses a population offset so it estimates mortality counts/rates

pkgload::load_all(".", quiet = TRUE)
library(sf)
library(ggplot2)
theme_set(theme_light())
library(dplyr)

# package data; modified from the geodaData R package:
dat <- ohio_df

# Build the areal domain directly from the county polygons.
domain <- make_areal_domain(
  data = dat,
  spatial_domain = ohio_sf,
  space_column = "county",
  id_column = "county"
)

# Fit a Poisson GLMM areal SAR model.
# The offset turns counts into mortality rates per person at risk, while the
# spatial and spatiotemporal random effects capture county clustering and
# county-by-year deviations.

fit <- sdmTMB(
  cases ~ 0 + as.factor(year) + pct_male,
  data = dat,
  mesh = domain,
  time = "year",
  family = poisson(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  offset = log(dat$pop)
)
sanity(fit)
fit

pred_full <- predict(fit, newdata = dat, type = "response", offset = log(dat$pop))

pred_dat <- data.frame(
  county = dat$county,
  year = dat$year,
  cases = dat$cases,
  pop = dat$pop,
  rate_obs = 1e5 * dat$cases / dat$pop,
  rate_fit_full = 1e5 * pred_full$est / dat$pop,
  omega_s = pred_full$omega_s,
  epsilon_st = pred_full$epsilon_st
)

pred_map <- dplyr::left_join(ohio_sf, pred_dat, by = "county")

ggplot(pred_map) +
  geom_sf(aes(fill = omega_s), colour = "grey40") +
  facet_wrap(~year) +
  labs(title = "Spatial random effects") +
  scale_fill_gradient2()

ggplot(pred_map) +
  geom_sf(aes(fill = epsilon_st), colour = "grey40") +
  facet_wrap(~year) +
  labs(title = "Spatiotemporal random effects") +
  scale_fill_gradient2()

ggplot(pred_map) +
  geom_sf(aes(fill = rate_fit_full), colour = "grey40") +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  labs(title = "Population-level fitted mortality rate per 100,000")

# Could we have done that with the SPDE?

# Yes, we could approximate county-level risk by evaluating a continuous latent field at each county centroid.

# That approximation is reasonable when:
#
# - polygons are small;
# - polygons are similarly sized;
# - the latent field is smooth over each polygon;
# - centroids are inside or near the effective population center;
# - the response is genuinely associated with a point-like location.
#
# It is weaker when:
#
# - polygons are large;
# - polygons vary a lot in size;
# - polygons are long, coastal, concave, or irregular;
# - populations are unevenly distributed within polygons;
# - the response is an aggregate count, rate, or average over the whole polygon.

# To use the SPDE, we collapse each county to a centroid and
# treat the observations as point-referenced data rather
# than areal units.

ohio_centroids <- suppressWarnings(sf::st_centroid(ohio_sf))
ohio_xy <- sf::st_coordinates(ohio_centroids)
ohio_spde_dat <- dplyr::left_join(
  ohio_df,
  data.frame(
    county = ohio_sf$county,
    X = ohio_xy[, "X"] / 1000,
    Y = ohio_xy[, "Y"] / 1000
  ),
  by = "county"
)

ohio_spde_mesh <- make_mesh(ohio_spde_dat, c("X", "Y"), cutoff = 20)
plot(ohio_spde_mesh)

fit_spde <- sdmTMB(
  cases ~ 0 + as.factor(year) + pct_male,
  data = ohio_spde_dat,
  mesh = ohio_spde_mesh,
  family = poisson(link = "log"),
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  offset = log(ohio_spde_dat$pop)
)
sanity(fit_spde)
fit_spde

pred_spde <- predict(
  fit_spde,
  newdata = ohio_spde_dat,
  type = "response",
  offset = log(ohio_spde_dat$pop)
)

pred_spde_dat <- data.frame(
  county = ohio_spde_dat$county,
  year = ohio_spde_dat$year,
  rate_fit_spde = 1e5 * pred_spde$est / ohio_spde_dat$pop
)

pred_map_spde <- dplyr::left_join(ohio_sf, pred_spde_dat, by = "county")

ggplot(pred_map_spde) +
  geom_sf(aes(fill = rate_fit_spde), colour = "grey40") +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  labs(title = "SPDE fitted mortality rate per 100,000")

combined <- pred_map_spde
combined$rate_fit_sar <- pred_map$rate_fit_full
ggplot(combined, aes(rate_fit_sar, rate_fit_spde)) + geom_point() +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1)

tidy(fit)
tidy(fit_spde)
