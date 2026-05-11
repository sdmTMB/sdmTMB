library(sf)
library(dplyr)

ohio_src <- geodaData::ohio_lung

ohio_sf <- ohio_src |>
  dplyr::select(NAME, geometry) |>
  dplyr::rename(county = NAME)
sf::st_crs(ohio_sf) <- sf::st_crs(32617)

# Pull the county-level attributes out of the sf object before reshaping.
ohio_df <- sf::st_drop_geometry(ohio_src)
ohio_df$county <- ohio_df$NAME

# Build a county-year table from the three available lung cancer mortality
# years. The response is the count of deaths, and the offset is the population
# at risk for that same county-year.
years <- c(1968L, 1978L, 1988L)
count_cols <- c("LM68", "LM78", "LM88")
pop_cols <- c("POPM68", "POPM78", "POPM88")
male_pop_cols <- c("POPMB68", "POPMB78", "POPMB88")
female_pop_cols <- c("POPFB68", "POPFB78", "POPFB88")

dat <- do.call(
  rbind,
  lapply(seq_along(years), function(i) {
    pop_male <- ohio_df[[male_pop_cols[[i]]]]
    pop_female <- ohio_df[[female_pop_cols[[i]]]]
    data.frame(
      county = ohio_df$county,
      year = years[[i]],
      cases = ohio_df[[count_cols[[i]]]],
      pop = ohio_df[[pop_cols[[i]]]],
      pct_male = pop_male / (pop_male + pop_female),
      stringsAsFactors = FALSE
    )
  })
)

# Keep the county id and year explicit. `year` stays numeric so it can be used
# as the time column, and `as.factor(year)` in the formula gives a separate
# baseline for each year.
dat$county <- as.character(dat$county)

ohio_df <- dat

usethis::use_data(ohio_df, ohio_sf, overwrite = TRUE)
