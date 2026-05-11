#' Example fish survey data
#'
#' @description
#' Various fish survey datasets.
#'
#' @format `pcod`: Trawl survey data for Pacific Cod in Queen Charlotte Sound. A
#'   data frame.
#' @rdname surveydata
"pcod"

#' @format `pcod_2011`: A version of `pcod` for years 2011 and after (smaller
#'   for speed). A data frame.
#' @rdname surveydata
"pcod_2011"

#' @format `pcod_mesh_2011`: A mesh pre-built for `pcod_2011` for examples. A
#'   list of class `sdmTMBmesh`.
#' @rdname surveydata
"pcod_mesh_2011"

#' @format `qcs_grid` A 2x2km prediction grid for Queen Charlotte Sound. A data
#'   frame.
#' @rdname surveydata
"qcs_grid"

#' @format `dogfish`: Trawl survey data for Pacific Spiny Dogfish on West Coast
#'   Vancouver Island. A data frame.
#' @rdname surveydata
"dogfish"

#' @format `yelloweye`: Survey data for Yelloweye Rockfish from the Hard Bottom
#'   Longline Survey (South) off West Coast Vancouver Island.
#' @rdname surveydata
"yelloweye"

#' @format `hbll_s_grid`: A survey domain grid to go with `yelloweye`. A data frame.
#' @rdname surveydata
"hbll_s_grid"

#' @format `wcvi_grid`: A survey domain grid to go with `dogfish`. A data frame.
#' @rdname surveydata
"wcvi_grid"

#' Ohio lung cancer mortality data, reshaped for spatiotemporal examples
#'
#' A derived version of `geodaData::ohio_lung`, reshaped into county-year
#' and county geometry formats for package examples and vignettes.
#'
#' @format `ohio_df`: A county-year data frame with columns:
#'   \describe{
#'   \item{county}{County name.}
#'   \item{year}{Observation year.}
#'   \item{cases}{Observed lung cancer deaths.}
#'   \item{pop}{Population at risk.}
#'   \item{pct_male}{Proportion of the population that is male.}
#'   }
#'
#' @source Derived from `geodaData::ohio_lung`, originally distributed in the
#'   CRAN package `geodaData` under CC0. The original dataset is based on the
#'   GeoDa Center Ohio lung cancer mortality dataset.
#'
#' @examples
#' data(ohio_df)
#' head(ohio_df)
#' data(ohio_sf)
#' head(ohio_sf)
#' @rdname ohio_sf
"ohio_df"

#' @format `ohio_sf`: An `sf` object with columns:
#'   \describe{
#'   \item{county}{County name.}
#'   \item{geometry}{County polygon geometry.}
#'   }
#' @rdname ohio_sf
"ohio_sf"
