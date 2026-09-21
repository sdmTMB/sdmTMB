#' Add UTM coordinates to a data frame
#'
#' Add UTM (Universal Transverse Mercator) coordinates to a data frame. This is
#' useful since geostatistical modeling should generally be performed in a
#' projected coordinate reference system with distance-based units. You can do
#' this yourself with the [sf::st_as_sf()], [sf::st_transform()], and
#' [sf::st_coordinates()] functions in the \pkg{sf} package.
#'
#' @param dat Data frame that contains longitude and latitude columns.
#' @param ll_names Longitude and latitude column names. **Note the order.**
#' @param ll_crs Input CRS value for `ll_names`.
#' @param utm_names Output column names for the UTM columns.
#' @param utm_crs Output CRS value for the UTM zone; tries to detect it with
#'   [get_crs()] but can be specified manually.
#' @param units UTM units.
#'
#' @details
#' **Note that longitudes west of the prime meridian should be encoded
#' as running from -180 to 0 degrees.**
#'
#' You may wish to work in km rather than the standard UTM meters so that the
#' range parameter estimate is not too small, which can cause computational
#' issues. This depends on the scale of your data.
#'
#' @return
#' A copy of the input data frame with new columns for UTM coordinates.
#' @export
#'
#' @examplesIf require("sf", quietly = TRUE)
#' d <- data.frame(lat = c(52.1, 53.4), lon = c(-130.0, -131.4))
#' get_crs(d, c("lon", "lat"))
#' add_utm_columns(d, c("lon", "lat"))
add_utm_columns <- function(dat,
                            ll_names = c("longitude", "latitude"),
                            ll_crs = 4326,
                            utm_names = c("X", "Y"),
                            utm_crs = get_crs(dat, ll_names),
                            units = c("km", "m")) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli_abort("The sf package must be installed to use this function.")
  }

  assert_that(length(ll_names) == 2L)
  assert_that(all(ll_names %in% names(dat)))
  units <- match.arg(units)
  if (any(utm_names %in% names(dat))) {
    cli_abort(c("`utm_names` were found in `names(dat)`.",
      "Remove them or choose different `utm_names`.")
    )
  }
  if (grepl("lat", ll_names[1]) || grepl("lon", ll_names[2])) {
    cli_warn("Make sure you didn't reverse the longitude and latitude in `ll_names`.")
  }

  x <- sf::st_as_sf(dat, crs = ll_crs, coords = ll_names)
  x <- sf::st_transform(x, utm_crs)
  x <- sf::st_coordinates(x)
  x <- as.data.frame(x)
  if (units == "km") {
    x$X <- x$X / 1000
    x$Y <- x$Y / 1000
  }
  dat[[utm_names[1]]] <- x$X
  dat[[utm_names[2]]] <- x$Y
  dat
}

#' @rdname add_utm_columns
#' @export
get_crs <- function(dat, ll_names = c("longitude", "latitude")) {
  lon <- dat[[ll_names[1]]]
  lat <- dat[[ll_names[2]]]
  # https://gis.stackexchange.com/a/190209
  zones <- round((183 + lon) / 6, 0)
  one_zone <- TRUE
  if (length(unique(zones)) > 1L) {
    warning("Multiple UTM zones detected.\n",
      "Proceeding with the most common value.\n",
      "You may wish to choose a different projection.",
      call. = FALSE
    )
    one_zone <- FALSE
  }
  check <- rev(sort(table(zones)))
  zone <- as.numeric(names(check)[[1]])

  lat_zones <- round((45 + lat) / 90, 0)
  if (length(unique(lat_zones)) > 1L) {
    warning("North and south latitudes detected.\n",
      "Proceeding with the most common value.\n",
      "You may wish to choose a different projection.",
      call. = FALSE
    )
    one_zone <- FALSE
  }
  check <- rev(sort(table(lat_zones)))
  lat_zone <- as.numeric(names(check)[[1]])

  crs_val <- 32700 - lat_zone * 100 + zone
  .message <- if (one_zone) "Detected" else "Proceeding with"
  message(
    .message, " UTM zone ", zone, if (lat_zone == 1) "N" else "S", "; CRS = ",
    crs_val, "."
  )
  message("Visit ", paste0("https://epsg.io/", crs_val, " to verify."))

  crs_val
}
