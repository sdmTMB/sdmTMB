make_nl_vertex_grid <- function(mesh, years = NULL) {
  loc <- as.data.frame(mesh$mesh$loc[, 1:2, drop = FALSE])
  names(loc) <- c("X", "Y")
  if (is.null(years)) {
    return(loc)
  }
  merge(data.frame(year = years), loc)
}

add_nl_grid_covariates <- function(grid, covariates = c("x1", "x2")) {
  for (i in seq_along(covariates)) {
    covariate <- covariates[[i]]
    if ("year" %in% names(grid)) {
      value <- sin(grid$year / (i + 1)) + grid$X / max(grid$X) +
        i * grid$Y / max(grid$Y)
    } else {
      value <- sin(grid$X / (i + 1)) + i * grid$Y / max(grid$Y)
    }
    grid[[covariate]] <- as.numeric(scale(value))
  }
  grid
}

make_nl_covariate_grid <- function(mesh, years = NULL, covariates = c("x1", "x2")) {
  add_nl_grid_covariates(make_nl_vertex_grid(mesh, years), covariates)
}
