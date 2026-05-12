# written to match the .cpp implementation in src/utils.h
.barrier_q_inlaspacetime_formula <- function(fem, ln_kappa, barrier_scaling) {
  range <- sqrt(8) / exp(ln_kappa)
  range_fraction <- if (length(barrier_scaling) > 1) barrier_scaling[2] else 0.1

  CC <- fem$C[[1]] + fem$C[[2]] * range_fraction^2
  Dmat <- fem$D[[1]] + fem$D[[2]] * range_fraction^2
  iC <- Matrix::Diagonal(nrow(fem$I), 1 / CC)

  ICI <- Matrix::t(fem$I) %*% iC %*% fem$I
  ICD <- Matrix::t(fem$I) %*% iC %*% Dmat
  DCI <- Matrix::t(Dmat) %*% iC %*% fem$I
  DCD <- Matrix::t(Dmat) %*% iC %*% Dmat

  pi2s2 <- 2 / pi
  param0 <- pi2s2 / range^2
  param1 <- pi2s2 / 8
  param2 <- param1
  param3 <- range^2 * pi2s2 / 64

  param0 * ICI + param1 * ICD + param2 * DCI + param3 * DCD
}

test_that("Barrier precision matrix matches INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("INLAspacetime")
  skip_if_not_installed("sf")

  pts <- expand.grid(
    x = seq(0, 10, length.out = 20),
    y = seq(0, 6, length.out = 12)
  )
  mesh <- make_mesh(pts, xy_cols = c("x", "y"), cutoff = 0.8)

  barrier_poly <- sf::st_polygon(list(rbind(
    c(4.6, -1.0),
    c(5.4, -1.0),
    c(5.4, 7.0),
    c(4.6, 7.0),
    c(4.6, -1.0)
  )))
  barrier_sf <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(barrier_poly, crs = sf::NA_crs_)
  )

  barrier_triangles <- unlist(
    fmesher::fm_contains(barrier_sf, mesh$mesh, type = "centroid")
  )
  expect_gt(length(barrier_triangles), 0)

  fem <- INLAspacetime::mesh2fem.barrier(
    mesh$mesh,
    barrier.triangles = barrier_triangles
  )

  for (range_fraction in c(0.1, 0.5)) {
    for (range in c(1.5, 2.5, 4.0)) {
      q_expected <- INLA::inla.barrier.q(
        fem,
        ranges = c(range, range * range_fraction),
        sigma = 1
      )
      q_observed <- .barrier_q_inlaspacetime_formula(
        fem,
        ln_kappa = log(sqrt(8) / range),
        barrier_scaling = c(1, range_fraction)
      )

      expect_equal(as.matrix(q_observed), as.matrix(q_expected), tolerance = 1e-12)
    }
  }
})
