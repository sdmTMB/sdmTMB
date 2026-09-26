test_that("check_bounds works", {
  check_bounds(c(1, 2, 3), lower = c(0, 0, 0), upper = c(5, 5, 5))
  expect_warning(check_bounds(c(0, 2, 3), lower = c(0, 0, 0), upper = c(5, 5, 5)))
  p <- c(1, 2)
  names(p) <- c("a", "b")
  expect_warning(check_bounds(p, c(1, 1.9), c(5, 5)), regexp = "lower")
  expect_warning(check_bounds(p, c(0, 0), c(1, 5)), regexp = "upper")
})

test_that("lower and upper work", {
  skip_on_cran()
  d <- subset(pcod, year == 2011)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
  expect_warning({
    suppressMessages({
      m <- sdmTMB(density ~ depth_scaled,
        data = d, mesh = pcod_spde, family = tweedie(link = "log"),
        control = sdmTMBcontrol(
          newton_loops = 0,
          lower = list(ln_phi = 0),
          upper = list(ln_phi = 2.5)))
    })},
    regexp = "upper")
  expect_equal(m$model$par[["ln_phi"]], 2.5, tolerance = 1e-6)
  # FIXME NEWTON LOOPS GOING OUTSIDE BOUNDS?

  expect_warning({
    suppressMessages({
      m <- sdmTMB(density ~ depth_scaled,
        data = d, mesh = pcod_spde, family = tweedie(link = "log"),
        control = sdmTMBcontrol(
          newton_loops = 0,
          lower = list(b_j = c(3, -0.46)),
          upper = list(b_j = c(3.5, -0.45))
        ))
    })}, regexp = "bound")
  expect_equal(m$model$par[[2]], -0.45, tolerance = 1e-6)
})

test_that("covariate diffusion spatial scale has default bounds from the mesh", {
  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  obj <- list(par = c(log_kappaS_nl = 0, log_kappaS_nl = 0, b_j = 1))
  bounds <- .nonlocal_log_kappaS_bounds(mesh$mesh)
  pick <- function(x) unname(x[names(x) == "log_kappaS_nl"])

  loc <- mesh$mesh$loc[, 1:2]
  tv <- mesh$mesh$graph$tv
  edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(3, 1)])
  min_edge <- min(sqrt(rowSums((loc[edges[, 1], ] - loc[edges[, 2], ])^2)))
  diagonal <- sqrt(sum(apply(loc, 2, function(x) diff(range(x)))^2))
  # RMSDK = 2 / kappaS_nl
  expect_equal(2 / exp(bounds), c(10 * diagonal, 0.5 * min_edge))

  lim <- set_limits(obj, lower = list(), upper = list(), mesh = mesh$mesh)
  expect_equal(pick(lim$lower), rep(bounds[[1]], 2L))
  expect_equal(pick(lim$upper), rep(bounds[[2]], 2L))
  expect_equal(unname(lim$lower["b_j"]), -Inf)

  lim <- set_limits(obj, lower = list(log_kappaS_nl = -1), upper = list(),
    mesh = mesh$mesh)
  expect_equal(pick(lim$lower), rep(-1, 2L))
  expect_equal(pick(lim$upper), rep(bounds[[2]], 2L))

  lim <- set_limits(obj, lower = list(), upper = list())
  expect_true(all(is.infinite(pick(lim$lower))))
})
