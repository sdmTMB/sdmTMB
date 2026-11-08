test_that("make_mesh returns informative error if some coordinates are NA", {
  d <- data.frame(x = c(NA, runif(10)), y = c(NA, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(1, runif(10)), y = c(NA, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(NA, runif(10)), y = c(1, runif(10)))
  expect_error({make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)}, regexp = "NA")

  d <- data.frame(x = c(1, runif(10)), y = c(1, runif(10)))
  mesh <- make_mesh(d, xy_cols = c("x", "y"), cutoff = 1)
})

test_that("make_mesh keeps data locations with cutoff and convex", {
  plain <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
  hull <- make_mesh(pcod, c("X", "Y"), cutoff = 10, convex = 10)
  # at least as many vertices as without the boundary, all data on the mesh
  expect_gte(hull$mesh$n, plain$mesh$n)
  expect_true(all(Matrix::rowSums(hull$A_st) > 0))
})

test_that("make_mesh errors if n_knots is used with type = 'cutoff'", {
  expect_error(
    make_mesh(pcod, c("X", "Y"), n_knots = 50, cutoff = 10, type = "cutoff"),
    regexp = "n_knots"
  )
})

test_that("make_mesh type = 'cutoff_search' finds the requested n_knots", {
  capture.output(
    mesh <- make_mesh(pcod, c("X", "Y"), type = "cutoff_search", n_knots = 50)
  )
  expect_identical(mesh$mesh$n, 50L)
})
