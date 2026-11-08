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

test_that("make_mesh kmeans caps n_knots at the number of unique locations", {
  d <- data.frame(x = rep(1:10, 10), y = rep(c(1, 3, 2, 5, 4, 7, 6, 9, 8, 10), 10))
  expect_warning(
    mesh <- make_mesh(d, c("x", "y"), n_knots = 50),
    regexp = "unique locations"
  )
  expect_identical(nrow(mesh$loc_centers), 9L)
})

test_that("make_mesh does not change the global random seed", {
  set.seed(1)
  a <- runif(1)
  set.seed(1)
  make_mesh(pcod, c("X", "Y"), n_knots = 50)
  expect_identical(runif(1), a)
})

test_that("make_mesh warns about large negative coordinates", {
  d <- data.frame(x = -seq(5e5, 6e5, length.out = 20), y = seq(1, 20))
  expect_warning(make_mesh(d, c("x", "y"), cutoff = 1e4), regexp = "fairly large")
})

test_that("make_mesh checks xy_cols", {
  expect_error(make_mesh(pcod, c("X", "Yy"), cutoff = 10), regexp = "Yy")
  expect_error(make_mesh(pcod, "X", cutoff = 10), regexp = "two column names")
  d <- pcod
  d$Y <- as.character(d$Y)
  expect_error(make_mesh(d, c("X", "Y"), cutoff = 10), regexp = "must be numeric")
})
