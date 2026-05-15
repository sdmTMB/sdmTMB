test_that("make_example returns expected structure for all diffusion types", {
  types <- c("space diffusion", "time diffusion", "space-time diffusion", "spacetime diffusion")

  for (tp in types) {
    ex <- make_example(type = tp, seed = 1)

    expect_true(is.list(ex))
    expect_true(all(c("data", "mesh", "covariate_diffusion", "time", "truth") %in% names(ex)))
    expect_s3_class(ex$mesh, "sdmTMBmesh")

    expect_true(is.data.frame(ex$data))
    expect_true(all(c("X", "Y", "year", "x1", "x1_truth", "observed") %in% names(ex$data)))
    expect_equal(nrow(ex$data), 1200)

    if (tp == "space diffusion") {
      expect_null(ex$time)
      expect_equal(deparse(ex$covariate_diffusion), "~space(x1)")
    } else if (tp == "time diffusion") {
      expect_equal(ex$time, "year")
      expect_equal(deparse(ex$covariate_diffusion), "~time(x1)")
    } else if (tp == "space-time diffusion") {
      expect_equal(ex$time, "year")
      expect_equal(deparse(ex$covariate_diffusion), "~space(x1) + time(x1)")
    } else {
      expect_equal(ex$time, "year")
      expect_equal(deparse(ex$covariate_diffusion), "~spacetime(x1)")
    }
  }
})
