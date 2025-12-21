test_that("Multi-likelihood family list validation works", {
  fam <- list(
    gauss = gaussian(),
    pois = poisson(link = "log")
  )
  res <- sdmTMB:::.validate_multi_family_list(fam)
  expect_identical(res$family_names, c("gauss", "pois"))
  expect_equal(res$family_enum[[1]], sdmTMB:::.enum_family("gaussian"))
  expect_equal(res$link_enum[[2]], sdmTMB:::.enum_link("log"))
})

test_that("Multi-likelihood family list validation rejects mix families", {
  fam <- list(mix = gamma_mix())
  expect_error(
    sdmTMB:::.validate_multi_family_list(fam),
    regexp = "_mix"
  )
})

test_that("Multi-likelihood family list validation requires names", {
  fam <- list(gaussian(), poisson())
  expect_error(
    sdmTMB:::.validate_multi_family_list(fam),
    regexp = "named list"
  )
})

test_that("Multi-likelihood distribution_column maps to e_i", {
  fam <- list(
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  dat <- data.frame(
    y = c(1.2, 3, 0, 2.4, 1, 0),
    dist = c("gaussian", "poisson", "binomial", "gaussian", "poisson", "binomial")
  )
  res <- sdmTMB:::.validate_multi_family_list(fam, data = dat, distribution_column = "dist")
  expect_equal(res$e_i, c(1L, 2L, 3L, 1L, 2L, 3L))
})
