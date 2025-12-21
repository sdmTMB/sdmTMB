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
