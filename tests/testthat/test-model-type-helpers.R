test_that("model type helpers detect regular and multi-family delta states", {
  regular_delta <- list(
    family = structure(
      list(delta = TRUE),
      class = "family"
    ),
    tmb_data = list(multi_family = 0L)
  )
  expect_false(sdmTMB:::.is_multi_family_model(regular_delta))
  expect_false(sdmTMB:::.has_multi_family_delta(regular_delta))
  expect_true(sdmTMB:::.is_delta_like_model(regular_delta))
  expect_true(sdmTMB:::.is_regular_delta_model(regular_delta))

  multi_non_delta <- list(
    family = list(gaussian = gaussian()),
    tmb_data = list(
      multi_family = 1L,
      delta_family_e = c(0L, 0L)
    )
  )
  expect_true(sdmTMB:::.is_multi_family_model(multi_non_delta))
  expect_false(sdmTMB:::.has_multi_family_delta(multi_non_delta))
  expect_false(sdmTMB:::.is_delta_like_model(multi_non_delta))
  expect_false(sdmTMB:::.is_regular_delta_model(multi_non_delta))

  multi_with_delta <- list(
    family = list(delta_gamma = delta_gamma(), poisson = poisson()),
    tmb_data = list(
      multi_family = 1L,
      delta_family_e = c(1L, 0L)
    )
  )
  expect_true(sdmTMB:::.is_multi_family_model(multi_with_delta))
  expect_true(sdmTMB:::.has_multi_family_delta(multi_with_delta))
  expect_true(sdmTMB:::.is_delta_like_model(multi_with_delta))
  expect_false(sdmTMB:::.is_regular_delta_model(multi_with_delta))
})

test_that("model type helpers handle missing fields safely", {
  bare <- list()
  expect_false(sdmTMB:::.is_multi_family_model(bare))
  expect_false(sdmTMB:::.has_multi_family_delta(bare))
  expect_false(sdmTMB:::.is_delta_like_model(bare))
  expect_false(sdmTMB:::.is_regular_delta_model(bare))
})
