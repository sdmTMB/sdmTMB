test_that("family specification validator enforces component invariants", {
  spec <- .compile_family_spec(
    list(count = poisson(), biomass = delta_gamma(type = "poisson-link")),
    data = data.frame(dist = c("count", "biomass")),
    distribution_column = "dist"
  )
  expect_invisible(.validate_family_spec(spec))

  invalid <- spec
  invalid$components <- invalid$components[invalid$components$component == 1L, ]
  expect_error(
    .validate_family_spec(invalid),
    "combine kind does not match active components"
  )

  invalid <- spec
  invalid$family_id_i[1L] <- NA_integer_
  expect_error(.validate_family_spec(invalid), "every row must have one valid family ID")
})

test_that("response preparation is row-family aware", {
  dat <- data.frame(
    y = c(0.2, 1, 2, 0),
    dist = c("encounter", "encounter", "count", "count")
  )
  spec <- .compile_family_spec(
    list(encounter = binomial(), count = poisson()),
    data = dat,
    distribution_column = "dist"
  )
  response <- .prepare_family_response(dat$y, c(10, NA, 1, 1), spec)

  expect_equal(response$y_i, c(2, 1, 2, 0))
  expect_equal(response$size, c(10, 1, 1, 1))
  expect_equal(response$weights, rep(1, 4))
  expect_equal(response$response[, 1L], response$y_i)
})

test_that("ordinary binomial response forms use the shared processor", {
  spec <- .compile_family_spec(binomial(), data = data.frame(y = 1:2))

  factor_response <- .prepare_family_response(factor(c("no", "yes")), NULL, spec)
  expect_equal(factor_response$y_i, c(0, 1))

  matrix_response <- .prepare_family_response(
    cbind(success = c(2, 3), failure = c(8, 7)), NULL, spec
  )
  expect_equal(matrix_response$y_i, c(2, 3))
  expect_equal(matrix_response$size, c(10, 10))
})

test_that("family TMB boundary owns auxiliary parameter sizing and mapping", {
  spec <- .compile_family_spec(
    list(count = nbinom2(), biomass = delta_gamma(type = "poisson-link")),
    data = data.frame(dist = c("count", "biomass")),
    distribution_column = "dist"
  )
  params <- .family_parameter_values(spec, FALSE, NULL)
  expect_length(params$ln_phi, 2L)
  expect_length(params$thetaf, 0L)

  tmb_params <- c(params, list(other = 0))
  mapped <- .map_family_parameters(
    map_all_params(tmb_params), tmb_params,
    has_dispformula = FALSE, estimate_student_df = FALSE
  )
  expect_null(mapped$ln_phi)
  expect_true(all(c("component_active", "family_code", "link_code") %in%
    names(.as_tmb_family_data(spec))))
})
