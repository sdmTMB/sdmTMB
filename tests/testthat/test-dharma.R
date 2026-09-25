test_that("simulate() and dharma_residuals() work", {
  skip_on_cran()

  mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 15)
  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = tweedie()
  )
  s <- simulate(fit, nsim = 100)
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  expect_warning(dharma_residuals(s, fit), regexp = "mle")

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    family = delta_gamma()
  )
  s <- simulate(fit, nsim = 100, type = 'mle-mvn')
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  dharma_residuals(s, fit)

  fit <- sdmTMB(density ~ 1,
    data = pcod, mesh = mesh,
    spatial = "off",
    family = delta_gamma(type = "poisson-link")
  )
  s <- simulate(fit, nsim = 100, type = "mle-mvn")
  expect_equal(ncol(s), 100)
  expect_equal(nrow(s), nrow(pcod))
  dharma_residuals(s, fit)
  dharma_residuals(s, fit, plot = FALSE)

  s <- simulate(fit, nsim = 100, params = "mle-mvn")
  expect_equal(ncol(s), 100)
  s <- simulate(fit, nsim = 100, params = "mle-mvn", re_form = NA)
  expect_equal(ncol(s), 100)

  # s <- simulate(fit, nsim = 100, type = 'mle-mvn')
  # dharma_residuals(s, fit, expected_distribution = "normal")
  # x <- dharma_residuals(s, fit, return_DHARMa = TRUE)
  # expect_s3_class(x, "DHARMa")
})

test_that("dharma_residuals() uses family_spec for multi-family observed responses", {
  skip_on_cran()
  skip_if_not_installed("DHARMa")

  make_fit <- function(family_list) {
    set.seed(401)
    x <- seq(-1, 1, length.out = 90)
    dist <- rep(c("gauss", "delta", "delta"), length.out = length(x))
    gauss_rows <- dist == "gauss"
    delta_rows <- !gauss_rows

    y <- numeric(length(x))
    y[gauss_rows] <- 1 + 0.4 * x[gauss_rows] + stats::rnorm(sum(gauss_rows), sd = 0.15)
    present <- stats::rbinom(sum(delta_rows), size = 1, prob = stats::plogis(-0.3 + 0.7 * x[delta_rows]))
    y_pos <- stats::rgamma(sum(delta_rows), shape = 6, scale = exp(0.2 + 0.3 * x[delta_rows]) / 6)
    y[delta_rows] <- ifelse(present == 1, y_pos, 0)

    sdmTMB(
      y ~ x,
      data = data.frame(y = y, x = x, dist = dist),
      spatial = "off",
      spatiotemporal = "off",
      family = family_list,
      distribution_column = "dist",
      control = sdmTMBcontrol(newton_loops = 0, getsd = FALSE)
    )
  }

  for (family_list in list(
    list(gauss = gaussian(), delta = delta_gamma()),
    list(delta = delta_gamma(), gauss = gaussian())
  )) {
    fit <- make_fit(family_list)
    sims <- simulate(fit, nsim = 20, type = "mle-mvn")
    dharma_obj <- dharma_residuals(sims, fit, return_DHARMa = TRUE, plot = FALSE)
    expected_y <- ifelse(!is.na(fit$response[, 2]), fit$response[, 2], fit$response[, 1])

    expect_equal(dharma_obj$observedResponse, expected_y)
  }
})
