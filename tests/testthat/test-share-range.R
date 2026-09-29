kappa_map <- function(...) get_kappa_map(range_group_labels(...))

test_that("get_kappa_map() works with share_range", {
  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- kappa_map(
  #   n_m = 2,
  #   spatial = c("off", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("off", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 4)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, NA, NA)))

  # x <- kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, TRUE)
  # )
  # expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, 1, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, NA, NA)))

  # non-delta:
  x <- kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(NA, NA)))

  #######

  x <- kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 2)))

  x <- kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(NA, NA)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))
})

test_that("range_groups maps ranges within and across components", {
  on <- c("on", "on")
  st <- c("iid", "iid")
  map <- function(g, spatial = on, spatiotemporal = st) {
    kappa_map(2, spatial, spatiotemporal, c(TRUE, TRUE), g)
  }
  # Defaults match share_range = TRUE
  expect_identical(map(NULL), factor(c(1, 1, 2, 2)))
  # Shared spatial, separate spatiotemporal ranges across components
  expect_identical(
    map(list(c(spatial = "a", spatiotemporal = "b"),
      c(spatial = "a", spatiotemporal = "c"))),
    factor(c(1, 2, 1, 3))
  )
  # All four shared
  expect_identical(
    map(list(c(spatial = "a", spatiotemporal = "a"),
      c(spatial = "a", spatiotemporal = "a"))),
    factor(c(1, 1, 1, 1))
  )
  # Spatiotemporal ranges shared across components; missing spatiotemporal
  # labels follow the spatial label
  expect_identical(
    map(list(c(spatial = "a", spatiotemporal = "st"),
      c(spatial = "b", spatiotemporal = "st"))),
    factor(c(1, 2, 3, 2))
  )
  expect_identical(
    map(list(c(spatial = "a"), c(spatial = "a"))),
    factor(c(1, 1, 1, 1))
  )
  # Labels of fields that are off are ignored
  expect_identical(
    map(list(c(spatial = "a", spatiotemporal = "b"),
      c(spatial = "c", spatiotemporal = "b")), spatiotemporal = c("iid", "off")),
    factor(c(1, 2, 3, 3))
  )
  expect_identical(
    map(list(c(spatial = "a"), c(spatial = "a")), spatial = c("on", "off"),
      spatiotemporal = c("off", "off")),
    factor(c(1, 1, NA, NA))
  )
  expect_error(map(c(spatial = "a")), "one element per model component")
  expect_error(map(list(c(spatial = "a"), c(foo = "a"))), "names")
  expect_error(map(list(c("a", "b"), NULL)), "names")
})

test_that("range_groups builds models that match on both backends", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  build <- function(...) {
    sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
      family = delta_gamma(), do_fit = FALSE, ...)
  }
  expect_error(build(share_range = FALSE,
    range_groups = list(c(spatial = "a"), c(spatial = "b"))),
    "only one of")

  # Separate labels reproduce share_range = FALSE
  f1 <- build(share_range = FALSE)
  f2 <- build(range_groups = list(c(spatial = "a", spatiotemporal = "b"),
    c(spatial = "c", spatiotemporal = "d")))
  expect_identical(f1$tmb_map$ln_kappa, f2$tmb_map$ln_kappa)
  expect_identical(f1$tmb_data$share_range, f2$tmb_data$share_range)

  # Spatial range shared across components, spatiotemporal ranges separate
  fit <- build(range_groups = list(c(spatial = "a", spatiotemporal = "b"),
    c(spatial = "a", spatiotemporal = "c")))
  expect_identical(fit$tmb_map$ln_kappa, factor(c(1, 2, 1, 3)))
  expect_identical(fit$tmb_data$share_range, c(0L, 0L))
  d <- fit$tmb_data
  d$normalize_in_r <- 0L
  p <- fit$tmb_params
  p$ln_tau_O[] <- 0.5
  p$ln_tau_E[] <- 1
  p$ln_kappa[] <- c(-2, -1.5, -2, -1.8)
  args <- list(d, p, fit$tmb_map, random = fit$tmb_random)
  cpp <- do.call(make_sdmTMB_adfun, args)
  rt <- do.call(make_sdmTMB_adfun, c(args, backend = "rtmb"))
  expect_length(cpp$par[names(cpp$par) == "ln_kappa"], 3L)
  expect_equal(rt$fn(rt$par), cpp$fn(cpp$par), tolerance = 1e-7)
  expect_equal(rt$gr(rt$par), cpp$gr(cpp$par), ignore_attr = TRUE,
    tolerance = 1e-6)
})

test_that("PC Matern range priors apply once per range group", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  build <- function(priors, ...) {
    sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
      priors = priors, do_fit = FALSE, ...)
  }
  pc <- pc_matern(range_gt = 10, sigma_lt = 5)

  # The sigma part applies to estimated fields; the range part goes to the
  # first field that is on and has a prior
  range_prior <- function(...) build(...)$tmb_data$range_prior
  sigma_prior <- function(...) build(...)$tmb_data$sigma_prior
  expect_identical(sigma_prior(sdmTMBpriors(matern_s = pc, matern_st = pc)),
    matrix(c(1L, 1L), 2L))
  expect_identical(sigma_prior(sdmTMBpriors(matern_s = pc, matern_st = pc),
    spatial = "off"), matrix(c(0L, 1L), 2L))
  expect_identical(sigma_prior(sdmTMBpriors(matern_st = pc)),
    matrix(c(0L, 1L), 2L))
  # With only spatially varying coefficients, matern_s sets their range prior
  svc <- list(sdmTMBpriors(matern_s = pc), spatial = "off",
    spatiotemporal = "off", spatial_varying = ~ 0 + depth_scaled)
  expect_identical(do.call(sigma_prior, svc), matrix(c(0L, 0L), 2L))
  expect_identical(do.call(range_prior, svc), matrix(c(1L, 0L), 2L))
  expect_identical(range_prior(sdmTMBpriors(matern_s = pc, matern_st = pc)),
    matrix(c(1L, 0L), 2L))
  expect_identical(range_prior(sdmTMBpriors(matern_st = pc)),
    matrix(c(0L, 1L), 2L))
  expect_identical(range_prior(sdmTMBpriors(matern_s = pc, matern_st = pc),
    share_range = FALSE), matrix(c(1L, 1L), 2L))
  expect_identical(range_prior(sdmTMBpriors(matern_s = pc, matern_st = pc),
    spatial = "off"), matrix(c(0L, 1L), 2L))
  expect_identical(range_prior(sdmTMBpriors(matern_s = pc),
    spatiotemporal = "off", family = delta_gamma(),
    range_groups = list(c(spatial = "a"), c(spatial = "a"))),
    matrix(c(1L, 0L, 0L, 0L), 2L))

  # Both backends add the range term once for a range shared across
  # components, including with the Stan Jacobian
  g <- list(c(spatial = "a"), c(spatial = "a"))
  for (bayesian in c(FALSE, TRUE)) {
    nll <- lapply(list(sdmTMBpriors(matern_s = pc), sdmTMBpriors()),
      function(priors) {
        fit <- build(priors, spatiotemporal = "off", family = delta_gamma(),
          range_groups = g, bayesian = bayesian)
        d <- fit$tmb_data
        d$normalize_in_r <- 0L
        p <- fit$tmb_params
        p$ln_tau_O[] <- c(0.5, 0.8)
        p$ln_kappa[] <- -2
        args <- list(d, p, fit$tmb_map, random = fit$tmb_random)
        c(tmb = do.call(make_sdmTMB_adfun, args)$fn(),
          rtmb = do.call(make_sdmTMB_adfun, c(args, backend = "rtmb"))$fn())
      })
    expected <- -rtmb_pc_matern(0.5, -2, pc, stan = bayesian) -
      rtmb_pc_matern(0.8, -2, pc, include_range = FALSE, stan = bayesian)
    expect_equal(nll[[1]] - nll[[2]], c(tmb = expected, rtmb = expected),
      tolerance = 1e-8)
  }
})

test_that("PC Matern priors skip fields that are off", {
  skip_on_cran()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
  pc <- pc_matern(range_gt = 10, sigma_lt = 5)
  # Component 2 has no spatial or spatiotemporal field
  for (bayesian in c(FALSE, TRUE)) {
    nll <- lapply(list(sdmTMBpriors(matern_s = pc, matern_st = pc),
      sdmTMBpriors()), function(priors) {
        fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
          time = "year", family = delta_gamma(), share_range = FALSE,
          spatial = list("on", "off"), spatiotemporal = list("iid", "off"),
          priors = priors, bayesian = bayesian, do_fit = FALSE)
        d <- fit$tmb_data
        d$normalize_in_r <- 0L
        p <- fit$tmb_params
        p$ln_tau_O[] <- 0.5
        p$ln_tau_E[] <- 0.8
        p$ln_kappa[] <- c(-2, -1.5, 0, 0)
        args <- list(d, p, fit$tmb_map, random = fit$tmb_random)
        c(tmb = do.call(make_sdmTMB_adfun, args)$fn(),
          rtmb = do.call(make_sdmTMB_adfun, c(args, backend = "rtmb"))$fn())
      })
    expected <- -rtmb_pc_matern(0.5, -2, pc, stan = bayesian) -
      rtmb_pc_matern(0.8, -1.5, pc, stan = bayesian)
    expect_equal(nll[[1]] - nll[[2]], c(tmb = expected, rtmb = expected),
      tolerance = 1e-8)
  }
})
