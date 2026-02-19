test_that("get_kappa_map() works", {
  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(NA, NA, NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("off", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("off", "on"),
  #   share_range = c(FALSE, FALSE)
  # )
  # expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, 3, 4)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "off"),
    share_range = c(TRUE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, NA, NA)))

  # x <- get_kappa_map(
  #   n_m = 2,
  #   spatial = c("on", "on"),
  #   spatiotemporal = c("on", "on"),
  #   share_range = c(FALSE, TRUE)
  # )
  # expect_identical(x, factor(c(1, 2, 3, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(TRUE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(NA, NA, 1, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "off"),
    spatiotemporal = c("on", "off"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 2, NA, NA)))

  # non-delta:
  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = TRUE
  )
  expect_identical(x, factor(c(NA, NA)))

  #######

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 2)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "on",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "on",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(1, 1)))

  x <- get_kappa_map(
    n_m = 1,
    spatial = "off",
    spatiotemporal = "off",
    share_range = FALSE
  )
  expect_identical(x, factor(c(NA, NA)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("off", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 3)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("off", "off"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, FALSE)
  )
  expect_identical(x, factor(c(1, 1, 2, 2)))

  x <- get_kappa_map(
    n_m = 2,
    spatial = c("on", "on"),
    spatiotemporal = c("on", "on"),
    share_range = c(FALSE, TRUE)
  )
  expect_identical(x, factor(c(1, 2, 3, 3)))
})

test_that("grouped RW parameter mapping coexists with share_range", {
  local_edition(2)

  d <- subset(pcod, year %in% c(2003, 2005, 2007))
  d <- d[seq_len(min(300L, nrow(d))), ]
  d$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 30)

  fit_sr <- sdmTMB(
    density ~ 1,
    data = d,
    time = "year",
    mesh = mesh,
    family = tweedie(link = "log"),
    spatial = "on",
    spatiotemporal = "rw",
    share_range = TRUE,
    do_fit = FALSE,
    experimental = list(groups = "grp")
  )

  fit_no_sr <- update(fit_sr, share_range = FALSE, do_fit = FALSE)

  expect_true(all(!is.na(as.numeric(fit_sr$tmb_map$ln_tau_U))))
  expect_true(all(!is.na(as.numeric(fit_no_sr$tmb_map$ln_tau_U))))
  expect_true(all(is.na(as.numeric(fit_sr$tmb_map$ln_tau_E))))
  expect_true(all(is.na(as.numeric(fit_no_sr$tmb_map$ln_tau_E))))
  expect_true(all(is.na(as.numeric(fit_sr$tmb_map$epsilon_st))))
  expect_true(all(is.na(as.numeric(fit_no_sr$tmb_map$epsilon_st))))
})
