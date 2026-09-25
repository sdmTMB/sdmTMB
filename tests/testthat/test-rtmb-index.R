rtmb_index_fits <- function(backend) {
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 30)
  control <- sdmTMBcontrol(backend = backend)
  list(
    tweedie = sdmTMB(density ~ 0 + as.factor(year), data = pcod_2011,
      mesh = mesh, time = "year", family = tweedie(),
      spatiotemporal = "iid", control = control),
    delta = sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
      time = "year", family = delta_gamma(), spatiotemporal = "off",
      control = control)
  )
}
tmb_index_fits <- fit_once(function() rtmb_index_fits("tmb"))
rtmb_fits <- fit_once(function() rtmb_index_fits("rtmb"))

# Prediction grid without 2013, so one time step is excluded.
rtmb_index_grid <- function() {
  nd <- replicate_df(qcs_grid_small, "year", unique(pcod_2011$year))
  nd[nd$year != 2013, ]
}

test_that("RTMB derived indices match every TMB report and sdreport row", {
  skip_on_cran()
  nd <- rtmb_index_grid()
  for (name in c("tweedie", "delta")) {
    fit <- tmb_index_fits()[[name]]
    data <- predict(fit, newdata = nd, return_tmb_data = TRUE)
    # Unequal areas and all three derived outputs together.
    data$area_i <- rep(c(2, 4), length.out = nrow(nd))
    data$proj_vector <- nd$depth
    data$calc_index_totals <- data$calc_weighted_avg <- data$calc_eao <- 1L
    p <- rtmb_test_parameters(fit$tmb_params, fit$tmb_map)
    p$eps_index <- seq(0.1, 0.3, length.out = data$n_t)
    expect_rtmb_matches_tmb(data, p, fit$tmb_map, fit$tmb_random, info = name)
  }
})

test_that("RTMB get_index(), get_cog(), and get_eao() match TMB", {
  skip_on_cran()
  nd <- rtmb_index_grid()
  cpp <- tmb_index_fits()
  rt <- rtmb_fits()
  # Bias correction is slow for COG and EAO and shares its code path with the
  # index, so it is only checked for the index.
  outputs <- function(fit) suppressMessages(list(
    get_index(fit, newdata = nd, area = 4, bias_correct = FALSE),
    get_cog(fit, newdata = nd, bias_correct = FALSE),
    get_eao(fit, newdata = nd, bias_correct = FALSE)
  ))
  for (name in names(cpp)) {
    expect_equal(outputs(rt[[name]]), outputs(cpp[[name]]), tolerance = 1e-5,
      info = name)
  }
  bc_index <- function(fit) suppressMessages(
    get_index(fit, newdata = nd, area = 4, bias_correct = TRUE))
  expect_equal(bc_index(rt$tweedie), bc_index(cpp$tweedie), tolerance = 1e-5)
})
