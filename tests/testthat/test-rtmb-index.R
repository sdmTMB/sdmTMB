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

test_that("index standard errors reuse the fit's fixed-effect Hessian", {
  skip_on_cran()
  nd <- rtmb_index_grid()
  for (backend in c("tmb", "rtmb")) {
    fit <- if (backend == "tmb") tmb_index_fits()$tweedie else rtmb_fits()$tweedie
    p <- predict(fit, newdata = nd, return_tmb_data = TRUE)
    p$calc_index_totals <- 1L
    pars <- get_pars(fit)
    pars$eps_index <- numeric(0)
    new_obj <- make_sdmTMB_adfun(p, pars, map = fit$tmb_map,
      random = fit$tmb_random, backend = backend)
    H <- fit_hessian_fixed(fit, new_obj, fit$model$par)
    expect_false(is.null(H), info = backend)
    expect_equal(H, solve(fit$sd_report$cov.fixed), info = backend)
  }
})

test_that("index SEs from the fit's joint precision match sdreport()", {
  skip_on_cran()
  nd <- rtmb_index_grid()
  mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 30)
  reml <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh, time = "year",
    spatiotemporal = "off", family = tweedie(), reml = TRUE,
    control = sdmTMBcontrol(backend = "rtmb"))
  fits <- c(rtmb_fits(), list(reml = reml, tmb = tmb_index_fits()$tweedie))
  calls <- list(index = get_index, cog = get_cog)
  for (name in names(fits)) for (f in names(calls)) {
    if (name == "delta" && f == "cog") next
    run <- function() suppressMessages(calls[[f]](fits[[name]], newdata = nd))
    fast <- local({
      local_mocked_bindings(index_sdreport = function(...) stop("fell back"))
      run()
    })
    slow <- local({
      local_mocked_bindings(joint_precision_report = function(...) NULL)
      run()
    })
    expect_equal(fast, slow, tolerance = 1e-6, info = paste(name, f))
  }
})
