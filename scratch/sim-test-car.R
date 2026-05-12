sim_car <- function(seed) {
  cat(seed, "\n")
  set.seed(seed)
  N <- 10
  grid <- expand.grid(x = 1:N, y = 1:N)
  n <- nrow(grid)
  distance <- as.matrix(dist(grid))
  W <- array(0, c(n, n))
  W[distance == 1] <- 1
  rownames(W) <- colnames(W) <- as.character(seq_len(n))

  D <- rowSums(W)
  alpha <- 0.8
  tau <- 2
  B0 <- 0.5

  Tau <- tau * (diag(rowSums(W)) - alpha * W)

  re <- mvtnorm::rmvnorm(1, mean = rep(0, nrow(W)), sigma = solve(Tau))

  # spatial data
  spat_data <- data.frame(
    car_region = 1:nrow(W),
    year = 1,
    lon = runif(nrow(W)), # dummy
    lat = runif(nrow(W))
  ) # dummy

  # library(ggplot2)
  # ggplot(spat_data, aes(lon, lat)) + geom_point()

  ln_phi <- log(0.2)
  df <- data.frame(
    car_region = 1:nrow(W),
    x = rnorm(nrow(W), 3, 1),
    resid = rnorm(nrow(W), 0, exp(ln_phi))
  )
  df <- dplyr::left_join(df, spat_data, by = "car_region")

  df$mu <- B0 + re[df$car_region]
  df$y <- B0 + re[df$car_region] + df$resid

  # ggplot(df, aes(lon, lat, colour = mu)) + geom_point() + scale_color_viridis_c()
  # ggplot(df, aes(lon, lat, colour = y)) + geom_point() + scale_color_viridis_c()

  # try with limits, no priors
  m <- sdmTMB(y ~ 1,
    data = df, time = "year",
    mesh = make_areal_domain(df, W, space_column = "car_region"),
    spatial_model = "car",
    spatiotemporal = "off",
    spatial = "on", silent = TRUE
  )

  est_ln_tau <- as.numeric(m$sd_report$par.fixed[which(names(m$sd_report$par.fixed) == "ln_tau_O")])
  est_tau <- exp(est_ln_tau)^2
  # expect_equal(est_tau, tau, tolerance = 0.1)

  td <- tidy(m, "ran_pars")
  alpha_est <- td$estimate[td$term == "alpha_car"]
  # expect_equal(alpha_est, alpha, tolerance = 0.07)

  # names(ln_phi) = "ln_phi"
  ln_phi <- m$sd_report$par.fixed[which(names(m$sd_report$par.fixed) == "ln_phi")]

  data.frame(ln_phi = ln_phi, est_tau = est_tau, alpha_est = alpha_est)
}

# library(future)
# theme_set(ggsidekick::theme_sleek())
# plan(multisession, workers = floor(availableCores() / 2))
# out <- furrr::future_map_dfr(seq_len(50), ~ sim_car(seed = .x), .options = furrr_options(seed = TRUE))

out <- purrr::map_dfr(seq_len(50), ~ sim_car(seed = .x))
par(mfrow = c(1, 3))
hist(exp(out$ln_phi));abline(v = (0.2), col = "red")
hist(log(out$est_tau));abline(v = log(2), col = "red")
hist(out$alpha_est);abline(v = 0.8, col = "red")

mean(out$ln_phi)
log(0.2)

mean(out$est_tau)
2

mean(out$alpha_est)
0.8
