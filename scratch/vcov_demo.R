library(sdmTMB)

mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
fit <- sdmTMB(
  density ~ log(depth) + I(log(depth)^2),
  data = pcod_2011,
  mesh = mesh,
  family = tweedie()
)

vcov_mat <- vcov(fit)
vcov_mat

coefs <- coef(fit)
coefs

sim_coefs <- MASS::mvrnorm(n = 200, mu = coefs, Sigma = vcov_mat)

hist(sim_coefs[,"log(depth)"])
plot(sim_coefs[, "(Intercept)"], sim_coefs[, "log(depth)"])
plot(sim_coefs[, "I(log(depth)^2)"], sim_coefs[, "log(depth)"])
