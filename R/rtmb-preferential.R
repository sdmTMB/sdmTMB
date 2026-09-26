# Preferential-sampling likelihood for the RTMB objective: a Bernoulli model
# for which cells of the sampling frame were visited,
#   logit(p) = Z gamma + b h + xi,
# where `h` is the main model's log expected catch at the frame rows, computed
# from the same latent effects as the catch likelihood, and `xi` is the
# optional sampling-only field (its density is in rtmb_latent_effects()).
# Returns the negative log likelihood and the predictor pieces by frame row.
rtmb_sampling <- function(par, theta, effects, prepared) {
  inputs <- prepared$preferential
  shared <- rtmb_linear_predictors(par, theta, effects, prepared, inputs$rows)
  target <- rtmb_log_mean(shared$eta, prepared$families[[1L]])
  fixed <- rtmb_product(inputs$Z, par$gamma_pref)
  preference <- par$b_pref * target
  field <- if (inputs$xi) {
    rtmb_product(inputs$rows$A_rows, effects$xi_s)
  } else {
    rep(0, length(target))
  }
  eta <- fixed + preference + field
  i <- inputs$observed
  # dbinom_robust() works on the logit scale, so it stays finite for
  # extreme predictors. Unknown (`NA`) indicators contribute nothing.
  nll <- -sum(RTMB::dbinom_robust(inputs$R[i], 1, eta[i], log = TRUE))
  list(nll = nll, target = target, fixed = fixed, preference = preference,
    field = field, eta = eta)
}

# Log expected response from link-scale predictors (including any offset) for
# the families preferential sampling supports. For a log-link family this is
# the predictor itself. For a conventional logit/log delta family it is
# log(plogis(eta1)) + eta2, the log of encounter probability times positive
# mean, with log(plogis(x)) = -log(1 + exp(-x)) computed stably.
rtmb_log_mean <- function(eta, family) {
  switch(family$combine,
    single = eta[, 1L],
    delta = -RTMB::logspace_add(0, -eta[, 1L]) + eta[, 2L],
    cli_abort("Internal error: unsupported family for the sampling target."))
}
