# Preferential-sampling likelihood for the RTMB objective: a Bernoulli model
# for which cells of the sampling frame were visited,
#   logit(p) = Z gamma + alpha_t + b h + u_t (h - hbar_t) + xi,
# where `h` is the main model's log expected catch at the frame rows, computed
# from the same latent effects as the catch likelihood, `b` is `b_pref`,
# `u_t` is an optional temporal deviation of the preference coefficient,
# `hbar_t` is the mean of `h` without its random fields over the frame rows of
# time step t, `alpha_t` is
# an optional temporal baseline deviation, and `xi` is the optional
# sampling-only field. The densities of the deviations and `xi` are in
# rtmb_latent_effects().
#
# The deviations multiply the centered target, so they change how sampling
# concentrates within a time step, not (roughly) its rate at the average
# target. An
# uncentered deviation would also shift that rate, which the time step's
# intercept pins down well: the marginal likelihood then penalizes any
# deviation, and the SD of the deviations collapses to 0 even when the
# slopes vary.
# Returns the negative log likelihood and the predictor pieces by frame row.
rtmb_sampling <- function(par, theta, effects, prepared) {
  inputs <- prepared$preferential
  shared <- rtmb_linear_predictors(par, theta, effects, prepared, inputs$rows)
  # Supported families have a log link (a log positive link for deltas, and
  # log links with the offset in component 1 for Poisson-link deltas), so the
  # combined link-scale value is the log expected catch.
  target <- rtmb_combined_link(shared$eta[, 1L], shared$eta[, prepared$n_m],
    prepared$families[[1L]])
  fixed <- rtmb_product(inputs$Z, par$gamma_pref)
  time <- inputs$rows$time
  preference <- par$b_pref * target
  b_t <- NULL
  if (inputs$coefficient != "none") {
    b_t <- par$b_pref + effects$b_pref_t
    # Center on the time step's mean target without the random fields.
    # Averaging the fields over the frame would make every row depend on
    # every mesh vertex of the time step, and the Laplace Hessian dense.
    no_fields <- rtmb_combined_link(shared$fe[, 1L], shared$fe[, prepared$n_m],
      prepared$families[[1L]])
    centered <- target -
      rtmb_product(inputs$time_mean, no_fields)[inputs$time_mean_index]
    preference <- preference + effects$b_pref_t[time] * centered
  }
  baseline <- if (inputs$baseline != "none") {
    effects$alpha_pref_t[time]
  } else {
    rep(0, length(target))
  }
  field <- if (inputs$xi) {
    rtmb_product(inputs$rows$A_rows, effects$xi_s)
  } else {
    rep(0, length(target))
  }
  eta <- fixed + baseline + preference + field
  i <- inputs$observed
  # dbinom_robust() works on the logit scale, so it stays finite for
  # extreme predictors. Unknown (`NA`) indicators contribute nothing.
  nll <- -sum(RTMB::dbinom_robust(inputs$R[i], 1, eta[i], log = TRUE))
  list(nll = nll, target = target, fixed = fixed, baseline = baseline,
    preference = preference, field = field, eta = eta, b_t = b_t,
    alpha_t = if (inputs$baseline != "none") effects$alpha_pref_t)
}

