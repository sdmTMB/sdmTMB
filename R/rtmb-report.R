# Register reports with the C++ template's names and dimensions. Post-fit
# methods index these names.
rtmb_report <- function(par, theta, prepared, effects, fitted, obs,
                        projected, derived, simulating) {
  n_m <- prepared$n_m
  any_field <- prepared$any_field
  spde_field <- any_field && !rtmb_areal(prepared$precision)
  areal_type <- prepared$precision$type
  in_simulation <- length(simulating) > 0L || inherits(obs$y_i, "simref")
  requested <- prepared$derived
  y_i <- obs$y_i
  # Report drawn responses, not the simulation reference.
  if (inherits(y_i, "simref")) y_i <- y_i$value
  report_y <- inherits(obs$y_i, "simref") || !prepared$simulate_obs
  diffusion_scales <- if (prepared$diffusion$n_terms > 0L) {
    rtmb_diffusion_scales(par, prepared$diffusion)
  }

  reports <- c(
    # Parameters
    if (prepared$dispersion_model) {
      list(ln_phi_i = obs$ln_phi_i, phi_i = obs$phi_i)
    } else if (length(theta$phi)) {
      list(phi = theta$phi)
    },
    if (prepared$time_varying) list(sigma_V = theta$sigma_V),
    list(re_cov_pars = par$re_cov_pars, re_b_pars = effects$re_b_pars),
    if (prepared$smooths) {
      list(b_smooth = effects$b_smooth, ln_smooth_sigma = par$ln_smooth_sigma)
    },
    theta$threshold,
    list(sigma_O = theta$sigma_O, sigma_E = exp(effects$log_sigma_E),
      sigma_Z = theta$sigma_Z),
    if (any_field) list(rho = theta$rho),
    if (spde_field) list(range = theta$range, log_range = log(theta$range)),
    if (areal_type == "sar") list(rho_sar = theta$rho_sar),
    if (areal_type == "car") list(alpha_car = theta$alpha_car),
    if (prepared$anisotropy) list(H = theta$H[[1L]]),
    if (prepared$anisotropy && n_m > 1L) list(H2 = theta$H[[2L]]),
    if (prepared$epsilon_trend) list(b_epsilon = par$b_epsilon),
    if (prepared$mixture) {
      list(p_extreme = theta$p_extreme, mix_ratio = theta$mix_ratio)
    },
    if (length(par$thetaf)) list(tweedie_p = theta$tweedie_p),
    if (length(par$ln_student_df)) list(student_df = theta$student_df),

    # Fitted rows
    list(eta_fixed_i = fitted$fixed, eta_i = fitted$eta,
      eta_smooth_i = fitted$smooth, eta_rw_i = fitted$rw,
      eta_iid_re_i = fitted$iid,
      covariate_diffusion_values = fitted$diffusion,
      omega_s_A = fitted$omega, epsilon_st_A_vec = fitted$epsilon,
      zeta_s_A = fitted$zeta, b_rw_t = effects$b_rw_t,
      jnll_obs = obs$jnll_obs, devresid = obs$devresid),
    if (report_y) list(y_i = y_i),
    # Like the C++ SIMULATE block, report latent fields in every simulation.
    if (in_simulation) {
      effects[c("omega_s", "epsilon_st", "zeta_s")]
    },
    if (prepared$rsr) {
      list(b_j_prime = rtmb_rsr_coefficients(par$b_j, fitted, prepared, 1L))
    },
    if (prepared$rsr && n_m > 1L) {
      list(b_j2_prime = rtmb_rsr_coefficients(par$b_j2, fitted, prepared, 2L))
    },

    # Projected rows and derived indices
    if (!is.null(projected)) {
      list(proj_fe = projected$fe, proj_eta = projected$eta,
        proj_rf = projected$rf, proj_omega_s_A = projected$omega,
        proj_epsilon_st_A_vec = projected$epsilon,
        proj_zeta_s_A = projected$zeta, proj_rw_i = projected$rw,
        proj_iid_re_i = projected$iid,
        proj_covariate_diffusion_values = projected$diffusion)
    },
    if (!is.null(projected) && n_m > 1L) {
      list(proj_fe_combined = projected$combined$fe,
        proj_eta_combined = projected$combined$eta,
        proj_response_combined = projected$combined$response)
    },
    if (requested[["total"]]) derived["link_total"],
    if (requested[["weighted_avg"]]) derived["weighted_avg"],
    if (requested[["eao"]]) derived[c("eao", "mean_dens")],
    diffusion_scales
  )
  # Standard errors for population-level or full projected predictions.
  se_name <- if (prepared$pop_pred) "proj_fe" else "proj_eta"
  # Reported values that also get standard errors; sdreport() is read by name.
  with_se <- c("sigma_O", "sigma_E", "sigma_Z", "sigma_V", "re_cov_pars",
    "re_b_pars", "b_j_prime", "b_j2_prime", names(theta$threshold),
    "b_epsilon", "rho_sar", "alpha_car", "range",
    "log_range", names(diffusion_scales), "phi", "tweedie_p", "student_df",
    "link_total", "weighted_avg", "eao",
    if (any(prepared$temporal & prepared$epsilon_ar1)) "rho",
    if (prepared$adreport_projection) c(se_name, paste0(se_name, "_combined")))
  adreports <- c(
    reports[intersect(names(reports), with_se)],
    list(log_sigma_O = theta$log_sigma_O, log_sigma_E = effects$log_sigma_E,
      log_sigma_Z = log(theta$sigma_Z)),
    if (prepared$dispersion_model && length(par$b_disp_k) == 1L) {
      list(`ln_phi_i(0)` = obs$ln_phi_i[[1L]], `phi_i(0)` = obs$phi_i[[1L]])
    },
    if (prepared$mixture) par[c("logit_p_extreme", "log_ratio_mix")],
    if (requested[["total"]]) derived["total"],
    if (requested[["eao"]]) derived["log_eao"]
  )
  rtmb_register_reports(reports, RTMB::REPORT)
  rtmb_register_reports(adreports, RTMB::ADREPORT)
  invisible(NULL)
}

# Call `REPORT()` or `ADREPORT()` on each element of a named list, so each
# value is recorded under its list name.
rtmb_register_reports <- function(values, fun) {
  for (name in names(values)) eval(as.call(list(fun, as.name(name))), values)
}

# Restricted spatial regression (Hanks et al. 2015): fixed effects adjusted
# for their projection onto the summed random fields at fitted rows.
rtmb_rsr_coefficients <- function(b, fitted, prepared, m) {
  fields <- fitted$omega[, m] + fitted$epsilon[, m] + fitted$svc[, m]
  b + rtmb_product(prepared$rsr_projection[[m]], fields)
}
