# Translate the TMB data list into named feature flags, one-based indices, and
# row inputs for the RTMB objective. Only this function and the builders it
# calls (`rtmb_*_inputs()` and the helpers below) read the C++ data list; the
# objective reads `prepared` alone, so the data list built by `fit.R` can
# change here without touching the other modules. Nothing here depends on
# parameters.
rtmb_prepare <- function(data) {
  n_m <- ncol(data$y_i)
  components <- seq_len(n_m)
  on <- function(flag) isTRUE(flag == 1L)
  families <- rtmb_family_inputs(data)
  any_field <- data$no_spatial == 0L
  include_spatial <- data$include_spatial[components] == 1L
  svc <- on(data$spatial_covariate)
  smooths <- on(data$has_smooths)
  time_varying_type <- if (on(data$ar1_time)) "ar1" else
    c("none", "rw", "rw0")[data$random_walk + 1L]
  prepared <- list(
    n_m = n_m,
    n_t = data$n_t,
    families = families,
    # Mixture families are single-family models; C++ checks family 1.
    mixture = families[[1L]]$family[[n_m]] %in% rtmb_mixture_families,
    dispersion_model = on(data$has_dispersion_model),

    # Random fields. Logical vectors have one entry per component.
    any_field = any_field,
    precision = rtmb_precision_inputs(data),
    anisotropy = any_field && on(data$anisotropy),
    include_spatial = include_spatial,
    spatial = any_field & include_spatial & data$omit_spatial_intercept == 0L,
    temporal = any_field & data$spatial_only[components] == 0L,
    epsilon_ar1 = data$ar1_fields[components] == 1L,
    epsilon_rw = data$rw_fields[components] == 1L,
    share_range = data$share_range[components] == 1L,
    epsilon_trend = on(data$est_epsilon_slope),
    epsilon_predictor = data$epsilon_predictor,
    # SVC fields enter the predictor whenever present; like the C++ template,
    # their density is only evaluated for components with a spatial field.
    svc = svc,
    svc_density = svc & include_spatial,

    # Other effects and predictor terms
    time_varying = time_varying_type != "none",
    time_varying_type = time_varying_type,
    iid_re = data$n_re_groups[components] > 0L,
    re_groups = lapply(components, function(m) rtmb_re_groups(data, m)),
    smooths = smooths,
    smooth_index = if (smooths) rtmb_smooth_index(data) else list(),
    threshold = c("none", "breakpt", "logistic")[data$threshold_func + 1L],
    diffusion = rtmb_diffusion_inputs(data),
    priors = rtmb_prior_inputs(data),

    # Simulation: which time steps, latent effects, and whether the response
    # is a simulation target (otherwise `y_i` reports response means).
    simulate_t = which(data$simulate_t == 1L),
    simulate_re = stats::setNames(data$sim_re == 1L, c("omega_s",
      "epsilon_st", "zeta_s", "re_b_pars", "b_rw_t", "b_smooth")),
    simulate_obs = on(data$sim_obs),

    # Reports
    rsr = on(data$do_rsr),
    pop_pred = on(data$pop_pred),
    adreport_projection = on(data$calc_se),
    derived = on(data$do_predict) & c(total = on(data$calc_index_totals),
      weighted_avg = on(data$calc_weighted_avg), eao = on(data$calc_eao)),
    index = list(link = rtmb_code_name(data$link_pred, .valid_link),
      time_include = data$proj_time_include == 1L)
  )
  prepared$fit <- rtmb_row_inputs(data, prepared$families, projection = FALSE)
  prepared$obs_groups <- rtmb_obs_groups(data, prepared$families,
    prepared$fit$family_id)
  if (prepared$rsr) {
    # (X'X)^-1 X' for restricted spatial regression; depends only on data.
    prepared$rsr_projection <- lapply(prepared$fit$X, function(X) {
      X <- as.matrix(X)
      solve(crossprod(X), t(X))
    })
  }
  if (on(data$do_predict)) {
    prepared$proj <- rtmb_row_inputs(data, prepared$families,
      projection = TRUE)
  }
  if (has_preferential(data)) {
    prepared$preferential <- rtmb_preferential_inputs(data, prepared$families)
  }
  prepared
}

has_preferential <- function(data) isTRUE(data$preferential$n_pref > 0L)

# Sampling-frame inputs: the observed indicators and sampling design, and
# `rows`, the frame's shared catch-predictor rows for
# rtmb_linear_predictors(). Like projection rows, these select from the
# unique frame locations. Terms that preferential sampling doesn't support
# yet (smoothers, SVCs, thresholds, time-varying, and diffusion) are rejected
# before this point and have no inputs here. The optional sampling field
# `xi` always uses an isotropic SPDE precision on the model mesh, even in the
# first multiphase fit, when the catch fields are off.
rtmb_preferential_inputs <- function(data, families) {
  pref <- data$preferential
  station_index <- pref$station_i + 1L
  rows <- list(
    X = pref$X_ij, offset = pref$offset_i,
    Zt = list(), include_iid = pref$include_iid == 1L,
    A_rows = pref$A_station[station_index, , drop = FALSE],
    A_station = pref$A_station, station_index = station_index,
    time = pref$year_i + 1L,
    family_id = rep(1L, pref$n_pref)
  )
  list(
    rows = rtmb_row_family_flags(rows, families),
    R = pref$R_i,
    observed = which(!is.na(pref$R_i)),
    Z = pref$Z_ij,
    xi = pref$spatial_xi == 1L,
    precision = if (pref$spatial_xi == 1L) {
      rtmb_precision_inputs(list(spatial_model = 0L, no_spatial = 0L,
        barrier = 0L, anisotropy = 0L, spde = data$spde))
    }
  )
}

# Names for integer codes from `R/enum.R`; unmatched codes give NA.
rtmb_code_name <- function(code, codes) names(codes)[match(code, codes)]

# One entry per observation family: how its components combine ("single",
# "delta", or "poisson_link"), the family and link name of each component, and
# one-based auxiliary parameter slots.
rtmb_family_inputs <- function(data) {
  slot <- function(x, f) if (x[[f]] >= 0L) x[[f]] + 1L else NA_integer_
  lapply(seq_along(data$combine_kind), function(f) {
    combine <- c("single", "delta", "poisson_link")[data$combine_kind[[f]] + 1L]
    list(
      combine = combine,
      active = data$component_active[f, ] == 1L,
      family = rtmb_code_name(data$family_code[f, ], .valid_family),
      link = rtmb_code_name(data$link_code[f, ], .valid_link),
      # Offsets enter a single family's only component and a standard delta
      # family's positive component. Poisson-link deltas apply them in the
      # response mean instead.
      offset_applies = c(combine == "single",
        combine == "delta")[seq_len(ncol(data$y_i))],
      phi = slot(data$ln_phi_slot, f),
      thetaf = slot(data$thetaf_slot, f),
      student_df = slot(data$ln_student_df_slot, f),
      gengamma_Q = slot(data$gengamma_Q_slot, f)
    )
  })
}

# Covariates and projection matrices for fitted or projected rows. The fitted
# spatial fields use `A_st` directly, as in the C++ template; projection rows
# select from the unique projection locations. Fitted rows also carry the
# response and observation inputs; projected rows the index area weights.
rtmb_row_inputs <- function(data, families, projection) {
  if (projection) {
    station_index <- data$proj_spatial_index + 1L
    out <- list(
      X = data$proj_X_ij, offset = data$proj_offset_i,
      diffusion_x = data$covariate_diffusion$proj_covariate_vertex_time,
      X_threshold = data$proj_X_threshold, Zs = data$proj_Zs,
      Xs = data$proj_Xs, z = data$proj_z_i, X_rw = data$proj_X_rw_ik,
      Zt = data$Zt_list_proj, include_iid = data$exclude_RE == 0L,
      A_rows = data$proj_mesh[station_index, , drop = FALSE],
      A_station = data$proj_mesh, station_index = station_index,
      time = data$proj_year + 1L,
      family_id = rep_len(data$proj_family_id + 1L, nrow(data$proj_X_ij[[1L]])),
      area = data$area_i, weight = data$proj_vector
    )
  } else {
    out <- list(
      X = data$X_ij, offset = data$offset_i,
      diffusion_x = data$covariate_diffusion$covariate_vertex_time,
      X_threshold = data$X_threshold, Zs = data$Zs, Xs = data$Xs,
      z = data$z_i, X_rw = data$X_rw_ik,
      Zt = data$Zt_list, include_iid = TRUE,
      A_rows = data$A_st, A_station = data$A_st,
      # `A_spatial_index` can be longer than the data when the mesh is unused
      # (e.g., the placeholder mesh with no spatial fields); as in C++, read
      # only one entry per row
      station_index = data$A_spatial_index[seq_len(nrow(data$y_i))] + 1L,
      time = data$year_i + 1L,
      family_id = data$obs_family_id + 1L,
      y = data$y_i, size = data$size, weights = data$weights_i,
      upr = rep_len(data$upr, nrow(data$y_i)), Xdisp = data$Xdisp_ij
    )
  }
  rtmb_row_family_flags(out, families)
}

# Per-row component activity and offset placement from each row's family.
rtmb_row_family_flags <- function(rows, families) {
  rows$active <- do.call(rbind, lapply(families, `[[`, "active"))[
    rows$family_id, , drop = FALSE]
  rows$offset_applies <- do.call(rbind,
    lapply(families, `[[`, "offset_applies"))[rows$family_id, , drop = FALSE]
  rows
}

# Fitted rows grouped by component and family. `observed` excludes missing
# responses, which contribute no likelihood but are still simulated.
rtmb_obs_groups <- function(data, families, family_id) {
  lapply(seq_len(ncol(data$y_i)), function(m) {
    groups <- lapply(seq_along(families), function(f) {
      if (!families[[f]]$active[[m]]) return(NULL)
      rows <- which(family_id == f)
      list(f = f, rows = rows, observed = rows[!is.na(data$y_i[rows, m])])
    })
    Filter(Negate(is.null), groups)
  })
}

# Parameter rows and coefficient positions for each IID random-effect group
# of component `m`. Both components of a delta model must have identical
# random-effect terms, so, like the C++ template, component-local start and
# end positions index the combined data frames.
rtmb_re_groups <- function(data, m) {
  first <- sum(data$n_re_groups[seq_len(m - 1L)])
  lapply(first + seq_len(data$n_re_groups[[m]]), function(g) {
    dimension <- data$re_cov_df_map[g, 2L]
    cov_rows <- seq.int(data$re_cov_df_map[g, 3L] + 1L,
      data$re_cov_df_map[g, 4L] + 1L)
    is_sd <- data$re_cov_df[cov_rows, 4L] == 1L
    sd_rows <- cov_rows[is_sd]
    sd_rows <- sd_rows[order(data$re_cov_df[sd_rows, 3L])]
    corr_rows <- cov_rows[!is_sd]
    if (dimension > 1L) {
      # TMB consumes correlation parameters by row; RTMB's lower-triangle
      # constructor fills by column. A four-dimensional group exposes it.
      lower <- which(lower.tri(diag(dimension)), arr.ind = TRUE)
      rowwise <- do.call(rbind, lapply(seq.int(2L, dimension),
        function(row) cbind(row, seq_len(row - 1L))))
      corr_rows <- corr_rows[match(paste(lower[, 1L], lower[, 2L]),
        paste(rowwise[, 1L], rowwise[, 2L]))]
    }
    levels <- seq.int(data$re_b_map[g, 2L] + 1L, data$re_b_map[g, 3L] + 1L)
    list(
      dimension = dimension, sd_rows = sd_rows, corr_rows = corr_rows,
      coefficients = lapply(levels, function(level) {
        seq.int(data$re_b_df[level, 1L] + 1L, data$re_b_df[level, 2L] + 1L)
      })
    )
  })
}

rtmb_smooth_index <- function(data) {
  lapply(seq_along(data$Zs), function(s) {
    data$b_smooth_start[[s]] + seq_len(ncol(data$Zs[[s]]))
  })
}

# Reject untranslated features before RTMB tapes an incomplete model.
rtmb_validate <- function(data, prepared, parameters, random, ...) {
  if (!all(names(list(...)) %in% c("intern", "inner.control"))) {
    cli::cli_abort("Additional MakeADFun options are not supported by the RTMB backend yet.")
  }
  # Random parameters must be translated effects. The first multiphase fit
  # integrates none of them.
  expected_random <- c(
    if (any(prepared$spatial)) "omega_s",
    if (any(prepared$temporal)) "epsilon_st",
    if (prepared$svc) "zeta_s",
    if (prepared$time_varying) "b_rw_t",
    if (any(prepared$iid_re)) "re_b_pars",
    if ("b_j" %in% random) "b_j",
    if ("b_j2" %in% random) "b_j2",
    if (prepared$smooths && "bs" %in% random) "bs",
    if (prepared$smooths) "b_smooth",
    if (isTRUE(prepared$preferential$xi)) "xi_s"
  )
  unexpected <- setdiff(random, expected_random)
  if (length(unexpected)) {
    cli::cli_abort("The RTMB backend can't integrate over {.val {unexpected}}.")
  }
  if (prepared$dispersion_model && length(parameters$ln_phi) == 0L) {
    cli::cli_abort("The RTMB backend needs a family with a dispersion parameter for a dispersion model.")
  }
  if (prepared$dispersion_model && !is.null(prepared$priors$phi)) {
    cli::cli_abort("The RTMB backend doesn't support a prior on phi with a dispersion model.")
  }
  active <- data$component_active == 1L
  tweedie <- data$family_code[active] == .valid_family[["tweedie"]]
  if (!is.null(prepared$priors$tweedie_p) && any(tweedie)) {
    cli::cli_abort("Priors not enabled for Tweedie p currently.")
  }
}
