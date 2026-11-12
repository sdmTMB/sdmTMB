#' Extract parameter simulations from the joint precision matrix
#'
#' `spread_sims()` returns a wide-format data frame. `gather_sims()` returns a
#' long-format data frame. The format matches the format in the \pkg{tidybayes}
#' `spread_draws()` and `gather_draws()` functions.
#'
#' @param object Output from [sdmTMB()].
#' @param nsim The number of simulation draws.
#'
#' @export
#' @rdname gather_sims
#'
#' @return
#' A data frame. `gather_sims()` returns a long-format data frame:
#'
#' * `.iteration`: the sample ID
#' * `.variable`: the parameter name
#' * `.value`: the parameter sample value
#'
#' `spread_sims()` returns a wide-format data frame:
#'
#' * `.iteration`: the sample ID
#' * columns for each parameter with a sample per row
#'
#' @examples
#' m <- sdmTMB(density ~ depth_scaled,
#'   data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie())
#' head(spread_sims(m, nsim = 10))
#' head(gather_sims(m, nsim = 10))
#' samps <- gather_sims(m, nsim = 1000)
#'
#' if (require("ggplot2", quietly = TRUE)) {
#'   ggplot(samps, aes(.value)) + geom_histogram() +
#'     facet_wrap(~.variable, scales = "free_x")
#' }

spread_sims <- function(object, nsim = 200) {
  .check_family_capability(object, "spread_sims")
  if (.object_has_two_components(object, caller = "`spread_sims()`")) {
    cli_abort("This function isn't yet set up for delta models.")
  }
  n_sims <- nsim
  tmb_sd <- sd_report_joint(object)
  samps <- draw_joint_pars(object, tmb_sd, n_sims)
  # rows of `samps` follow `last.par.best`; keep the fixed effects, including
  # those REML integrates out
  pn <- names(object$tmb_obj$env$last.par.best)
  keep <- !seq_along(pn) %in% object$tmb_obj$env$random |
    (object$reml & pn %in% c("b_j", "bs"))
  samps <- samps[keep, , drop = FALSE]
  pn <- pn[keep]
  # tidy() lists the b_j terms first, then any threshold and smoother terms
  terms <- tidy(object, silent = TRUE)$term
  b_names <- estimated_names(terms[seq_along(object$parlist$b_j)], object$tmb_map$b_j)
  bs_names <- if (object$tmb_data$has_smooths) {
    rownames(print_smooth_effects(object, silent = TRUE)$smooth_effects)
  }
  pn[pn == "b_j"] <- b_names
  if (sum(pn == "bs") == length(bs_names)) pn[pn == "bs"] <- bs_names
  row.names(samps) <- pn
  out <- as.data.frame(t(samps))
  is_areal <- is_areal_fit(object)
  is_car <- is_car_fit(object)
  has_ln_kappa <- "ln_kappa" %in% names(out)

  if (has_ln_kappa && !is_areal) {
    out$range <- sqrt(8) / exp(out$ln_kappa)
  }
  if ("ln_phi" %in% names(out)) {
    out$phi <- exp(out$ln_phi)
  }
  if ("thetaf" %in% names(out)) {
    out$tweedie_p <- stats::plogis(out$thetaf) + 1
  }
  if ("ar1_phi" %in% names(out)) {
    out$ar1_rho <- 2 * stats::plogis(out$ar1_phi) - 1
  }
  if ("logit_rho_sar" %in% names(out)) {
    if (is_car) {
      out$alpha_car <- stats::plogis(out$logit_rho_sar)
    } else {
      out$rho_sar <- 2 * stats::plogis(out$logit_rho_sar) - 1
    }
  }
  if ("ln_tau_O" %in% names(out)) {
    if (is_areal || !has_ln_kappa) {
      out$sigma_O <- exp(-out$ln_tau_O)
    } else {
      out$sigma_O <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_O + 2 * out$ln_kappa))
    }
  }
  if ("ln_tau_E" %in% names(out)) {
    if (is_areal || !has_ln_kappa) {
      out$sigma_E <- exp(-out$ln_tau_E)
    } else {
      out$sigma_E <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_E + 2 * out$ln_kappa))
    }
  }
  if ("ln_tau_Z" %in% names(out)) {
    if (is_areal || !has_ln_kappa) {
      out$sigma_Z <- exp(-out$ln_tau_Z)
    } else {
      out$sigma_Z <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_Z + 2 * out$ln_kappa))
    }
  }
  if ("ln_tau_O_trend" %in% names(out)) {
    if (is_areal || !has_ln_kappa) {
      out$sigma_O_trend <- exp(-out$ln_tau_O_trend)
    } else {
      out$sigma_O_trend <- 1 / sqrt(4 * pi * exp(2 * out$ln_tau_O_trend + 2 * out$ln_kappa))
    }
  }
  out$ln_kappa <- out$ln_tau_O <- out$ln_tau_E <- out$ln_tau_O_trend <-
    out$ln_tau_Z <- out$ar1_phi <- out$thetaf <- out$ln_phi <- out$logit_rho_sar <- NULL
  data.frame(.iteration = seq_len(n_sims), out)
}

#' @export
#' @rdname gather_sims
gather_sims <- function(object, nsim = 200) {

  n_sims <- nsim

  out_wide <- spread_sims(object, n_sims)
  out_wide$.iteration <- NULL
  out <- stats::reshape(out_wide, direction = "long", varying = list(names(out_wide)),
    idvar = ".iteration", timevar = "variable_num")
  names(out)[2] <- ".value"
  row.names(out) <- NULL
  par_names <- data.frame(variable_num = unique(out$variable), .variable = names(out_wide))
  out <- base::merge(out, par_names)
  out$variable_num <- NULL
  out[ , c(".iteration", ".variable", ".value"), drop = FALSE]
}

# Names of the estimated elements of a parameter given the names of all its
# elements: a mapped parameter has one estimated element per `map` level
estimated_names <- function(names, map) {
  if (!length(map)) return(names)
  names[match(seq_len(nlevels(map)), as.integer(map))]
}

# `object$sd_report`, recomputing the joint precision if it wasn't stored
# (e.g., removed by strip_sdmTMB()). Fits without random effects have none.
sd_report_joint <- function(object) {
  sd_report <- object$sd_report
  if (is.null(sd_report[["jointPrecision"]]) && !has_no_random_effects(object)) {
    message("Rerunning RTMB::sdreport() with `getJointPrecision = TRUE`.")
    reinitialize(object)
    sd_report <- RTMB::sdreport(object$tmb_obj, getJointPrecision = TRUE)
  }
  sd_report
}

# Draws (parameters x draws) of all parameters around `last.par.best`: from
# the joint precision, or from `cov.fixed` without random effects
draw_joint_pars <- function(object, sd_report, n_sims) {
  mu <- object$tmb_obj$env$last.par.best
  if (!has_no_random_effects(object)) {
    return(rmvnorm_prec(mu, sd_report, n_sims))
  }
  L <- t(chol(sd_report$cov.fixed))
  z <- matrix(stats::rnorm(length(mu) * n_sims), ncol = n_sims)
  mu + as.matrix(L %*% z)
}

rmvnorm_prec <- function(mu, tmb_sd, n_sims) {
  L <- Matrix::Cholesky(tmb_sd[["jointPrecision"]], super = TRUE)
  rmvnorm_chol(mu, L, n_sims)
}

# Draws from N(mu, Q^-1) given a sparse Cholesky factor `L` of Q
rmvnorm_chol <- function(mu, L, n_sims) {
  z <- matrix(stats::rnorm(length(mu) * n_sims), ncol = n_sims)
  z <- Matrix::solve(L, z, system = "Lt")
  z <- Matrix::solve(L, z, system = "Pt")
  z <- as.matrix(z)
  mu + z
}
