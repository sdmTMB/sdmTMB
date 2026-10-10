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
#' Spatially varying coefficient SDs are returned as `sigma_Z`. If any
#' coefficient has its own Matérn range (see `range_groups` in [sdmTMB()]),
#' their ranges are returned as `range_Z`, and `range` is omitted if there are
#' no spatial or spatiotemporal fields. With several coefficients, these names
#' are suffixed with the coefficient name, e.g., `sigma_Z_depth_scaled`.
#' Ranges that are fixed (mapped off) are not returned.
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
  if (!"jointPrecision" %in% names(object$sd_report)) {
    cli_abort("TMB::sdreport() must be run with the joint precision returned.")
  }

  if (.object_has_two_components(object, caller = "`spread_sims()`")) {
    cli_abort("This function isn't yet set up for delta models.")
  }
  n_sims <- nsim
  tmb_sd <- object$sd_report
  samps <- rmvnorm_prec(object$tmb_obj$env$last.par.best, tmb_sd, n_sims)
  pars <- c(tmb_sd$par.fixed, tmb_sd$par.random)
  pn <- names(pars)

  pn <- c(pn[pn == "b_j"], pn[pn == "b_j2"], pn[!pn %in% c("b_j", "b_j2")]) # if REML, must move b_j to beginning:
  par_fixed <- names(tmb_sd$par.fixed)
  if (object$reml) par_fixed <- c("b_j", par_fixed)
  pars <- pars[pn %in% par_fixed]
  samps <- samps[pn %in% par_fixed, , drop = FALSE]
  pn <- pn[pn %in% par_fixed]
  .formula <- object$split_formula[[1]]$form_no_bars # TODO DELTA HARDCODED TO 1!
  if (isFALSE(object$mgcv)) {
    fe_names <- colnames(model.matrix(.formula, object$data))
  } else {
    fe_names <- colnames(model.matrix(mgcv::gam(.formula, data = object$data)))
  }
  fe_names <- tidy(object, silent = TRUE)$term
  row.names(samps) <- pn
  row.names(samps)[row.names(samps) == "b_j"] <- fe_names
  out <- as.data.frame(t(samps))
  is_areal <- is_areal_fit(object)
  is_car <- is_car_fit(object)
  has_ln_kappa <- "ln_kappa" %in% names(out)
  par_list <- object$tmb_obj$env$parList()
  # Draws of entry `r` of parameter `par`: entries tied by the map share one
  # column, and entries mapped off take their fixed value
  draws <- function(par, r) {
    cols <- which(names(out) == par)
    k <- if (is.null(object$tmb_map[[par]])) r else
      as.integer(object$tmb_map[[par]])[r]
    if (is.na(k)) rep(as.vector(par_list[[par]])[r], n_sims) else out[[cols[k]]]
  }
  # Areal fields have no range: their SD is exp(-ln_tau)
  use_kappa <- !is_areal && !is.null(par_list$ln_kappa)
  ln_kappa <- function(r) if (use_kappa) draws("ln_kappa", r)
  kappa_estimated <- function(r) {
    has_ln_kappa && (is.null(object$tmb_map$ln_kappa) ||
      !is.na(object$tmb_map$ln_kappa[r]))
  }
  matern_sd <- function(ln_tau, ln_kappa) {
    if (!use_kappa) return(exp(-ln_tau))
    1 / sqrt(4 * pi * exp(2 * ln_tau + 2 * ln_kappa))
  }

  # The spatial range row is mapped off if no field uses it (e.g., only SVCs
  # with their own ranges)
  if (has_ln_kappa && !is_areal && kappa_estimated(1L)) {
    out$range <- sqrt(8) / exp(ln_kappa(1L))
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
    out$sigma_O <- matern_sd(out$ln_tau_O, ln_kappa(1L))
  }
  if ("ln_tau_E" %in% names(out)) {
    out$sigma_E <- matern_sd(out$ln_tau_E, ln_kappa(2L))
  }
  sims_z <- list()
  if ("ln_tau_Z" %in% names(out)) {
    # One SD per SVC and, if any SVC has its own range, one estimated range
    # per SVC (as in `tidy()`): named `sigma_Z` and `range_Z` for a single SVC
    # and suffixed with the coefficient name (e.g., `sigma_Z_depth_scaled`)
    # for several
    fields <- range_fields(object)
    svc <- fields[fields$type == "svc" & fields$component == 1L, ]
    n_z <- nrow(svc)
    suffix <- if (n_z > 1L) paste0("_", svc$field) else ""
    svc_ranges <- use_kappa && any(svc$kappa_row != 1L)
    for (z in seq_len(n_z)) {
      sims_z[[paste0("sigma_Z", suffix[z])]] <-
        matern_sd(draws("ln_tau_Z", z), ln_kappa(svc$kappa_row[z]))
    }
    for (z in seq_len(n_z)) {
      r <- svc$kappa_row[z]
      if (svc_ranges && kappa_estimated(r)) {
        sims_z[[paste0("range_Z", suffix[z])]] <- sqrt(8) / exp(ln_kappa(r))
      }
    }
  }
  if ("ln_tau_O_trend" %in% names(out)) {
    out$sigma_O_trend <- matern_sd(out$ln_tau_O_trend, ln_kappa(1L))
  }
  # Remove internal columns (duplicated names included) only after all the
  # transformations
  out <- out[!names(out) %in% c("ln_kappa", "ln_tau_Z")]
  if (length(sims_z)) out <- cbind(out, as.data.frame(sims_z))
  out$ln_tau_O <- out$ln_tau_E <- out$ln_tau_O_trend <-
    out$ar1_phi <- out$thetaf <- out$ln_phi <- out$logit_rho_sar <- NULL
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
