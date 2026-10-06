#' Length-stratified age sampling
#'
#' Describe a length-stratified age sample for `sdmTMB(length_stratified =
#' ...)`. In such a sample, ages are read for a fixed number of fish per length
#' bin, so aged fish over-represent lengths in bins where few fish were
#' measured, and length-at-age fitted to them as a random sample is biased.
#' With this, the likelihood of each aged fish's length is conditioned on its
#' age and on it having been selected for ageing (Candy et al. 2007):
#'
#' \deqn{f(l_i \mid a_i, \textrm{aged}) = \frac{\pi_{k(l_i)} f(l_i \mid a_i)}
#'   {\sum_k \pi_k \Pr(l \in \textrm{bin } k \mid a_i)},}
#'
#' where \eqn{\pi_k} is the fraction of measured fish that were aged in length
#' bin \eqn{k} of the fish's sampling unit and \eqn{f(l_i \mid a_i)} is the
#' model's length distribution for fish \eqn{i}. The numerator's
#' \eqn{\pi_{k(l_i)}} does not depend on parameters, so only the denominator,
#' the probability the fish would have been aged, enters the fit.
#'
#' Only the RTMB backend and the [gaussian()], [lognormal()], [Gamma()], and
#' [gengamma()] families are supported. Predictions are of the population, not
#' of the aged sample. Residuals, response simulations, and DHARMa diagnostics
#' are not yet supported for these fits. Gaussian and lognormal bin
#' probabilities use stable log-tail calculations. Gamma and generalized-gamma
#' support is provisional: their ordinary CDF differences can lose precision
#' in extreme tails and cause non-finite likelihoods or gradients. Generalized
#' gamma can also be unstable near its lognormal limit (`Q = 0`).
#'
#' This implements Candy's conditional-on-age method, not the joint age-length
#' empirical proportion (EP) likelihood. Sampling fractions are treated as
#' fixed, and empty length strata are excluded from the selection denominator.
#'
#' @param unit Name of the column in both `data` and `counts` identifying the
#'   sampling unit within which ages were sampled by length bin (e.g., a survey
#'   set or a survey-year).
#' @param breaks Increasing length bin edges on the scale of the response. The
#'   lowest and highest bins are open-ended. If lengths are recorded rounded
#'   (e.g., to the nearest cm), place edges between recorded values (e.g., 9.5,
#'   11.5) so that bins match true lengths.
#' @param counts A data frame with the `unit` column, `length` (any length in
#'   the bin, e.g., a recorded length or the bin's lower edge), `n_measured`,
#'   and `n_aged`. Rows with the same unit and bin are summed. Counts must be
#'   finite, nonnegative integers, with `n_aged <= n_measured`. Every sampling
#'   unit in `data` must occur in `counts`, and every bin containing fitted
#'   fish must have positive measured and aged counts. Omitted bins within a
#'   represented unit are treated as empty, with sampling fraction zero,
#'   following Candy's method as described by Perreault et al. (2020).
#'
#' @return An object to pass to `sdmTMB(length_stratified = ...)`.
#' @references
#' Candy, S.G., Constable, A.J., Lamb, T., and Williams, R. 2007. A von
#' Bertalanffy growth model for toothfish at Heard Island fitted to
#' length-at-age data and compared to observed growth from mark-recapture
#' studies. CCAMLR Science, 14: 43–66.
#'
#' Perreault, A.M.J., Zheng, N., and Cadigan, N.G. 2020. Estimation of growth
#' parameters based on length-stratified age samples. Canadian Journal of
#' Fisheries and Aquatic Sciences, 77(3): 439–450.
#' \doi{10.1139/cjfas-2019-0129}
#' @export
lsas <- function(unit, breaks, counts) {
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    cli_abort("`unit` must be the name of a column.")
  }
  if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) ||
      is.unsorted(breaks, strictly = TRUE)) {
    cli_abort("`breaks` must be at least two increasing numbers.")
  }
  needed <- c(unit, "length", "n_measured", "n_aged")
  if (!is.data.frame(counts) || !all(needed %in% names(counts))) {
    cli_abort("`counts` must be a data frame with columns {.val {needed}}.")
  }
  if (anyNA(counts[[unit]])) cli_abort("`counts` sampling units can't have missing values.")
  if (!is.numeric(counts$length) || any(!is.finite(counts$length))) {
    cli_abort("`counts$length` must contain finite numbers.")
  }
  for (column in c("n_measured", "n_aged")) {
    z <- counts[[column]]
    if (!is.numeric(z) || any(!is.finite(z) | z < 0 | z != floor(z))) {
      cli_abort("`counts${column}` must contain finite, nonnegative integers.")
    }
  }
  if (any(counts$n_aged > counts$n_measured)) {
    cli_abort("`counts$n_aged` can't exceed `counts$n_measured`.")
  }
  structure(list(unit = unit, breaks = breaks, counts = counts),
    class = "sdmTMB_lsas")
}

# Length bin of each length; lengths beyond the edges fall in the open-ended
# end bins.
.lsas_bin <- function(length, breaks) {
  n_bins <- length(breaks) - 1L
  pmin(pmax(findInterval(length, breaks), 1L), n_bins)
}

# Data for the RTMB objective: the unit-by-bin matrix of sampling fractions,
# the interior bin edges, and each fish's unit.
.lsas_tmb_data <- function(x, data, y, family_spec, backend) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, "sdmTMB_lsas")) {
    cli_abort("`length_stratified` must be created with `lsas()`.")
  }
  if (backend != "rtmb") {
    cli_abort("`length_stratified` needs `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  family <- family_spec$family$family
  if (length(family) != 1L || !family %in% c("gaussian", "lognormal", "Gamma", "gengamma")) {
    cli_abort("`length_stratified` needs the `gaussian()`, `lognormal()`, `Gamma()`, or `gengamma()` family.")
  }
  if (family %in% c("Gamma", "gengamma")) {
    cli_inform(c("i" = paste0(
      "Length-stratified sampling with `", family, "()` is provisional: ",
      "bin probabilities use ordinary CDF differences, which can lose precision ",
      "in extreme tails and cause non-finite likelihoods or gradients.",
      if (family == "gengamma") " Generalized gamma can also be unstable near Q = 0."
    )))
  }
  if (!x$unit %in% names(data)) {
    cli_abort("Column {.val {x$unit}} is missing from `data`.")
  }
  if (anyNA(data[[x$unit]])) cli_abort("Sampling units in `data` can't have missing values.")
  units <- unique(data[[x$unit]])
  if (any(!units %in% x$counts[[x$unit]])) {
    cli_abort("Every sampling unit in `data` must occur in `counts`.")
  }
  unit_i <- match(data[[x$unit]], units)

  counts <- x$counts[x$counts[[x$unit]] %in% units, , drop = FALSE]
  n_bins <- length(x$breaks) - 1L
  unit_c <- factor(match(counts[[x$unit]], units), levels = seq_along(units))
  bin_c <- factor(.lsas_bin(counts$length, x$breaks), levels = seq_len(n_bins))
  measured <- unclass(stats::xtabs(counts$n_measured ~ unit_c + bin_c))
  aged <- unclass(stats::xtabs(counts$n_aged ~ unit_c + bin_c))
  pi <- ifelse(measured > 0, aged / measured, 0)
  dimnames(pi) <- NULL

  if (any(pi[cbind(unit_i, .lsas_bin(y, x$breaks))] == 0)) {
    cli_abort("Some fish in `data` fall in a unit and length bin with no aged fish in `counts`.")
  }
  list(pi = pi, cuts = x$breaks[-c(1L, n_bins + 1L)], unit_i = unit_i)
}

# Log probability that each fitted fish `i` was selected for ageing given its
# age: the sum over length bins of the bin's sampling fraction times the
# probability of the bin under the fish's length distribution `s`.
rtmb_lsas_log_selection <- function(s, family, lsas, i) {
  "[<-" <- RTMB::ADoverload("[<-")
  pi <- lsas$pi[lsas$unit_i[i], , drop = FALSE]
  # WIP: RTMB::pgamma has no AD log-tail interface. Retain the original
  # calculation for Gamma/gengamma, with an explicit notice at fit setup.
  if (family %in% c("Gamma", "gengamma")) {
    selected <- 0
    below <- 0
    for (k in seq_along(lsas$cuts)) {
      cdf <- rtmb_lsas_cdf(lsas$cuts[[k]], s, family)
      selected <- selected + pi[, k] * (cdf - below)
      below <- cdf
    }
    return(log(selected + pi[, ncol(pi)] * (1 - below)))
  }
  edges <- c(-Inf, lsas$cuts, Inf)
  selected <- rep(-Inf, length(i))
  initialized <- rep(FALSE, length(i))
  for (k in seq_len(ncol(pi))) {
    rows <- which(pi[, k] > 0)
    if (!length(rows)) next
    # Positive families have no mass at or below zero.
    if (family != "gaussian" && edges[k + 1L] <= 0) next
    state <- lapply(s, rtmb_pick, i = rows)
    term <- log(pi[rows, k]) +
      rtmb_lsas_log_interval(edges[k], edges[k + 1L], state, family)
    old <- initialized[rows]
    selected[rows[!old]] <- term[!old]
    if (any(old)) {
      selected[rows[old]] <- RTMB::logspace_add(selected[rows[old]], term[old])
    }
    initialized[rows] <- TRUE
  }
  selected
}

# Choose the tail before subtracting, so a rounded-to-one CDF does not erase
# a small bin probability. AD branching avoids evaluating an unused
# logspace_sub(0, 0) branch.
rtmb_lsas_log_interval <- function(lower, upper, s, family) {
  if (is.infinite(lower) || (family != "gaussian" && lower <= 0)) {
    if (is.infinite(upper)) return(0 * s$mu)
    return(rtmb_lsas_cdf(upper, s, family, log.p = TRUE))
  }
  if (is.infinite(upper)) {
    return(rtmb_lsas_cdf(lower, s, family, lower.tail = FALSE, log.p = TRUE))
  }
  lo <- rtmb_lsas_cdf(lower, s, family, log.p = TRUE)
  hi <- rtmb_lsas_cdf(upper, s, family, log.p = TRUE)
  slo <- rtmb_lsas_cdf(lower, s, family, lower.tail = FALSE, log.p = TRUE)
  shi <- rtmb_lsas_cdf(upper, s, family, lower.tail = FALSE, log.p = TRUE)
  RTMB::Vectorize(function(lo, hi, slo, shi) {
    "if" <- RTMB::ADoverload("if")
    if (hi < log(0.5)) RTMB::logspace_sub(hi, lo)
    else RTMB::logspace_sub(slo, shi)
  })(lo, hi, slo, shi)
}

# CDF at length `q`, parameterized as in `rtmb_obs_families`. Gaussian and
# lognormal tails are evaluated directly; Gamma/gengamma remain provisional.
rtmb_lsas_cdf <- function(q, s, family, lower.tail = TRUE, log.p = FALSE) {
  if (family != "gaussian" && q <= 0) {
    p <- if (lower.tail) 0 else 1
    return(rep(if (log.p) log(p) else p, length(s$mu)))
  }
  p <- switch(family,
    gaussian = return(RTMB::pnorm((q - s$mu) / s$phi,
      lower.tail = lower.tail, log.p = log.p)),
    lognormal = return(RTMB::pnorm((log(q) - log(s$mu) + s$phi^2 / 2) / s$phi,
      lower.tail = lower.tail, log.p = log.p)),
    # WIP: these ordinary CDFs are used only by the provisional path above.
    Gamma = RTMB::pgamma(q, shape = s$phi, scale = s$mu / s$phi),
    gengamma = {
      k <- s$Q^-2
      p <- RTMB::pgamma(k * exp(rtmb_gengamma_qw(q, s$mu, s$phi, s$Q)), shape = k)
      0.5 + sign(s$Q) * (p - 0.5)
    })
  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

.check_lsas_diagnostics <- function(object, caller) {
  if (!is.null(object$tmb_data$lsas)) {
    cli_abort(c(
      "{caller} is not yet supported for `length_stratified` fits.",
      "i" = "Residuals and response simulations must account for selection for ageing. Population predictions are available with `predict()`."
    ))
  }
}
