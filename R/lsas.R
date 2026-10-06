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
#' the probability the fish would have been aged, enters the fit. The
#' resulting log likelihood is therefore not comparable (e.g., by AIC) with
#' that of a model fitted without `length_stratified`.
#'
#' This implements Candy's conditional-on-age method, not the joint age-length
#' empirical proportion (EP) likelihood. Sampling fractions are treated as
#' fixed. Bins with no measured fish in a unit have sampling fraction zero
#' and are excluded from that unit's selection probability, as described by
#' Perreault et al. (2020), so very fine bins with sparse tails truncate the
#' modelled length distribution.
#'
#' Only the RTMB backend and the [gaussian()] and [lognormal()] families are
#' supported. Predictions are of the population, not of the aged sample.
#' Residuals, response simulations, and DHARMa diagnostics are not yet
#' supported for these fits.
#'
#' @param measured A data frame with one row per measured fish, including
#'   those not aged.
#' @param unit Name of the column in both `measured` and the fitted `data`
#'   identifying the sampling unit within which ages were sampled by length
#'   bin (e.g., a survey set or a survey-year).
#' @param breaks Increasing length bin edges on the scale of the response. The
#'   lowest and highest bins are open-ended. If lengths are recorded rounded
#'   (e.g., to the nearest cm), place edges between recorded values (e.g., 9.5,
#'   11.5) so that bins match true lengths.
#' @param length Name of the length column in `measured`.
#' @param aged Name of the logical column in `measured` that is `TRUE` for
#'   fish with an age (typically `!is.na(age)`) and `FALSE` for fish that were
#'   measured but not aged, including any whose age could not be read.
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
#' @examples
#' # 1000 measured fish in two years; up to 5 aged per 5 cm bin and year
#' set.seed(1)
#' measured <- data.frame(year = rep(1:2, each = 500), age = rpois(1000, 3) + 1)
#' measured$length <- rlnorm(1000, log(40 * (1 - exp(-0.4 * measured$age))), 0.1)
#' bin <- findInterval(measured$length, seq(0, 60, by = 5))
#' measured$aged <- ave(seq_len(1000), measured$year, bin, FUN = seq_along) <= 5
#' measured$age[!measured$aged] <- NA # ages are only known for aged fish
#' design <- length_strata(measured, unit = "year", breaks = seq(0, 60, by = 5))
#'
#' fit <- sdmTMB(length ~ 0 + factor(age), data = measured[measured$aged, ],
#'   family = lognormal(), spatial = "off", length_stratified = design)
length_strata <- function(measured, unit, breaks, length = "length",
                          aged = "aged") {
  if (!is.data.frame(measured)) cli_abort("`measured` must be a data frame.")
  for (column in c(unit, length, aged)) {
    if (!is.character(column) || length(column) != 1L ||
        !column %in% names(measured)) {
      cli_abort("{.val {column}} must name a column of `measured`.")
    }
  }
  if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) ||
      is.unsorted(breaks, strictly = TRUE)) {
    cli_abort("`breaks` must be at least two increasing numbers.")
  }
  if (anyNA(measured[[unit]])) {
    cli_abort("Sampling units in `measured` can't have missing values.")
  }
  if (!is.numeric(measured[[length]]) || any(!is.finite(measured[[length]]))) {
    cli_abort("Lengths in `measured` must be finite numbers.")
  }
  if (!is.logical(measured[[aged]]) || anyNA(measured[[aged]])) {
    cli_abort("{.val {aged}} in `measured` must be `TRUE` or `FALSE`.")
  }

  units <- unique(measured[[unit]])
  unit_f <- factor(match(measured[[unit]], units), levels = seq_along(units))
  bin_f <- factor(.lsas_bin(measured[[length]], breaks),
    levels = seq_len(length(breaks) - 1L))
  is_aged <- measured[[aged]]
  n_measured <- table(unit_f, bin_f)
  n_aged <- table(unit_f[is_aged], bin_f[is_aged])
  # Empty bins have no aged fish, so their fraction is zero.
  frac <- matrix(n_aged / pmax(n_measured, 1), nrow = length(units))
  structure(list(unit = unit, units = units, breaks = breaks, frac = frac),
    class = "sdmTMB_length_strata")
}

# Length bin of each length; lengths beyond the edges fall in the open-ended
# end bins.
.lsas_bin <- function(length, breaks) {
  n_bins <- length(breaks) - 1L
  pmin(pmax(findInterval(length, breaks), 1L), n_bins)
}

# Data for the RTMB objective: the unit-by-bin matrix of sampling fractions,
# the bin edges, and each fish's unit.
.lsas_tmb_data <- function(x, data, y, family_spec, backend) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, "sdmTMB_length_strata")) {
    cli_abort("`length_stratified` must be created with `length_strata()`.")
  }
  if (backend != "rtmb") {
    cli_abort("`length_stratified` needs `sdmTMBcontrol(backend = \"rtmb\")`.")
  }
  family <- family_spec$family$family
  if (length(family) != 1L || !family %in% c("gaussian", "lognormal")) {
    cli_abort("`length_stratified` needs the `gaussian()` or `lognormal()` family.")
  }
  if (!x$unit %in% names(data)) {
    cli_abort("Column {.val {x$unit}} is missing from `data`.")
  }
  if (anyNA(data[[x$unit]])) cli_abort("Sampling units in `data` can't have missing values.")
  unit_i <- match(data[[x$unit]], x$units)
  if (anyNA(unit_i)) {
    cli_abort("Every sampling unit in `data` must occur in `measured`.")
  }
  if (any(x$frac[cbind(unit_i, .lsas_bin(y, x$breaks))] == 0)) {
    cli_abort("Some fish in `data` fall in a unit and length bin with no aged fish in `measured`.")
  }
  list(frac = x$frac, breaks = x$breaks, unit_i = unit_i)
}

# Log probability that each fitted fish `i` was selected for ageing given its
# age: the log sum over length bins of the bin's sampling fraction times the
# bin's probability under the fish's length distribution `s`. Both families
# are normal on the modelled scale, so bin edges become z-scores.
rtmb_lsas_log_selection <- function(s, family, lsas, i) {
  "[<-" <- RTMB::ADoverload("[<-")
  frac <- lsas$frac[lsas$unit_i[i], , drop = FALSE]
  edges <- c(-Inf, lsas$breaks[-c(1L, length(lsas$breaks))], Inf)
  centre <- s$mu
  if (family == "lognormal") {
    edges <- log(pmax(edges, 0)) # bins at or below zero have no mass
    centre <- log(s$mu) - s$phi^2 / 2
  }
  selected <- rep(-Inf, length(i))
  for (k in seq_len(ncol(frac))) {
    rows <- which(frac[, k] > 0)
    if (!length(rows) || edges[k + 1L] == -Inf) next
    log_p <- rtmb_lsas_log_interval(edges[k], edges[k + 1L],
      rtmb_pick(centre, rows), rtmb_pick(s$phi, rows))
    selected[rows] <- RTMB::logspace_add(selected[rows], log(frac[rows, k]) + log_p)
  }
  selected
}

# Log probability that a normal variable with mean `centre` and SD `sd` falls
# between data edges `lower` and `upper`, either of which can be infinite.
# A finite interval's probability is symmetric in its standardized centre, so
# reflecting it into the lower tail keeps the difference of log CDFs accurate
# far into either tail without branching on parameters.
rtmb_lsas_log_interval <- function(lower, upper, centre, sd) {
  if (is.infinite(lower) && is.infinite(upper)) return(0 * centre)
  if (is.infinite(lower)) return(RTMB::pnorm((upper - centre) / sd, log.p = TRUE))
  if (is.infinite(upper)) return(RTMB::pnorm((centre - lower) / sd, log.p = TRUE))
  mid <- abs((lower + upper) / 2 - centre) / sd
  half_width <- (upper - lower) / (2 * sd)
  RTMB::logspace_sub(RTMB::pnorm(half_width - mid, log.p = TRUE),
    RTMB::pnorm(-half_width - mid, log.p = TRUE))
}

.check_lsas_diagnostics <- function(object, caller) {
  if (!is.null(object$tmb_data$lsas)) {
    cli_abort(c(
      "{caller} is not yet supported for `length_stratified` fits.",
      "i" = "Residuals and response simulations must account for selection for ageing. Population predictions are available with `predict()`."
    ))
  }
}
