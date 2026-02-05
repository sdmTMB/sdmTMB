#' Calculate range edges via simulation from the joint precision matrix
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Calculate range edges as density-weighted quantiles along a spatial axis.
#' Range edges are calculated as the positions along a user-supplied spatial
#' axis (e.g., latitude, coastal distance) where the cumulative proportion of
#' density equals specified quantiles (e.g., 0.01 and 0.99 for the lower and
#' upper 1% range edges). Uncertainty is calculated via simulation from the
#' joint precision matrix.
#'
#' @details This function implements a similar approach to VAST's range edge
#'   calculations, following methods from Fredston et al. (2021) and similar
#'   studies. The method:
#'   \enumerate{
#'     \item Orders spatial locations by position along the specified axis
#'     \item Calculates cumulative proportion of total density along that axis
#'     \item Finds positions where cumulative proportion equals target quantiles
#'     \item Uses simulation from the joint precision to quantify uncertainty
#'   }
#'
#'   To find the exact position where the cumulative proportion equals a target
#'   quantile, the function uses linear interpolation between adjacent grid
#'   points. This provides more accurate range edge estimates than selecting the
#'   closest grid point, especially on coarser grids or for extreme quantiles
#'   (e.g., 0.01, 0.99).
#'
#' @param obj [predict.sdmTMB()] output with `nsim > 0`. The prediction object
#'   should include predictions on a spatial grid that covers the area of interest.
#' @param axis Numeric vector of the same length as the prediction data,
#'   representing the spatial axis along which to calculate range edges
#'   (e.g., latitude, coastal distance values). This should align with
#'   the rows of the prediction matrix.
#' @param quantiles Numeric vector of quantiles to calculate. Default is
#'   `c(0.025, 0.975)` for lower and upper 1% range edges. Common alternatives
#'   include `c(0.01, 0.99)` for 1% edges or `c(0.05, 0.5, 0.95)` to
#'   include the median.
#' @param level The confidence level for uncertainty intervals.
#' @param return_sims Logical. Return simulation draws? The default (`FALSE`)
#'   returns a quantile summary of the simulation draws.
#'
#' @references
#' Fredston, A. L., Pinsky, M., Selden, R. L., Szuwalski, C., Thorson, J. T.,
#' Gaines, S. D., & Halpern, B. S. (2021). Range edges of North American
#' marine species are tracking temperature over decades. Global Change Biology,
#' 27(13), 3145-3156. \doi{10.1111/gcb.15614}
#'
#' @return
#' A data frame. If `return_sims = FALSE`:
#'
#' * name of time column (e.g., `year`) that was supplied to [sdmTMB()] time argument
#' * `quantile`: the quantile value (from `quantiles` argument)
#' * `est`: estimated range edge position
#' * `lwr`: lower confidence interval
#' * `upr`: upper confidence interval
#' * `se`: standard error
#'
#' If `return_sims = TRUE`, simulation draws from range edge positions in long format:
#'
#' * name of time column (e.g., `year`)
#' * `quantile`: the quantile value
#' * `.value`: simulated range edge position
#' * `.iteration`: simulation number
#'
#' @export
#' @examples
#' \donttest{
#' # Fit a spatiotemporal model
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 100)
#' m <- sdmTMB(
#'   density ~ 0 + as.factor(year),
#'   data = pcod, mesh = mesh, family = tweedie(link = "log"),
#'   time = "year", spatiotemporal = "iid", spatial = "on"
#' )
#'
#' # Create prediction grid
#' nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
#'
#' # Get predictions with simulations
#' p <- predict(m, newdata = nd, nsim = 100)
#'
#' # Calculate range edges along latitude (Y coordinate)
#' edges <- get_range_edge(p, axis = nd$Y)
#' edges
#'
#' # Plot range edges over time
#' if (require("ggplot2", quietly = TRUE)) {
#'   ggplot(edges, aes(year, est, colour = as.factor(quantile))) +
#'     geom_line() +
#'     geom_ribbon(aes(ymin = lwr, ymax = upr, fill = as.factor(quantile)),
#'       alpha = 0.2
#'     ) +
#'     labs(y = "Latitude", colour = "Quantile", fill = "Quantile")
#' }
#'
#' # Get simulation draws for further analysis
#' edges_sims <- get_range_edge(p, axis = nd$Y, return_sims = TRUE)
#' }
get_range_edge <- function(obj,
                           axis,
                           quantiles = c(0.025, 0.975),
                           level = 0.95,
                           return_sims = FALSE) {
  # Validation
  assert_that(is.matrix(obj),
    !is.null(attr(obj, "time")),
    msg = paste0(
      "`obj` should be matrix output from `predict.sdmTMB()` ",
      "with `nsim > 0`."
    )
  )
  assert_that(is.numeric(axis),
    msg = "`axis` should be a numeric vector."
  )
  assert_that(length(axis) == nrow(obj),
    msg = paste0(
      "`axis` should be the same length as the number of rows in `obj` (",
      nrow(obj), ")."
    )
  )
  assert_that(is.numeric(quantiles), all(quantiles > 0), all(quantiles < 1),
    msg = "`quantiles` should be numeric values between 0 and 1."
  )
  assert_that(is.logical(return_sims))
  assert_that(level > 0 && level < 1)

  # Extract attributes
  .time_attr <- attr(obj, "time")
  .link <- attr(obj, "link")

  if (is.null(.link)) {
    cli_warn(c(
      "No link attribute found in prediction object.",
      "Assuming log link. If this is incorrect, the range edges will be wrong.",
      "Consider re-running predict() with a newer version of sdmTMB."
    ))
    .link <- "log"
  }

  # Get time values from rownames
  .t <- as.numeric(rownames(obj))
  if (is.null(.t) || any(is.na(.t))) {
    cli_abort(c(
      "Could not extract time values from prediction matrix rownames.",
      "This may indicate the prediction object is malformed."
    ))
  }
  yrs <- sort(unique(.t))

  # Calculate range edges for each simulation
  n_sims <- ncol(obj)
  n_times <- length(yrs)
  n_quantiles <- length(quantiles)

  # Initialize results array: [time, quantile, simulation]
  edges_array <- array(NA_real_,
    dim = c(n_times, n_quantiles, n_sims),
    dimnames = list(
      time = yrs,
      quantile = quantiles,
      sim = seq_len(n_sims)
    )
  )

  # Pre-compute time indices and orderings once (not per-simulation)
  time_cache <- vector("list", n_times)
  for (t_idx in seq_along(yrs)) {
    t <- yrs[t_idx]
    time_idx <- which(.t == t)
    axis_t <- axis[time_idx]
    ord <- order(axis_t)
    time_cache[[t_idx]] <- list(
      time_idx = time_idx,
      axis_t = axis_t,
      ord = ord,
      axis_ordered = axis_t[ord]
    )
  }

  # Loop over simulations
  for (sim in seq_len(n_sims)) {
    # Get predictions for this simulation
    pred_sim <- obj[, sim]

    # Loop over time points
    for (t_idx in seq_along(yrs)) {
      # Use cached values instead of recomputing
      cache <- time_cache[[t_idx]]
      pred_t <- pred_sim[cache$time_idx]

      # Convert to natural space using inverse link
      dens_t <- apply_inverse_link(pred_t, .link)

      # Order by axis position (use cached ordering)
      axis_ordered <- cache$axis_ordered
      dens_ordered <- dens_t[cache$ord]

      # Calculate cumulative proportion of density
      total_dens <- sum(dens_ordered)
      if (total_dens > 0) {
        cum_dens <- cumsum(dens_ordered)
        cum_prop <- cum_dens / total_dens

        # Find range edge for each quantile
        for (q_idx in seq_along(quantiles)) {
          q <- quantiles[q_idx]

          # Find the position where cumulative proportion crosses the quantile
          # Use linear interpolation between adjacent points
          cross_idx <- which(cum_prop >= q)[1]

          if (is.na(cross_idx)) {
            # Quantile is beyond the range
            edges_array[t_idx, q_idx, sim] <- NA_real_
          } else if (cross_idx == 1) {
            # Quantile is at or before first point
            edges_array[t_idx, q_idx, sim] <- axis_ordered[1]
          } else {
            # Linear interpolation between cross_idx - 1 and cross_idx
            idx1 <- cross_idx - 1
            idx2 <- cross_idx

            prop1 <- cum_prop[idx1]
            prop2 <- cum_prop[idx2]
            axis1 <- axis_ordered[idx1]
            axis2 <- axis_ordered[idx2]

            # Interpolate
            weight <- (q - prop1) / (prop2 - prop1)
            edge_pos <- axis1 + weight * (axis2 - axis1)

            edges_array[t_idx, q_idx, sim] <- edge_pos
          }
        }
      } else {
        # No density, can't calculate edges
        edges_array[t_idx, , sim] <- NA_real_
      }
    }
  }

  # Format output using pre-allocation instead of growing lists
  if (return_sims) {
    # Return long format with all simulations
    result <- data.frame(
      time = rep(yrs, each = n_quantiles * n_sims),
      quantile = rep(rep(quantiles, each = n_sims), times = n_times),
      .value = as.vector(aperm(edges_array, c(3, 2, 1))),
      .iteration = rep(seq_len(n_sims), times = n_times * n_quantiles)
    )
    names(result)[1] <- .time_attr
    return(result)
  } else {
    # Compute summary statistics efficiently using array operations
    est_mat <- apply(edges_array, c(1, 2), stats::median, na.rm = TRUE)
    lwr_mat <- apply(edges_array, c(1, 2), stats::quantile,
      probs = (1 - level) / 2, na.rm = TRUE
    )
    upr_mat <- apply(edges_array, c(1, 2), stats::quantile,
      probs = 1 - (1 - level) / 2, na.rm = TRUE
    )
    se_mat <- apply(edges_array, c(1, 2), stats::sd, na.rm = TRUE)

    # Pre-allocate and fill (transpose to get correct row ordering)
    result <- data.frame(
      time = rep(yrs, each = n_quantiles),
      quantile = rep(quantiles, times = n_times),
      est = as.vector(t(est_mat)),
      lwr = as.vector(t(lwr_mat)),
      upr = as.vector(t(upr_mat)),
      se = as.vector(t(se_mat))
    )
    names(result)[1] <- .time_attr
    rownames(result) <- NULL
    result
  }
}

# Helper function to apply inverse link transformation
apply_inverse_link <- function(eta, link) {
  switch(link,
    "log" = exp(eta),
    "logit" = plogis(eta),
    "identity" = eta,
    "inverse" = 1 / eta,
    "cloglog" = 1 - exp(-exp(eta)),
    "response" = eta, # already in natural space
    {
      cli_warn(paste0("Unknown link '", link, "'. Returning untransformed values."))
      eta
    }
  )
}
