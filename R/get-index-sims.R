#' Calculate a population index via simulation from the joint precision matrix
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Calculate a population index via simulation from the joint precision matrix.
#' Compared to [get_index()], this version can be faster if bias correction was
#' turned on in [get_index()] while being approximately equivalent. **This is an
#' experimental function.** This function usually works reasonably well, but we
#' make no guarantees. It is recommended to use [get_index()] with `bias_correct
#' = TRUE` for final inference.
#'
#' @details Can also be used to produce an index from a model fit with
#'   \pkg{tmbstan}.
#'
#' @details This function does nothing more than summarize and reshape the
#'   matrix of simulation draws into a data frame.
#'
#' @param obj [predict.sdmTMB()] output with `nsim > 0`.
#' @param level The confidence level.
#' @param return_sims Logical. Return simulation draws? The default (`FALSE`) is
#'   a quantile summary of those simulation draws.
#' @param area A vector of grid cell/polygon areas for each year-grid cell (row
#'   of data) in `obj`. Adjust this if cells are not of unit area or not all
#'   the same area (e.g., some cells are partially over land/water). Note that
#'   the area vector is added as `log(area)` to the raw values in `obj`. In
#'   other words, the function assumes a log link, which typically makes sense.
#' @param est_function Function to summarize the estimate (the expected value).
#'   `mean()` would be an alternative to `median()`.
#' @param area_function Function to apply area weighting.
#'   Assuming a log link, the default `function(x, area) x + log(area)` makes sense.
#'   If in natural space, `function(x, area) x * area` makes sense.
#' @param agg_function Function to aggregate samples within each time slice.
#'   Assuming a log link, the default `function(x) sum(exp(x))` makes sense.
#'   If in natural space, `function(x) sum(x)` makes sense.
#'
#' @seealso [get_index()]
#'
#' @return
#' A data frame. If `return_sims = FALSE`:
#'
#' * the name of the time column supplied to [sdmTMB()] (e.g., `year`)
#' * `est`: estimate
#' * `lwr`: lower confidence interval value
#' * `upr`: upper confidence interval value
#' * `log_est`: log estimate
#' * `se`: standard error on the log estimate
#'
#' If `return_sims = TRUE`, a long-format data frame of index draws:
#'
#' * the name of the time column supplied to [sdmTMB()] (e.g., `year`)
#' * `.value`: sample value
#' * `.iteration`: sample number
#'
#' @export
#' @examples
#' \donttest{
#' m <- sdmTMB(density ~ 0 + as.factor(year),
#'   data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie(link = "log"),
#'   time = "year"
#' )
#' qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))
#' p <- predict(m, newdata = qcs_grid_2011, nsim = 100)
#' x <- get_index_sims(p)
#' x_sims <- get_index_sims(p, return_sims = TRUE)
#'
#' if (require("ggplot2", quietly = TRUE)) {
#'   ggplot(x, aes(year, est, ymin = lwr, ymax = upr)) +
#'     geom_line() +
#'     geom_ribbon(alpha = 0.4)
#'   ggplot(x_sims, aes(as.factor(year), .value)) +
#'     geom_violin()
#' }
#'
#' # Demo custom functions if working in natural space:
#' ind <- get_index_sims(
#'   exp(p),
#'   agg_function = function(x) sum(x),
#'   area_function = function(x, area) x * area
#' )
#' }
get_index_sims <- function(obj,
                           level = 0.95,
                           return_sims = FALSE,
                           area = rep(1, nrow(obj)),
                           est_function = stats::median,
                           area_function = function(x, area) x + log(area),
                           agg_function = function(x) sum(exp(x))) {
  assert_that(is.matrix(obj), !is.null(attr(obj, "time")),
    msg = paste0("`obj` should be matrix output from `predict.sdmTMB()` ",
      "with the `nsim > 0` or a matrix with a `time` attribute."))
  assert_that(is.logical(return_sims))
  assert_that(is.function(est_function))
  assert_that(is.function(agg_function))
  assert_that(level > 0 && level < 1)
  assert_that(length(area) == nrow(obj))
  assert_that(sum(is.na(area)) == 0L)
  assert_that(all(area >= 0))

  cli_inform(c("We generally recommend using `get_index(..., bias_correct = TRUE)`",
    "rather than `get_index_sims()`."))

  if(length(attributes(obj))>3) {
  if (!(attr(obj,"link") == "log")) {
    cli_warn(c("Default `agg_function` and `area_function` apply to ",
      "values provided that are in log space. This matrix may be of ",
      "`type` = `response` or in a different link space."))
  }
  } else {
    cli_warn(c("Your predictions appear to be lacking a link attribute, so be ",
               "sure to check that the area_ and agg_ functions are",
               "appropriate."))
  }
  .time_attr <- attr(obj, "time")
  obj <- apply(obj, 2L, function(x) area_function(x, area))

  .t <- as.numeric(rownames(obj))
  yrs <- sort(unique(.t))
  yr_indexes <- lapply(yrs, function(x) which(.t %in% x))
  out1 <- lapply(yr_indexes, function(x) {
    apply(obj[x, , drop = FALSE], 2L, agg_function)
  })
  if (return_sims) {
    out2 <- lapply(seq_along(out1), function(i) {
      ret <- data.frame(
        .time = yrs[i], .value = out1[[i]],
        .iteration = seq_along(out1[[i]])
      )
      stats::setNames(ret, c(.time_attr, ".value", ".iteration"))
    })
    return(do.call("rbind", out2))
  } else {
    out <- lapply(out1, function(x) {
      data.frame(
        est = est_function(x),
        lwr = stats::quantile(x, probs = (1 - level) / 2),
        upr = stats::quantile(x, probs = 1 - (1 - level) / 2),
        se = stats::sd(log(x)),
        log_est = mean(log(x))
      )
    })
    out <- do.call("rbind", out)
    out[[.time_attr]] <- yrs
    out <- out[, c(
      .time_attr, "est", "lwr", "upr", "log_est", "se"
    ), drop = FALSE]
    return(`row.names<-`(out, NULL))
  }
}
