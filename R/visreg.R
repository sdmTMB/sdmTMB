#' Plot sdmTMB models with the \pkg{visreg} package
#'
#' sdmTMB models fit with regular (non-delta) families can be passed to
#' [visreg::visreg()] or [visreg::visreg2d()] directly. Examples are shown
#' below. Delta models can use the helper functions `visreg_delta()` or
#' `visreg2d_delta()` described here.
#'
#' @param object Fit from [sdmTMB()].
#' @param model First or second delta-model component.
#' @param ... Arguments passed to [visreg::visreg()] or
#'   [visreg::visreg2d()].
#'
#' @details
#' Note that the residuals are currently randomized quantile residuals,
#' *not* the deviance residuals usually used for GLMs with \pkg{visreg}.
#'
#' @return
#' A plot from the visreg package. If `plot = FALSE`, the plotted data are
#' returned invisibly. This is useful if you want to make your own plot.
#'
#' @export
#' @rdname visreg_delta
#'
#' @examples
#' if (require("ggplot2", quietly = TRUE) &&
#'   require("visreg", quietly = TRUE)) {
#'
#' \donttest{
#'   fit <- sdmTMB(
#'     density ~ s(depth_scaled),
#'     data = pcod_2011,
#'     spatial = "off",
#'     family = tweedie()
#'   )
#'   visreg::visreg(fit, xvar = "depth_scaled")
#'
#'   visreg::visreg(fit, xvar = "depth_scaled", scale = "response")
#'   v <- visreg::visreg(fit, xvar = "depth_scaled")
#'   head(v$fit)
#'   # now use ggplot2 etc. if desired
#'
#'   # Delta model example:
#'   fit_dg <- sdmTMB(
#'     density ~ s(depth_scaled, year, k = 8),
#'     data = pcod_2011, mesh = pcod_mesh_2011,
#'     spatial = "off",
#'     family = delta_gamma()
#'   )
#'   visreg_delta(fit_dg, xvar = "depth_scaled", model = 1, gg = TRUE)
#'   visreg_delta(fit_dg, xvar = "depth_scaled", model = 2, gg = TRUE)
#'   visreg_delta(fit_dg,
#'     xvar = "depth_scaled", model = 1,
#'     scale = "response", gg = TRUE
#'   )
#'   visreg_delta(fit_dg,
#'     xvar = "depth_scaled", model = 2,
#'     scale = "response"
#'   )
#'   visreg_delta(fit_dg,
#'     xvar = "depth_scaled", model = 2,
#'     scale = "response", gg = TRUE, rug = FALSE
#'   )
#'   visreg2d_delta(fit_dg,
#'     xvar = "depth_scaled", yvar = "year",
#'     model = 2, scale = "response"
#'   )
#'   visreg2d_delta(fit_dg,
#'     xvar = "depth_scaled", yvar = "year",
#'     model = 1, scale = "response", plot.type = "persp"
#'   )
#'   visreg2d_delta(fit_dg,
#'     xvar = "depth_scaled", yvar = "year",
#'     model = 2, scale = "response", plot.type = "gg"
#'   )
#'   }
#' }
visreg_delta <- function(object, ..., model = c(1, 2)) {
  object$visreg_model <- check_model_arg(model)
  dat <- object$data[!is.na(object$tmb_data$y_i[, model]), , drop = FALSE]
  visreg::visreg(fit = object, data = dat, ...)
}

#' @export
#' @rdname visreg_delta
visreg2d_delta <- function(object, ..., model = c(1, 2)) {
  object$visreg_model <- check_model_arg(model)
  dat <- object$data[!is.na(object$tmb_data$y_i[, model]), , drop = FALSE]
  visreg::visreg2d(fit = object, data = dat, ...)
}

check_model_arg <- function(model) {
  assert_that(is.numeric(model[[1]]))
  model <- as.integer(model[[1]])
  assert_that(model %in% c(1L, 2L))
  model
}
