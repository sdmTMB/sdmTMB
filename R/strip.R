#' Shrink a fitted model before saving it
#'
#' Removes objects from a fitted model that sdmTMB rebuilds as needed, so the
#' model takes less space when saved with [saveRDS()]. The stripped model
#' works the same as one read back from disk: its objective function is taped
#' again the first time a method needs it.
#'
#' @param object A model fit with [sdmTMB()].
#' @param strip_joint_precision Also remove the joint precision matrix of the
#'   fixed and random effects? This usually gives the largest saving, but makes
#'   some methods slower. See Details.
#'
#' @details
#' This removes the handle to the in-memory RTMB tape, which doesn't survive
#' a save and reload and can be the largest part of a saved model, and the
#' cached sparse Cholesky factor of the random-effects Hessian, which sdmTMB
#' rebuilds cheaply when needed.
#'
#' The joint precision matrix is often the largest part of what remains, but
#' removing it (`strip_joint_precision = TRUE`) has a cost later:
#'
#' * Simulation-based methods (e.g., `predict(..., nsim = )`,
#'   [get_index_sims()], [spread_sims()], [gather_sims()], and [project()])
#'   recompute it with [RTMB::sdreport()] each time they are called. It isn't
#'   stored back in the model, so every call costs about as much as the
#'   `sdreport()` step of the original fit, which can take minutes for large
#'   spatiotemporal models.
#' * [get_index()], [get_cog()], and [get_eao()] can't use a shortcut based on
#'   the matrix and instead build a new objective function and run
#'   [RTMB::sdreport()] on it.
#'
#' Other methods (e.g., `predict()` without `nsim`, [tidy()], and
#' [residuals()]) are unaffected, as are models fit without the matrix (see
#' `get_joint_precision` in [sdmTMBcontrol()]).
#'
#' The original `object` is not changed.
#'
#' @returns The fitted model with the removed objects.
#' @export
#' @examples
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
#' fit <- sdmTMB(density ~ 1, data = pcod_2011, mesh = mesh,
#'   family = tweedie())
#' fit_small <- strip_sdmTMB(fit)
#' object.size(serialize(fit, NULL))
#' object.size(serialize(fit_small, NULL))
#' p <- predict(fit_small)
strip_sdmTMB <- function(object, strip_joint_precision = FALSE) {
  assert_that(inherits(object, "sdmTMB"))
  assert_that(isTRUE(strip_joint_precision) || isFALSE(strip_joint_precision))
  # the TMB object is an environment; a serialize round trip gives an
  # independent copy with null pointers, the same as reading it from disk
  obj <- unserialize(serialize(object$tmb_obj, NULL))
  env <- obj$env
  if (exists("L.created.by.newton", envir = env, inherits = FALSE)) {
    rm("L.created.by.newton", envir = env)
  }
  attr(env$ADFun, "rcpp") <- NULL
  object$tmb_obj <- obj
  if (strip_joint_precision) object$sd_report$jointPrecision <- NULL
  object
}
