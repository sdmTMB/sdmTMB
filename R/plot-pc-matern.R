#' Plot PC Matérn priors
#'
#' @param range_gt A value one expects the spatial or spatiotemporal range is
#'   **g**reater **t**han with `1 - range_prob` probability.
#' @param sigma_lt A value one expects the spatial or spatiotemporal marginal
#'   standard deviation (`sigma_O` or `sigma_E` internally) is **l**ess **t**han
#'   with `1 - sigma_prob` probability.
#' @param range_prob Probability. See description for `range_gt`.
#' @param sigma_prob Probability. See description for `sigma_lt`.
#' @param range_lims Plot range variable limits.
#' @param sigma_lims Plot sigma variable limits.
#' @param plot Logical controlling whether plot is drawn (defaults to `TRUE`).
#'
#' @seealso
#' [pc_matern()]
#'
#' @return
#' A plot from [graphics::image()].
#' Invisibly returns the underlying matrix data. The rows are the sigmas. The
#' columns are the ranges. Column and row names are provided.
#'
#' @export
#' @examples
#' plot_pc_matern(range_gt = 5, sigma_lt = 1)
#' plot_pc_matern(range_gt = 5, sigma_lt = 10)
#' plot_pc_matern(range_gt = 5, sigma_lt = 1, sigma_prob = 0.2)
#' plot_pc_matern(range_gt = 5, sigma_lt = 1, range_prob = 0.2)
plot_pc_matern <- function(range_gt,
                           sigma_lt,
                           range_prob = 0.05,
                           sigma_prob = 0.05,
                           range_lims = c(range_gt * 0.1, range_gt * 10),
                           sigma_lims = c(0, sigma_lt * 2),
                           plot = TRUE) {

  assert_that(range_prob > 0 & range_prob < 1)
  assert_that(sigma_prob > 0 & sigma_prob < 1)
  assert_that(range_gt > 0)
  assert_that(sigma_lt > 0)
  ranges <- seq(range_lims[1], range_lims[2], length.out = 200)
  sigmas <- seq(sigma_lims[1], sigma_lims[2], length.out = 201)
  out <- sapply(ranges, function(.range) {
    sapply(sigmas, function(.sigma) {
      inla_docs_log(.range, .sigma,
        matern_range = range_gt,
        matern_SD = sigma_lt,
        range_prob = range_prob,
        SD_prob = sigma_prob
      )
    })
  })
  if (isTRUE(plot)) {
    graphics::image(
      ranges, sigmas,
      z = t(exp(out)),
      col = .viridis100,
      xlab = "Range", ylab = "Sigma"
    )
    graphics::abline(v = range_gt, col = "white")
    graphics::abline(h = sigma_lt, col = "white")
  }
  rownames(out) <- sigmas
  colnames(out) <- ranges
  invisible(out)
}

inla_docs_log <- function(range = 1, sigma = 0.6, matern_range = 5,
                          range_prob = 0.05, matern_SD = 2, SD_prob = 0.05) {
  d <- 2
  dhalf <- d / 2
  lam1 <- -log(range_prob) * matern_range^dhalf
  lam2 <- -log(SD_prob) / matern_SD
  log(dhalf) + log(lam1) + log(range^(-1 - dhalf)) - lam1 * range^-dhalf +
    log(lam2) - lam2 * sigma
}

.viridis100 <-
  c(
    "#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF",
    "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF",
    "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF",
    "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF",
    "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF",
    "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF",
    "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF",
    "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF",
    "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF",
    "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF",
    "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF",
    "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF",
    "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF",
    "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF",
    "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF",
    "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF",
    "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF",
    "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF",
    "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF",
    "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF"
  )
