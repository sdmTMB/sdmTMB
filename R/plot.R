#' Plot anisotropy from an sdmTMB model
#'
#' Anisotropy is when spatial correlation is directionally dependent. In
#' [sdmTMB()], the default spatial correlation is isotropic, but anisotropy can
#' be enabled with `anisotropy = TRUE`. These plotting functions help visualize
#' that estimated anisotropy.
#'
#' @param object An object from [sdmTMB()].
#' @param return_data Logical. Return a data frame? `plot_anisotropy()` only.
#' @param model Which model if a delta model (only for `plot_anisotropy2()`;
#'   `plot_anisotropy()` always plots both).
#'
#' @return
#' `plot_anisotropy()`: One or more ellipses illustrating the estimated
#' anisotropy. The ellipses are centered at coordinates of zero in the space of
#' the X-Y coordinates being modeled. The ellipses show the spatial and/or
#' spatiotemporal range (distance at which correlation is effectively
#' independent) in any direction from zero. Uses \pkg{ggplot2}. If anisotropy
#' was turned off when fitting the model, `NULL` is returned instead of a
#' \pkg{ggplot2} object.
#'
#' `plot_anisotropy2()`: A plot of eigenvectors illustrating the estimated
#' anisotropy. A list of the plotted data is invisibly returned. Uses base
#' graphics. If anisotropy was turned off when fitting the model, `NULL` is
#' returned instead of a plot object.
#' @references Code adapted from VAST R package
#' @importFrom rlang .data
#' @examplesIf ggplot2_installed()
#' mesh <- make_mesh(pcod_2011, c("X", "Y"), n_knots = 80, type = "kmeans")
#' fit <- sdmTMB(
#'   data = pcod_2011,
#'   formula = density ~ 1,
#'   mesh = mesh,
#'   family = tweedie(),
#'   share_range = FALSE,
#'   anisotropy = TRUE #<
#' )
#' plot_anisotropy(fit)
#' plot_anisotropy2(fit)

#' @export
#' @rdname plot_anisotropy
plot_anisotropy <- function(object, return_data = FALSE) {
  stopifnot(inherits(object, "sdmTMB"))
  if (!check_for_H(object)) return(NULL)
  delta <- .is_regular_delta_model(object)

  # Calculate anisotropy components for model 1
  comp1 <- calculate_anisotropy_components(object, m = 1)
  eig <- comp1$eig
  maj1_s <- comp1$maj_s
  min1_s <- comp1$min_s
  maj1_st <- comp1$maj_st
  min1_st <- comp1$min_st

  # Calculate anisotropy components for model 2 if delta
  if (delta) {
    comp2 <- calculate_anisotropy_components(object, m = 2)
    eig2 <- comp2$eig
    maj2_s <- comp2$maj_s
    min2_s <- comp2$min_s
    maj2_st <- comp2$maj_st
    min2_st <- comp2$min_st
  }

  rss <- function(V) sqrt(sum(V[1]^2 + V[2]^2))
  get_angle <- function(m) {
    a <- -1 * (atan(m[1] / m[2]) / (2 * pi) * 360 - 90)
    a * (pi / 180)
  }

  angle1_s <- get_angle(maj1_s)
  angle1_st <- get_angle(maj1_st)

  if (delta) {
    angle2_s <- get_angle(maj2_s)
    angle2_st <- get_angle(maj2_st)
    dat <- data.frame(
      angle = c(angle1_s, angle1_st, angle2_s, angle2_st),
      a = c(rss(maj1_s), rss(maj1_st), rss(maj2_s), rss(maj2_st)),
      b = c(rss(min1_s), rss(min1_st), rss(min2_s), rss(min2_st)),
      maj1 = c(maj1_s, maj1_st, maj2_s, maj2_st),
      min1 = c(min1_s, min1_st, min2_s, min2_st),
      model = rep(object$family$family, each = 2L),
      model_num  = rep(seq(1L, 2L), each = 2L),
      random_field = rep(c("spatial", "spatiotemporal"), 2L),
      stringsAsFactors = FALSE
    )
    dat$model <- factor(dat$model, levels = object$family$family)
    for (i in seq(1L, 2L)) {
      if (object$spatiotemporal[i] == "off") {
        x <- dat$random_field == "spatiotemporal" & dat$model_num == i
        dat <- dat[!x, , drop = FALSE]
      }
    }
    for (i in seq(1L, 2L)) {
      if (object$spatial[i] == "off") {
        x <- dat$random_field == "spatial" & dat$model_num == i
        dat <- dat[!x, , drop = FALSE]
      }
    }
  } else {
    dat <- data.frame(
      angle = c(angle1_s, angle1_st),
      a = c(rss(maj1_s), rss(maj1_st)),
      b = c(rss(min1_s), rss(min1_st)),
      maj1 = c(maj1_s, maj1_st),
      min1 = c(min1_s, min1_st),
      model = object$family$family,
      random_field = rep(c("spatial", "spatiotemporal"), 1L),
      stringsAsFactors = FALSE
    )
    if (object$spatiotemporal == "off") {
      x <- dat$random_field == "spatiotemporal"
      dat <- dat[!x, , drop = FALSE]
    }
    if (object$spatial == "off") {
      x <- dat$random_field == "spatial"
      dat <- dat[!x, , drop = FALSE]
    }
  }

  if (return_data) return(dat)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli_abort("ggplot2 must be installed to use this function.")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    cli_abort("ggforce must be installed to use this function.")
  }
  g <- ggplot2::ggplot(dat,
    ggplot2::aes(
      x0 = 0, y0 = 0,
      a = .data$a, b = .data$b,
      angle = .data$angle,
      colour = `if`(delta, .data$model, NULL),
      linetype = .data$random_field
    )
  ) +
    ggforce::geom_ellipse() +
    ggplot2::coord_fixed() +
    ggplot2::labs(linetype = "Random field", colour = "Model",
      x = object$spde$xy_cols[1], y = object$spde$xy_cols[2]) +
    ggplot2::scale_colour_brewer(palette = "Dark2")
  g
}

#' @export
#' @rdname plot_anisotropy
plot_anisotropy2 <- function(object, model = 1) {
  stopifnot(inherits(object, "sdmTMB"))
  if (!check_for_H(object)) return(NULL)
  report <- object$tmb_obj$report(object$tmb_obj$env$last.par.best)
  if (model == 1) eig <- eigen(report$H)
  if (model == 2) eig <- eigen(report$H2)
  dat <- data.frame(
    x0 = c(0, 0),
    y0 = c(0, 0),
    x1 = eig$vectors[1, , drop = TRUE] * eig$values,
    y1 = eig$vectors[2, , drop = TRUE] * eig$values
  )
  plot(0,
    xlim = range(c(dat$x0, dat$x1)),
    ylim = range(c(dat$y0, dat$y1)),
    type = "n", asp = 1, xlab = "", ylab = ""
  )
  graphics::arrows(dat$x0, dat$y0, dat$x1, dat$y1)
  invisible(list(eig = eig, dat = dat, H = report$H))
}

#' Plot a smooth term from an sdmTMB model
#'
#' **Deprecated: use `visreg::visreg()`. See [visreg_delta()] for examples.**
#'
#' @param object An [sdmTMB()] model.
#' @param select The smoother term to plot.
#' @param n The number of equally spaced points to evaluate the smoother along.
#' @param level The confidence level.
#' @param ggplot Logical: use the \pkg{ggplot2} package?
#' @param rug Logical: add rug lines along the lower axis?
#' @param return_data Logical: return the predicted data instead of making a plot?
#' @export
#'
#' @details
#' Note:
#' * Any numeric predictor is set to its mean
#' * Any factor predictor is set to its first-level value
#' * The time element (if present) is set to its minimum value
#' * The x and y coordinates are set to their mean values
#'
#' @return
#' A plot of a smoother term.
#' @examples
#'   d <- subset(pcod, year >= 2000 & density > 0)
#'   pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 30)
#'   m <- sdmTMB(
#'     data = d,
#'     formula = log(density) ~ s(depth_scaled) + s(year, k = 5),
#'     mesh = pcod_spde
#'   )
#'   plot_smooth(m)
plot_smooth <- function(object, select = 1, n = 100, level = 0.95,
                        ggplot = FALSE, rug = TRUE, return_data = FALSE) {
  msg <- c(
    "This function may be deprecated.",
    "Consider using `visreg::visreg()` or `visreg_delta()`.",
    "See ?visreg_delta() for examples."
  )
  cli_inform(msg)
  se <- TRUE
  if (isTRUE(object$delta))
    cli_abort("This function doesn't work with delta models yet")

  assert_that(inherits(object, "sdmTMB"))
  assert_that(is.logical(ggplot))
  assert_that(is.logical(return_data))
  assert_that(is.logical(se))
  assert_that(is.numeric(n))
  assert_that(is.numeric(level))
  assert_that(length(level) == 1L)
  assert_that(length(select) == 1L)
  assert_that(length(n) == 1L)
  assert_that(is.numeric(select))
  assert_that(level > 0 & level < 1)
  assert_that(n < 500)

  if (ggplot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli_abort("ggplot2 not installed")
    }
  }

  sm <- parse_smoothers(object$formula[[1]], object$data)
  sm_names <- unlist(lapply(sm$Zs, function(x) attr(x, "s.label")))
  sm_names <- gsub("\\)$", "", gsub("s\\(", "", sm_names))

  fe_names <- colnames(object$tmb_data$X_ij)
  fe_names <- fe_names[!fe_names == "offset"] # FIXME
  fe_names <- fe_names[!fe_names == "(Intercept)"]

  all_names <- c(sm_names, fe_names)
  if (select > length(sm_names)) {
    cli_abort("`select` is greater than the number of smooths")
  }
  sel_name <- sm_names[select]
  non_select_names <- all_names[!all_names %in% sel_name]

  x <- object$data[[sel_name]]
  nd <- data.frame(x = seq(min(x), max(x), length.out = n))
  names(nd)[1] <- sel_name

  dat <- object$data
  .t <- terms(object$formula[[1]])
  .t <- labels(.t)
  checks <- c("^as\\.factor\\(", "^factor\\(")
  for (ch in checks) {
    if (any(grepl(ch, .t))) { # any factors from formula? if so, explicitly switch class
      ft <- grep(ch, .t)
      for (i in ft) {
        x <- gsub(ch, "", .t[i])
        x <- gsub("\\)$", "", x)
        dat[[x]] <- as.factor(dat[[x]])
      }
    }
  }
  dat[, object$spde$xy_cols] <- NULL
  dat[[object$time]] <- NULL
  for (i in seq_len(ncol(dat))) {
    if (names(dat)[i] != sel_name) {
      if (is.factor(dat[, i, drop = TRUE])) {
        nd[[names(dat)[[i]]]] <- sort(dat[, i, drop = TRUE])[[1]] # TODO note!
      } else {
        nd[[names(dat)[[i]]]] <- mean(dat[, i, drop = TRUE], na.rm = TRUE) # TODO note!
      }
    }
  }
  nd[object$time] <- min(object$data[[object$time]], na.rm = TRUE) # TODO note!
  nd[[object$spde$xy_cols[1]]] <- mean(object$data[[object$spde$xy_cols[1]]], na.rm = TRUE) # TODO note!
  nd[[object$spde$xy_cols[2]]] <- mean(object$data[[object$spde$xy_cols[2]]], na.rm = TRUE) # TODO note!

  p <- predict(object, newdata = nd, se_fit = se, re_form = NA)
  if (return_data) {
    return(p)
  }
  inv <- object$family$linkinv
  qv <- stats::qnorm(1 - (1 - level) / 2)

  if (!ggplot) {
    if (se) {
      lwr <- inv(p$est - qv * p$est_se)
      upr <- inv(p$est + qv * p$est_se)
      ylim <- range(c(lwr, upr))
    } else {
      ylim <- range(p$est)
    }
    plot(nd[[sel_name]], inv(p$est),
      type = "l", ylim = ylim,
      xlab = sel_name, ylab = paste0("s(", sel_name, ")")
    )
    if (se) {
      graphics::lines(nd[[sel_name]], lwr, lty = 2)
      graphics::lines(nd[[sel_name]], upr, lty = 2)
    }
    if (rug) rug(object$data[[sel_name]])
  } else {
    g <- ggplot2::ggplot(p, ggplot2::aes(.data[[sel_name]], inv(.data$est),
      ymin = inv(.data$est - qv * .data$est_se), ymax = inv(.data$est + qv * .data$est_se)
    )) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(alpha = 0.4) +
      ggplot2::labs(x = sel_name, y = paste0("s(", sel_name, ")"))
    if (rug) {
      g <- g +
        ggplot2::geom_rug(
          data = object$data, mapping = ggplot2::aes(x = .data[[sel_name]]),
          sides = "b", inherit.aes = FALSE, alpha = 0.3
        )
    }
    return(g)
  }
}

# Calculate anisotropic eigenvectors and ranges for a given model
# Returns list with: eig, range_s, range_st, maj_s, min_s, maj_st, min_st
calculate_anisotropy_components <- function(x, m = 1L) {
  # Get report and extract H matrices
  report <- x$tmb_obj$report(x$tmb_obj$env$last.par.best)
  delta <- .is_regular_delta_model(x)

  # Extract range values from sd_report
  est_rep <- as.list(x$sd_report, "Estimate", report = TRUE)
  range_values <- as.numeric(est_rep$range[, m])
  range_s <- range_values[1]
  range_st <- if (length(range_values) > 1) range_values[2] else range_values[1]

  # Get eigenvalues/vectors from H matrix (or H2 for delta model 2)
  H <- if (delta && m == 2) report$H2 else report$H
  eig <- eigen(H)

  # Calculate major and minor axis vectors for spatial field
  maj_s <- eig$vectors[, 1, drop = TRUE] * eig$values[1] * range_s
  min_s <- eig$vectors[, 2, drop = TRUE] * eig$values[2] * range_s

  # Calculate major and minor axis vectors for spatiotemporal field
  maj_st <- eig$vectors[, 1, drop = TRUE] * eig$values[1] * range_st
  min_st <- eig$vectors[, 2, drop = TRUE] * eig$values[2] * range_st

  list(
    eig = eig,
    range_s = range_s,
    range_st = range_st,
    maj_s = maj_s,
    min_s = min_s,
    maj_st = maj_st,
    min_st = min_st
  )
}

check_for_H <- function(obj) {
  H <- any(grepl(
    pattern = "ln_H_input",
    x = names(obj$sd_report$par.fixed),
    ignore.case = TRUE
  ))
  if (!H) {
    cli::cli_inform("`anisotropy = FALSE` in `sdmTMB()`; no anisotropy figure is available.")
    # FIXME in the future plot the isotropic covariance instead of NULL?
  }
  H
}
