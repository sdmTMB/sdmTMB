#' Construct an SPDE mesh for sdmTMB
#'
#' Construct an SPDE mesh for use with sdmTMB.
#'
#' @param data A data frame.
#' @param xy_cols A character vector of x and y column names contained in
#'   `data`. These should generally be in a projected CRS with distance-based units. For
#'   a helper function to convert to UTMs, see [add_utm_columns()].
#' @param type Method to create the mesh. Also see `mesh` argument to supply
#'   your own mesh.
#' @param cutoff An optional cutoff if type is `"cutoff"`. The minimum allowed
#'   triangle edge length.
#' @param n_knots The number of desired knots if `type` is not `"cutoff"`.
#' @param seed Random seed. Affects [stats::kmeans()] determination of knot
#'   locations if `type = "kmeans"`.
#' @param mesh An optional mesh created via \pkg{fmesher} instead of using the above
#'   convenience options.
#' @param fmesher_func Which \pkg{fmesher} function to use. Options include
#'   [fmesher::fm_rcdt_2d_inla()] and [fmesher::fm_mesh_2d_inla()], along with
#'   versions without the `_inla` suffix.
#' @param convex If specified, passed to [fmesher::fm_nonconvex_hull()].
#'   Distance to extend non-convex hull from data.
#' @param concave If specified, passed to [fmesher::fm_nonconvex_hull()].
#'   "Minimum allowed reentrant curvature". Defaults to `convex`.
#' @param ... Other arguments to pass to `fmesher_func`. Common arguments
#'   include `offset` and `max.edge`. See examples below.
#'
#' @return
#' `make_mesh()`: A list of class `sdmTMBmesh`. The element `mesh` is the output
#' from `fmesher_func` (default is [fmesher::fm_mesh_2d_inla()]). See
#' `mesh$mesh$n` for the number of vertices.
#'
#' @export
#'
#' @examples
#' # Extremely simple cutoff:
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 5, type = "cutoff")
#' plot(mesh)
#'
#' # Using a k-means algorithm to assign vertices:
#' mesh <- make_mesh(pcod, c("X", "Y"), n_knots = 50, type = "kmeans")
#' plot(mesh)
#'
#' \donttest{
#' # But, it's better to develop more tailored meshes:
#'
#' # Pass arguments via '...' to fmesher::fm_mesh_2d_inla():
#' mesh <- make_mesh(
#'   pcod, c("X", "Y"),
#'   fmesher_func = fmesher::fm_mesh_2d_inla,
#'   cutoff = 8, # minimum triangle edge length
#'   max.edge = c(20, 40), # inner and outer max triangle lengths
#'   offset = c(5, 40) # inner and outer border widths
#' )
#' plot(mesh)
#'
#' # Or define a mesh directly with fmesher (formerly in INLA):
#' inla_mesh <- fmesher::fm_mesh_2d_inla(
#'   loc = cbind(pcod$X, pcod$Y), # coordinates
#'   max.edge = c(25, 50), # max triangle edge length; inner and outer meshes
#'   offset = c(5, 25),  # inner and outer border widths
#'   cutoff = 5 # minimum triangle edge length
#' )
#' mesh <- make_mesh(pcod, c("X", "Y"), mesh = inla_mesh)
#' plot(mesh)
#' }

make_mesh <- function(data, xy_cols,
                      type = c("kmeans", "cutoff", "cutoff_search"),
                      cutoff,
                      n_knots,
                      seed = 42,
                      mesh = NULL,
                      fmesher_func = fmesher::fm_rcdt_2d_inla,
                      convex = NULL, concave = convex,
                      ...) {

  if (nrow(data) == 0L) {
    cli_abort("The data frame supplied to `data` has no rows of data.")
  }
  if (missing(xy_cols) || is.numeric(xy_cols) || is.numeric(data)) {
    msg <- paste0(
      "It looks like you are using an old format of make_mesh(). ",
      "The function now uses `data` and `xy_cols` arguments ",
      "to enable carrying through the x and y column names ",
      "to the predict function. Please update your code."
    )
    cli_abort(msg)
  }

  all_x_non_na <- sum(is.na(data[[xy_cols[[1]]]])) == 0L
  all_y_non_na <- sum(is.na(data[[xy_cols[[2]]]])) == 0L
  if (!all_x_non_na || !all_y_non_na) {
    msg <- c("Some coordinates in `xy_cols` were NA.", "
      Remove or fix these rows before proceeding.")
    cli_abort(msg)
  }

  if (max(data[[xy_cols[1]]] > 1e4) || max(data[[xy_cols[2]]] > 1e4)) {
    msg <- paste0(
      "The x or y column values are fairly large. ",
      "This can cause estimation problems since the spatial range ",
      "is dependent on the scale of the coordinates. ",
      "Consider scaling the x and y coordinates. ",
      "For example, try working in UTM km instead of UTM m by divided by 1000."
      )
    cli_warn(msg)
  }
  type <- match.arg(type)
  if (!missing(n_knots) && type == "cutoff") {
    msg <- paste0(
      "`n_knots` specified but `type = 'cutoff'. ",
      "The default mesh type now uses a `cutoff` ",
      "intead of the number of knots. You can obtain the previous default ",
      "mesh type with `type = 'kmeans'`."
      )
    cli_abort(msg)
  }

  if (!missing(cutoff) && missing(n_knots)) {
    type <- "cutoff"
  }
  if (missing(cutoff) && type == "cutoff" && is.null(mesh)) {
    cli_abort("You need to specify the `cutoff` argument.")
  }
  if (missing(n_knots) && type != "cutoff" && is.null(mesh)) {
    cli_abort("You need to specify the `n_knots` argument.")
  }
  loc_xy <- as.matrix(data[, xy_cols, drop = FALSE])
  loc_centers <- NA

  if (is.null(mesh)) {
    if (type == "kmeans") {
      if (n_knots >= nrow(loc_xy)) {
        cli_warn(
          "Reducing `n_knots` to be one less than the number of data points."
        )
        n_knots <- nrow(loc_xy) - 1
      }
      set.seed(seed)
      knots <- stats::kmeans(x = loc_xy, centers = n_knots)
      loc_centers <- knots$centers
      if (!is.null(convex) || !is.null(concave)) {
        nch <- fmesher::fm_nonconvex_hull(loc_xy, convex = convex, concave = concave)
        mesh <- fmesher_func(loc_centers, refine = list(), boundary = nch, extend = list(), ...)
      } else {
        mesh <- fmesher_func(loc_centers, refine = list(), extend = list(), ...)
      }
    } else if (type == "cutoff") {
      if (!is.null(convex) || !is.null(concave)) {
        nch <- fmesher::fm_nonconvex_hull(loc_xy, convex = convex, concave = concave)
        mesh <- fmesher_func(loc_centers, refine = list(), cutoff = cutoff, boundary = nch, extend = list(), ...)
      } else {
        mesh <- fmesher_func(loc_xy, refine = list(), cutoff = cutoff, extend = list(), ...)
      }
    } else {
      mesh <- binary_search_knots(loc_xy, n_knots = n_knots,
        refine = list(), fmesher_func = fmesher_func, ...)
    }
  } else {
    knots <- list()
    loc_centers <- NA
  }
  # spde_inla <- INLA::inla.spde2.matern(mesh)
  spde <- fmesher::fm_fem(mesh)
  # identical(spde2$c0, spde_inla$param.inla$M0)
  # identical(as.matrix(spde2$g1), as.matrix(spde_inla$param.inla$M1))
  # identical(as.matrix(spde2$g2), as.matrix(spde_inla$param.inla$M2))
  A <- fmesher::fm_basis(mesh, loc = loc_xy)

  fake_data <- data
  fake_data[["sdm_spatial_id"]] <- seq(1, nrow(fake_data))

  if (mesh$n > 1000L) {
    msg <- paste0(
      "This mesh has > 1000 vertices. Mesh complexity has the single largest ",
      "influence on fitting speed. Consider whether you require a mesh this ",
      "complex, especially for initial model exploration. ",
      "Check `your_mesh$mesh$n` to view the number of vertices."
    )
    cli_inform(msg)
  }

  structure(list(
    loc_xy = loc_xy, xy_cols = xy_cols, mesh = mesh, spde = spde,
    loc_centers = loc_centers, A_st = A,
    sdm_spatial_id = fake_data$sdm_spatial_id
  ), class = "sdmTMBmesh")}

binary_search_knots <- function(loc_xy,
                                n_knots,
                                min = 1e-4,
                                max = 1e4,
                                length = 1e6,
                                fmesher_func = fmesher::fm_mesh_2d_inla,
                                refine = list(), ...) {
  vec <- exp(seq(log(min), log(max), length.out = length))
  L <- 0
  R <- length(vec)
  realized_knots <- Inf
  while (L <= R) {
    m <- floor((L + R) / 2)
    # mesh <- INLA::inla.mesh.create(loc_xy, refine = refine, cutoff = vec[m])
    mesh <- fmesher_func(loc_xy, refine = refine, extend = list(), cutoff = vec[m], ...)
    realized_knots <- mesh$n
    pretty_cutoff <- sprintf("%.2f", round(vec[m], 2))
    cat("cutoff =", pretty_cutoff, "| knots =", realized_knots)
    if (realized_knots > n_knots) {
      L <- m + 1
      cat(" |", cli::symbol$arrow_down, "\n")
    } else if (realized_knots < n_knots) {
      R <- m - 1
      cat(" |", cli::symbol$arrow_up, "\n")
    } else {
      cat(" |", cli::symbol$tick, "\n")
      return(mesh)
    }
  }
  cat("cutoff =", pretty_cutoff, "| knots =", mesh$n, "\n")
  mesh
}

#' @param x Output from [make_mesh()].
#' @param ... Passed to [graphics::plot()].
#'
#' @importFrom graphics points
#' @return `plot.sdmTMBmesh()`: A plot of the mesh and data points. To make your
#'   own \pkg{ggplot2} version, pass `your_mesh$mesh` to `inlabru::gg()`.
#' @rdname make_mesh
#' @export
plot.sdmTMBmesh <- function(x, ...) {
    plot(x$mesh, main = NA, edge.color = "grey60", asp = 1, ...)
    points(x$loc_xy, pch = 21, cex = 0.3, col = "#00000080")
    points(x$loc_centers, pch = 20, col = "red")
  # }
}

# from TMB examples repository:
make_anisotropy_spde <- function(spde, anistropy = TRUE) {
  if (anistropy) {
    inla_mesh <- spde$mesh
    Dset <- 1:2
    inla_mesh <- spde$mesh
    # Triangle info
    TV <- inla_mesh$graph$tv # Triangle to vertex indexing
    V0 <- inla_mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
    V1 <- inla_mesh$loc[TV[, 2], Dset]
    V2 <- inla_mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1 # E = edge for each triangle
    E1 <- V0 - V2
    E2 <- V1 - V0
    # Calculate Areas
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    for (i in seq_len(length(Tri_Area))) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2

    ret <- list(
      n_s = spde$mesh$n,
      n_tri = nrow(TV),
      Tri_Area = Tri_Area,
      E0 = E0,
      E1 = E1,
      E2 = E2,
      TV = TV - 1,
      G0 = spde$spde$c0,
      G0_inv = as(Matrix::diag(1 / Matrix::diag(spde$spde$c0)), "TsparseMatrix")
    )
  } else {
    ret <- list(
      n_s = 0L, n_tri = 0L, Tri_Area = rep(0, 1), E0 = matrix(0, 1),
      E1 = matrix(0, 1), E2 = matrix(0, 1), TV = matrix(0, 1),
      G0 = Matrix::Matrix(0, 1, 1, doDiag = FALSE), G0_inv = Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    )
  }
  ret
}

make_barrier_spde <- function(spde) {
  if ("spde_barrier" %in% names(spde)) {
    C0 <- spde$spde_barrier$C[[1]]
    C1 <- spde$spde_barrier$C[[2]]
    D0 <- spde$spde_barrier$D[[1]]
    D1 <- spde$spde_barrier$D[[2]]
    .I <- spde$spde_barrier$I
  } else {
    C0 <- rep(1, 2)
    C1 <- rep(1, 2)
    D0 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    D1 <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    .I <- Matrix::Matrix(0, 1, 1, doDiag = FALSE)
  }
  list(C0 = C0, C1 = C1, D0 = D0, D1 = D1, I = .I)
}
