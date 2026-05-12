#' Create an areal spatial domain for SAR/CAR models
#'
#' Build a reusable areal domain object with adjacency weights and unit labels.
#' The returned object can be supplied to `mesh` in [sdmTMB()] in areal
#' SAR/CAR workflows.
#'
#' @param data A data frame used for membership validation against the domain.
#' @param spatial_domain A named `igraph` object or a square numeric matrix /
#'   sparse matrix with matching row and column names. Can also be an `sf`/`sfc`
#'   polygon object when `id_column` identifies areal units.
#' @param space_column Column name in `data` / `newdata` that identifies areal
#'   unit membership.
#' @param id_column Optional column name in an `sf` polygon `spatial_domain`
#'   containing areal unit IDs. If omitted for `sf` polygons, stable IDs are
#'   generated.
#' @param adjacency Polygon adjacency type for `sf` polygon input: `"rook"` for
#'   shared edges or `"queen"` for any touching boundary.
#'
#' @return A list with class `c("sdmTMBareal", "sdmTMBdomain")`.
#' @export
make_areal_domain <- function(data, spatial_domain, space_column,
                              id_column = NULL,
                              adjacency = c("rook", "queen")) {
  if (!inherits(data, "data.frame")) {
    cli::cli_abort("`data` must be a data frame.")
  }
  if (!is.character(space_column) || length(space_column) != 1L || is.na(space_column) || !nzchar(space_column)) {
    cli::cli_abort("`space_column` must be a single non-empty column name.")
  }
  adjacency <- match.arg(adjacency)

  from_igraph <- FALSE
  if (inherits(spatial_domain, "igraph")) {
    from_igraph <- TRUE
    if (!requireNamespace("igraph", quietly = TRUE)) {
      cli::cli_abort("`igraph` must be installed to use an igraph `spatial_domain` input.")
    }
    unit_names <- igraph::V(spatial_domain)$name
    if (is.null(unit_names) || anyNA(unit_names) || any(!nzchar(unit_names))) {
      cli::cli_abort(c(
        "Named igraph input requires non-empty vertex names.",
        "i" = "Set names with `igraph::V(graph)$name <- ...`."
      ))
    }
    if (anyDuplicated(unit_names)) {
      cli::cli_abort("igraph vertex names must be unique.")
    }
    W <- .as_general_dgC(igraph::as_adjacency_matrix(spatial_domain, sparse = TRUE))
    Matrix::diag(W) <- 0
    dimnames(W) <- list(as.character(unit_names), as.character(unit_names))
  } else if (inherits(spatial_domain, "sf") || inherits(spatial_domain, "sfc")) {
    W <- .sf_areal_adjacency_matrix(spatial_domain, id_column = id_column, adjacency = adjacency)
  } else {
    if (!(inherits(spatial_domain, "Matrix") || is.matrix(spatial_domain))) {
      cli::cli_abort("`spatial_domain` must be a named igraph, an sf/sfc polygon object, or a numeric matrix / sparse matrix.")
    }
    W <- if (inherits(spatial_domain, "Matrix")) {
      .as_general_dgC(spatial_domain)
    } else {
      if (!is.numeric(spatial_domain)) {
        cli::cli_abort("Matrix `spatial_domain` input must be numeric.")
      }
      .as_general_dgC(Matrix::Matrix(spatial_domain, sparse = TRUE))
    }
    if (!is.numeric(W@x)) {
      cli::cli_abort("Matrix `spatial_domain` input must be numeric.")
    }
  }

  .validate_areal_matrix_structure(W, from_igraph = from_igraph)
  unit_names <- as.character(rownames(W))

  .validate_areal_matrix_values(W, stage = "before normalization")
  W_raw <- methods::as(W, "dgCMatrix")
  W <- .normalize_areal_W(W_raw)
  .validate_areal_matrix_values(W, stage = "after normalization")
  .validate_areal_row_sums(W)
  .validate_areal_data_membership(
    data = data,
    unit_names = unit_names,
    space_column = space_column
  )

  structure(
    list(
      W = W,
      W_raw = W_raw,
      n_s = as.integer(length(unit_names)),
      unit_names = unit_names,
      space_column = space_column
    ),
    class = c("sdmTMBareal", "sdmTMBdomain")
  )
}

#' Create an areal grid domain by overlaying point data
#'
#' Build or use an `sf` polygon grid, overlay point observations onto grid
#' cells, construct grid-cell adjacency, and return labelled data plus an
#' `sdmTMBareal` domain for SAR/CAR models.
#'
#' @param data A data frame containing point coordinates.
#' @param xy_cols Character vector of length 2 naming coordinate columns in
#'   `data`.
#' @param spatial_domain Optional `sf`/`sfc` polygon object. If `cellsize` or
#'   `n` is supplied, this object is treated as a boundary from which a grid is
#'   generated and clipped. Otherwise it is treated as the grid itself.
#' @param cellsize Optional cell size passed to [sf::st_make_grid()].
#' @param n Optional grid dimensions passed to [sf::st_make_grid()]. Ignored if
#'   `cellsize` is supplied.
#' @param square Logical passed to [sf::st_make_grid()]. Use `FALSE` for
#'   hexagonal grids.
#' @param space_column Column name to add to returned data.
#' @param adjacency Polygon adjacency type: `"rook"` for shared edges or
#'   `"queen"` for any touching boundary.
#' @param crs Optional coordinate reference system for `data` coordinates,
#'   passed to [sf::st_as_sf()].
#'
#' @return A list with `data`, `grid`, and `domain` elements. `domain` can be
#'   supplied to [sdmTMB()] via the `mesh` argument.
#' @export
make_areal_grid <- function(data, xy_cols, spatial_domain = NULL,
                            cellsize = NULL, n = NULL, square = TRUE,
                            space_column = "grid_cell",
                            adjacency = c("rook", "queen"),
                            crs = NA) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg sf} must be installed to use `make_areal_grid()`.")
  }
  if (!inherits(data, "data.frame")) {
    cli::cli_abort("`data` must be a data frame.")
  }
  if (!is.character(xy_cols) || length(xy_cols) != 2L || anyNA(xy_cols) || any(!nzchar(xy_cols))) {
    cli::cli_abort("`xy_cols` must be a character vector of length 2.")
  }
  if (!all(xy_cols %in% names(data))) {
    cli::cli_abort("All `xy_cols` must be columns in `data`.")
  }
  if (!is.character(space_column) || length(space_column) != 1L || is.na(space_column) || !nzchar(space_column)) {
    cli::cli_abort("`space_column` must be a single non-empty column name.")
  }
  adjacency <- match.arg(adjacency)

  pts <- sf::st_as_sf(data, coords = xy_cols, crs = crs, remove = FALSE)
  if (is.null(spatial_domain)) {
    spatial_domain <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(pts)))
  }
  spatial_domain <- .as_sf_polygons(spatial_domain, arg = "spatial_domain")
  pts <- .align_sf_crs(pts, spatial_domain)

  make_grid <- !is.null(cellsize) || !is.null(n)
  if (make_grid) {
    grid_args <- list(x = spatial_domain, square = square)
    if (!is.null(cellsize)) {
      grid_args$cellsize <- cellsize
    } else {
      grid_args$n <- n
    }
    grid_geom <- do.call(sf::st_make_grid, grid_args)
    grid_geom <- sf::st_intersection(grid_geom, sf::st_union(spatial_domain))
    if (any(sf::st_geometry_type(grid_geom) == "GEOMETRYCOLLECTION")) {
      grid_geom <- sf::st_collection_extract(grid_geom, "POLYGON")
    }
    grid <- sf::st_sf(geometry = grid_geom)
  } else {
    grid <- spatial_domain
  }

  centers <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grid)))
  grid[[xy_cols[1L]]] <- centers[, 1L]
  grid[[xy_cols[2L]]] <- centers[, 2L]
  grid[[space_column]] <- sprintf("cell_%03d", seq_len(nrow(grid)))
  hits <- sf::st_intersects(pts, grid)
  if (any(lengths(hits) == 0L)) {
    cli::cli_abort("Some points in `data` did not overlay any areal grid cell.")
  }
  multi <- which(lengths(hits) > 1L)
  if (length(multi)) {
    cli::cli_warn(c(
      "Some points overlay multiple areal grid cells; using the first match.",
      "i" = "This can happen for points exactly on polygon boundaries."
    ))
  }

  out_data <- data
  out_data[[space_column]] <- grid[[space_column]][vapply(hits, `[[`, integer(1), 1L)]
  domain <- make_areal_domain(
    data = out_data,
    spatial_domain = grid,
    space_column = space_column,
    id_column = space_column,
    adjacency = adjacency
  )

  list(
    data = out_data,
    grid = grid,
    domain = domain
  )
}

.validate_areal_matrix_structure <- function(W, from_igraph = FALSE) {
  if (!inherits(W, "sparseMatrix")) {
    cli::cli_abort("Internal areal matrix must be sparse.")
  }
  if (nrow(W) != ncol(W)) {
    cli::cli_abort("`spatial_domain` matrix must be square.")
  }
  rn <- rownames(W)
  cn <- colnames(W)
  if (is.null(rn) || is.null(cn)) {
    cli::cli_abort("`spatial_domain` matrix must have row and column names.")
  }
  if (!identical(rn, cn)) {
    cli::cli_abort("`spatial_domain` row and column names must be identical and in the same order.")
  }
  if (anyNA(rn) || any(!nzchar(rn))) {
    cli::cli_abort("`spatial_domain` unit names cannot contain missing or empty values.")
  }
  if (anyDuplicated(rn)) {
    cli::cli_abort("`spatial_domain` unit names must be unique.")
  }
  if (!from_igraph && any(Matrix::diag(W) != 0)) {
    cli::cli_abort("`spatial_domain` matrix must have a zero diagonal (no self-neighbors).")
  }
}

.validate_areal_matrix_values <- function(W, stage) {
  if (any(!is.finite(W@x))) {
    cli::cli_abort("`spatial_domain` matrix has non-finite values {stage}.")
  }
  if (any(W@x < 0)) {
    cli::cli_abort("`spatial_domain` matrix has negative weights {stage}.")
  }
  if (any(Matrix::diag(W) != 0)) {
    cli::cli_abort("`spatial_domain` matrix must have a zero diagonal {stage}.")
  }
}

.normalize_areal_W <- function(W) {
  rs <- Matrix::rowSums(W)
  islands <- which(rs == 0)
  if (length(islands) > 0L) {
    island_names <- rownames(W)[islands]
    cli::cli_warn(c(
      "Found areal units with no neighbors (islands).",
      "i" = paste(island_names, collapse = ", "),
      "i" = "These units will be treated as conditionally independent spatial effects."
    ))
  }

  non_islands <- which(rs > 0)
  if (!length(non_islands)) return(W)

  W[non_islands, ] <- Matrix::Diagonal(x = 1 / rs[non_islands]) %*%
    W[non_islands, , drop = FALSE]
  methods::as(W, "dgCMatrix")
}

.validate_areal_row_sums <- function(W, tol = 1e-8) {
  rs <- Matrix::rowSums(W)
  non_islands <- which(rs > 0)
  if (!length(non_islands)) return(invisible(NULL))

  bad <- non_islands[abs(rs[non_islands] - 1) > tol]
  if (length(bad)) {
    cli::cli_abort(c(
      "Normalized areal matrix rows do not sum to 1 for non-island units.",
      "i" = paste(rownames(W)[bad], collapse = ", ")
    ))
  }
  invisible(NULL)
}

.validate_areal_data_membership <- function(data, unit_names, space_column) {
  if (!space_column %in% names(data)) {
    cli::cli_abort("Column {.field {space_column}} was not found in `data`.")
  }
  values <- data[[space_column]]
  if (anyNA(values)) {
    cli::cli_abort("`data[[{space_column}]]` contains missing areal memberships.")
  }
  values <- as.character(values)
  missing_units <- setdiff(unique(values), unit_names)
  if (length(missing_units)) {
    cli::cli_abort(c(
      "Some `data[[{space_column}]]` values are not present in the domain.",
      "i" = paste(missing_units, collapse = ", ")
    ))
  }
  invisible(NULL)
}

.as_general_dgC <- function(x) {
  if (!inherits(x, "sparseMatrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  if (inherits(x, "symmetricMatrix")) {
    x <- methods::as(x, "generalMatrix")
  }
  methods::as(x, "dgCMatrix")
}

.as_sf_polygons <- function(x, arg = "spatial_domain") {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg sf} must be installed to use `sf` areal domains.")
  }
  if (inherits(x, "sfc")) {
    x <- sf::st_sf(geometry = x)
  }
  if (!inherits(x, "sf")) {
    cli::cli_abort("`{arg}` must be an `sf` or `sfc` polygon object.")
  }
  geom_type <- as.character(sf::st_geometry_type(x))
  if (!all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    cli::cli_abort("`{arg}` must contain only polygon geometries.")
  }
  x
}

.align_sf_crs <- function(x, target) {
  if (is.na(sf::st_crs(target)) || is.na(sf::st_crs(x))) {
    return(x)
  }
  if (!identical(sf::st_crs(x), sf::st_crs(target))) {
    x <- sf::st_transform(x, sf::st_crs(target))
  }
  x
}

.sf_unit_names <- function(x, id_column = NULL) {
  if (is.null(id_column)) {
    return(sprintf("area_%03d", seq_len(nrow(x))))
  }
  if (!is.character(id_column) || length(id_column) != 1L || is.na(id_column) || !nzchar(id_column)) {
    cli::cli_abort("`id_column` must be a single non-empty column name.")
  }
  if (!id_column %in% names(x)) {
    cli::cli_abort("`id_column` {.field {id_column}} was not found in `spatial_domain`.")
  }
  unit_names <- as.character(x[[id_column]])
  if (anyNA(unit_names) || any(!nzchar(unit_names))) {
    cli::cli_abort("`id_column` cannot contain missing or empty values.")
  }
  if (anyDuplicated(unit_names)) {
    cli::cli_abort("`id_column` values must be unique.")
  }
  unit_names
}

.sf_areal_adjacency_matrix <- function(spatial_domain, id_column = NULL,
                                       adjacency = c("rook", "queen")) {
  adjacency <- match.arg(adjacency)
  sf_domain <- .as_sf_polygons(spatial_domain)
  unit_names <- .sf_unit_names(sf_domain, id_column = id_column)

  neighbors <- switch(adjacency,
    rook = sf::st_relate(sf_domain, sf_domain, pattern = "F***1****"),
    queen = sf::st_touches(sf_domain, sf_domain)
  )
  from <- rep(seq_along(neighbors), lengths(neighbors))
  to <- unlist(neighbors, use.names = FALSE)
  keep <- from != to
  from <- from[keep]
  to <- to[keep]
  Matrix::sparseMatrix(
    i = from,
    j = to,
    x = 1,
    dims = c(length(unit_names), length(unit_names)),
    dimnames = list(unit_names, unit_names)
  )
}

is_areal_domain <- function(x) inherits(x, "sdmTMBareal")

is_areal_fit <- function(object) {
  if (!is.null(object$tmb_data) && "spatial_model" %in% names(object$tmb_data)) {
    if (object$tmb_data$spatial_model %in% c(1L, 2L)) {
      return(TRUE)
    }
  }
  !is.null(object$spde) && is_areal_domain(object$spde)
}

is_car_fit <- function(object) {
  !is.null(object$tmb_data) &&
    "spatial_model" %in% names(object$tmb_data) &&
    identical(object$tmb_data$spatial_model, 2L)
}

domain_n_s <- function(mesh) {
  if (is_areal_domain(mesh)) mesh$n_s else nrow(mesh$mesh$loc)
}

domain_obs_index <- function(mesh, data) {
  match(as.character(data[[mesh$space_column]]), mesh$unit_names)
}

domain_pred_index <- function(mesh, newdata) {
  match(as.character(newdata[[mesh$space_column]]), mesh$unit_names)
}

areal_projection_matrix <- function(mesh, data) {
  idx <- domain_obs_index(mesh, data)
  if (anyNA(idx)) {
    missing_units <- unique(as.character(data[[mesh$space_column]])[is.na(idx)])
    cli::cli_abort(c(
      "Some areal units were not found in the domain.",
      "i" = paste(missing_units, collapse = ", ")
    ))
  }
  Matrix::sparseMatrix(
    i = seq_len(nrow(data)),
    j = idx,
    x = 1,
    dims = c(nrow(data), mesh$n_s)
  )
}

prepare_spatial_domain <- function(mesh, data, mesh_missing, share_range = TRUE,
                                   anisotropy = FALSE, covariate_diffusion, priors = sdmTMBpriors(),
                                   normalize, experimental = NULL,
                                   share_range_user = share_range,
                                   spatial_model = c("spde", "sar", "car")) {
  spatial_model <- match.arg(tolower(spatial_model[1L]), c("spde", "sar", "car"))
  is_areal <- !mesh_missing && is_areal_domain(mesh)
  spde <- mesh

  if (spatial_model == "spde" && is_areal) {
    cli::cli_abort("`spatial_model = \"spde\"` requires an SPDE mesh from `make_mesh()`. Use `spatial_model = \"sar\"` or `\"car\"` with an areal domain.")
  }
  if (spatial_model %in% c("sar", "car") && !is_areal) {
    cli::cli_abort("`spatial_model = \"{spatial_model}\"` requires an areal domain from `make_areal_domain()` or `make_areal_grid()`.")
  }

  if (is_areal) {
    space_column <- spde$space_column
    if (!space_column %in% names(data)) {
      cli::cli_abort("Column {.field {space_column}} from the areal domain was not found in `data`.")
    }
    if (anyNA(data[[space_column]])) {
      cli::cli_abort("Areal domain column {.field {space_column}} contains missing values.")
    }
    areal_idx <- domain_obs_index(spde, data)
    if (anyNA(areal_idx)) {
      cli::cli_abort("Some values in `data[[{space_column}]]` have no match in the areal domain.")
    }

    if (isTRUE(anisotropy)) {
      cli::cli_abort("`anisotropy` is not supported with areal domains.")
    }
    if ("spde_barrier" %in% names(spde)) {
      cli::cli_abort("Barrier models are not supported with areal domains.")
    }
    if (!is.null(covariate_diffusion)) {
      cli::cli_abort("`covariate_diffusion` is not supported with areal domains.")
    }
    share_range_user <- unlist(share_range_user)
    if (any(!share_range_user)) {
      cli::cli_abort("`share_range = FALSE` is not supported with areal domains in v1.")
    }
    if (any(c(!is.na(priors$matern_s[1:2]), !is.na(priors$matern_st[1:2])))) {
      cli::cli_abort("PC Matern priors are not supported with areal domains.")
    }
    if (!is.null(experimental) &&
        "epsilon_model" %in% names(experimental) &&
        !is.null(experimental$epsilon_model)) {
      cli::cli_abort("`experimental$epsilon_model` is not supported with areal domains.")
    }
    if (spatial_model == "car" && !Matrix::isSymmetric(spde$W_raw)) {
      cli::cli_abort("`spatial_model = \"car\"` requires a symmetric areal adjacency matrix.")
    }

    return(list(
      type = "areal",
      spatial_model = spatial_model,
      spde = spde,
      n_s = spde$n_s,
      A_st = areal_projection_matrix(spde, data),
      A_spatial_index = seq_len(nrow(data)) - 1L,
      W_ss = if (spatial_model == "car") spde$W_raw else spde$W,
      normalize_in_r = 0L,
      barrier = 0L,
      barrier_scaling = c(1, 1),
      anisotropy = 0L,
      spde_struct = dummy_spde(),
      spde_aniso_struct = dummy_spde_aniso(),
      spde_barrier_struct = dummy_spde_barrier()
    ))
  }

  barrier <- "spde_barrier" %in% names(spde)
  if (barrier && anisotropy) {
    cli::cli_warn("Using a barrier mesh; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }
  if (any(c(!is.na(priors$matern_s[1:2]), !is.na(priors$matern_st[1:2]))) && anisotropy) {
    cli::cli_warn("Using PC Matern priors; therefore, anistropy will be disabled.")
    anisotropy <- FALSE
  }
  if (!"A_st" %in% names(spde)) {
    cli::cli_abort("`mesh` was created with an old version of `make_mesh()`.")
  }

  list(
    type = "spde",
    spatial_model = spatial_model,
    spde = spde,
    n_s = nrow(spde$mesh$loc),
    A_st = spde$A_st,
    A_spatial_index = spde$sdm_spatial_id - 1L,
    W_ss = dummy_sparse_1x1(),
    normalize_in_r = as.integer(normalize),
    barrier = as.integer(barrier),
    barrier_scaling = if (barrier) spde$barrier_scaling else c(1, 1),
    anisotropy = as.integer(anisotropy),
    spde_struct = get_spde_matrices(spde),
    spde_aniso_struct = make_anisotropy_spde(spde, anisotropy),
    spde_barrier_struct = make_barrier_spde(spde)
  )
}

dummy_sparse_1x1 <- function() {
  Matrix::sparseMatrix(i = 1L, j = 1L, x = 0, dims = c(1L, 1L))
}

dummy_spde <- function() {
  m <- dummy_sparse_1x1()
  list(M0 = m, M1 = m, M2 = m)
}

dummy_spde_aniso <- function() {
  m <- Matrix::Matrix(0, 1, 1, sparse = TRUE)
  list(
    n_s = 0L,
    n_tri = 0L,
    Tri_Area = rep(0, 1),
    E0 = matrix(0, 1),
    E1 = matrix(0, 1),
    E2 = matrix(0, 1),
    TV = matrix(0, 1),
    G0 = m,
    G0_inv = m
  )
}

dummy_spde_barrier <- function() {
  m <- dummy_sparse_1x1()
  list(C0 = rep(1, 2), C1 = rep(1, 2), D0 = m, D1 = m, I = m)
}
