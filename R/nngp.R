# Mesh-based nearest-neighbour Gaussian process (NNGP) prior for random
# fields (`spatial_model = "nngp"`, RTMB backend only).
#
# The field is still defined at mesh vertices and interpolated with the usual
# A matrices; only its prior changes. Vertices are visited in a max-min
# order, and each vertex is conditioned on up to `m` nearby preceding
# vertices (Datta et al. 2016; Guinness 2018), which defines a proper
# Gaussian distribution:
#
#   w_i | w[N_i] ~ Normal(sum(b_i * w[N_i]), sigma * sqrt(f_i)),
#   b_i = solve(R[N_i, N_i], R[N_i, i]),  f_i = 1 - sum(R[N_i, i] * b_i),
#
# where R is the parent correlation. The ordering and neighbour search follow
# tinyVAST's NNGP (Jim Thorson), using GPvecchia and GpGp.
#
# `make_nngp_data()` builds the data-only graph once per fit;
# `nngp_factor()` turns it into weights and conditional SDs for one range;
# `nngp_nll()` and `simulate_nngp()` evaluate and draw fields with them.

# Fixed parent nugget keeping local correlation matrices well conditioned:
# R_delta = (R + delta I) / (1 + delta), so the diagonal stays at one.
nngp_nugget <- 1e-10

# Parent correlation at inverse range `kappa`. The exponential (Matérn
# nu = 1/2) is scaled so that the correlation at sdmTMB's conventional range
# sqrt(8) / kappa is exp(-2), close to the Matérn nu = 1 correlation there.
# The Matérn nu = 1 is the covariance the SPDE approximates; zero distances
# (data) are set to one directly since K_1(0) is infinite.
nngp_correlation <- function(distance, kappa, covariance = "exponential") {
  if (covariance == "exponential") return(exp(-kappa * distance / sqrt(2)))
  "[<-" <- RTMB::ADoverload("[<-")
  out <- distance * 0 + 1
  nonzero <- distance > 0
  x <- kappa * distance[nonzero]
  if (length(x)) out[nonzero] <- x * RTMB::besselK(x, 1)
  out
}

# Ordering and predecessor graph for vertex coordinates `coords`. Lists and
# rows are indexed by ordered position; `order`, `neighbours`, and `index` hold
# one-based original vertex indices, and distances use the original
# coordinates. Row i of `index` holds the neighbours of ordered vertex i,
# padded with vertex 1 (given weight zero by `nngp_factor()`). `covariance`
# is the parent correlation passed to `nngp_correlation()`.
make_nngp_data <- function(coords, nngp_neighbours,
                           covariance = c("exponential", "matern")) {
  covariance <- match.arg(covariance)
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || ncol(coords) != 2L || !all(is.finite(coords))) {
    cli_abort("NNGP vertex coordinates must be a finite two-column matrix.")
  }
  if (anyDuplicated(coords)) {
    cli_abort("NNGP mesh vertices must have unique coordinates.")
  }
  n <- nrow(coords)
  if (n == 1L) {
    order <- 1L
    neighbours <- list(integer(0))
  } else {
    rlang::check_installed(c("GPvecchia", "GpGp", "withr"),
      reason = "for `spatial_model = \"nngp\"`.")
    # GpGp perturbs coordinates at random in its neighbour search; a fixed
    # seed makes the graph reproducible without changing the caller's RNG.
    withr::with_seed(1L, {
      order <- GPvecchia::order_maxmin_exact(coords)
      nn <- GpGp::find_ordered_nn(coords[order, , drop = FALSE],
        nngp_neighbours)
    })
    # Column 1 is the vertex itself; missing predecessors are NA.
    neighbours <- lapply(seq_len(n), function(i) {
      j <- nn[i, -1L]
      order[j[!is.na(j)]]
    })
  }
  distance_to <- lapply(seq_len(n), function(i) {
    sqrt(colSums((t(coords[neighbours[[i]], , drop = FALSE]) -
      coords[order[[i]], ])^2))
  })
  distance_between <- lapply(neighbours, function(j) {
    unname(as.matrix(stats::dist(coords[j, , drop = FALSE])))
  })
  k <- lengths(neighbours)
  index <- matrix(1L, n, max(1L, k))
  for (i in which(k > 0L)) index[i, seq_len(k[[i]])] <- neighbours[[i]]
  list(order = as.integer(order), neighbours = neighbours, index = index,
    distance_to = distance_to, distance_between = distance_between,
    m = as.integer(nngp_neighbours), covariance = covariance)
}

# Conditional weights and SDs (for unit parent variance) at inverse range
# `kappa`, with `weights` aligned to `graph$index`. The RTMB field code uses
# this in place of a precision matrix.
nngp_factor <- function(graph, kappa) {
  "[<-" <- RTMB::ADoverload("[<-")
  k <- lengths(graph$neighbours)
  weights <- matrix(0, nrow(graph$index), ncol(graph$index))
  variance <- rep(1, length(k))
  regularize <- 1 / (1 + nngp_nugget)
  covariance <- graph$covariance
  for (i in which(k > 0L)) {
    R <- (nngp_correlation(graph$distance_between[[i]], kappa, covariance) +
      diag(nngp_nugget, k[[i]])) * regularize
    r <- nngp_correlation(graph$distance_to[[i]], kappa, covariance) * regularize
    b <- RTMB::solve(R, r)
    weights[i, seq_len(k[[i]])] <- b
    variance[i] <- 1 - sum(r * b)
  }
  structure(list(order = graph$order, index = graph$index, weights = weights,
    sd = sqrt(variance)), class = "nngp_factor")
}

# Negative log density of fields `x` (a vector, or one field per column)
# with parent SD `sigma`.
nngp_nll <- function(x, factor, sigma) {
  n <- length(factor$order)
  dim(x) <- c(n, length(x) / n)
  nll <- 0
  for (j in seq_len(ncol(x))) {
    field <- x[, j]
    neighbours <- field[factor$index]
    dim(neighbours) <- dim(factor$index)
    mean <- (factor$weights * neighbours) %*% rep(1, ncol(neighbours))
    nll <- nll - sum(RTMB::dnorm(field[factor$order], mean,
      sigma * factor$sd, log = TRUE))
  }
  nll
}

# `n` independent fields (columns) with parent SD `sigma`, drawn forward in
# the NNGP order. Plain numbers only.
simulate_nngp <- function(factor, sigma, n = 1L) {
  x <- matrix(0, length(factor$order), n)
  for (i in seq_along(factor$order)) {
    mean <- colSums(factor$weights[i, ] * x[factor$index[i, ], , drop = FALSE])
    x[factor$order[[i]], ] <- mean + sigma * factor$sd[[i]] * stats::rnorm(n)
  }
  x
}
