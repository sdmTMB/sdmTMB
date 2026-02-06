suppressPackageStartupMessages({
  library(pkgload)
  library(Matrix)
  library(ggplot2)
})

# Use local package code from this repo
pkgload::load_all(".", quiet = TRUE)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
# Use one year to keep visuals simple and avoid time-slice indexing.
dat <- as.data.frame(sdmTMB::pcod_2011)
dat <- dat[is.finite(dat$depth) & !is.na(dat$year), , drop = FALSE]
dat <- dat[dat$year == min(dat$year), , drop = FALSE]
dat$depth_scaled <- as.numeric(scale(dat$depth))

mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 5)
A <- mesh$A_st

x_obs <- dat$depth_scaled
obs_x <- dat$X
obs_y <- dat$Y
vertex_x <- mesh$mesh$loc[, 1]
vertex_y <- mesh$mesh$loc[, 2]

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# Return nearest observation index for each query location.
nearest_obs_index <- function(query_x, query_y, obs_x, obs_y) {
  d2 <- (outer(query_x, obs_x, "-"))^2 + (outer(query_y, obs_y, "-"))^2
  max.col(-d2, ties.method = "first")
}

# Evaluate bilinear interpolation on a regular grid at query points.
# z[iy, ix] corresponds to (xg[ix], yg[iy]).
bilinear_eval <- function(query_x, query_y, xg, yg, z) {
  ix <- findInterval(query_x, xg, rightmost.closed = TRUE)
  iy <- findInterval(query_y, yg, rightmost.closed = TRUE)

  # Clamp to valid cell indices [1, n-1]
  ix <- pmax(1L, pmin(ix, length(xg) - 1L))
  iy <- pmax(1L, pmin(iy, length(yg) - 1L))

  x1 <- xg[ix]
  x2 <- xg[ix + 1L]
  y1 <- yg[iy]
  y2 <- yg[iy + 1L]

  q11 <- z[cbind(iy, ix)]
  q21 <- z[cbind(iy, ix + 1L)]
  q12 <- z[cbind(iy + 1L, ix)]
  q22 <- z[cbind(iy + 1L, ix + 1L)]

  wx <- (query_x - x1) / (x2 - x1)
  wy <- (query_y - y1) / (y2 - y1)

  q11 * (1 - wx) * (1 - wy) +
    q21 * wx * (1 - wy) +
    q12 * (1 - wx) * wy +
    q22 * wx * wy
}

# IDW interpolation using sf + gstat (point -> vertex).
idw_vertex_values <- function(obs_x, obs_y, obs_val, vertex_x, vertex_y, power = 2) {
  obs_df <- data.frame(X = obs_x, Y = obs_y, value = obs_val)
  pred_df <- data.frame(X = vertex_x, Y = vertex_y)

  data_sf <- sf::st_as_sf(obs_df, coords = c("X", "Y"))
  pred_sf <- sf::st_as_sf(pred_df, coords = c("X", "Y"))

  idw_model <- gstat::gstat(
    formula = value ~ 1,
    data = data_sf,
    set = list(idp = power)
  )

  invisible(capture.output({
    predicted <- stats::predict(idw_model, pred_sf)
  }))
  as.numeric(predicted$var1.pred)
}

# -----------------------------------------------------------------------------
# Method 1: normalized A^T back-projection (current approach)
# -----------------------------------------------------------------------------
# This is not a true matrix inverse. It's:
#   vertex = (A^T x) / (A^T 1)
# where division is element-wise on supported vertices.
n_i <- length(x_obs)
support <- as.vector(Matrix::t(A) %*% rep(1, n_i))
numerator <- as.vector(Matrix::t(A) %*% x_obs)

x_vertex_at <- rep(0, length(support))
keep <- support > 0
x_vertex_at[keep] <- numerator[keep] / support[keep]

# -----------------------------------------------------------------------------
# Method 2: nearest-neighbour point -> vertex assignment
# -----------------------------------------------------------------------------
idx_nn_vertex <- nearest_obs_index(vertex_x, vertex_y, obs_x, obs_y)
x_vertex_nn <- x_obs[idx_nn_vertex]

# -----------------------------------------------------------------------------
# Method 3: homemade bilinear interpolation
# -----------------------------------------------------------------------------
# Build a regular grid over observation coordinates.
# Fill grid-node values by nearest observed value, then evaluate bilinear at vertices.
nx <- 50L
ny <- 50L
xg <- seq(min(obs_x), max(obs_x), length.out = nx)
yg <- seq(min(obs_y), max(obs_y), length.out = ny)
grid_xy <- expand.grid(X = xg, Y = yg)

idx_nn_grid <- nearest_obs_index(grid_xy$X, grid_xy$Y, obs_x, obs_y)
z_grid <- matrix(
  x_obs[idx_nn_grid],
  nrow = ny,
  ncol = nx,
  byrow = TRUE
)

x_vertex_bilinear <- bilinear_eval(vertex_x, vertex_y, xg, yg, z_grid)

# -----------------------------------------------------------------------------
# Project all methods back to observation points with A %*% vertex
# -----------------------------------------------------------------------------
methods <- list(
  "A^T normalized" = x_vertex_at,
  "Nearest neighbor" = x_vertex_nn,
  "Bilinear (regular grid)" = x_vertex_bilinear
)
if (requireNamespace("sf", quietly = TRUE) && requireNamespace("gstat", quietly = TRUE)) {
  methods[["IDW (gstat)"]] <- idw_vertex_values(
    obs_x = obs_x,
    obs_y = obs_y,
    obs_val = x_obs,
    vertex_x = vertex_x,
    vertex_y = vertex_y,
    power = 2
  )
} else {
  cat("Skipping IDW method: install packages 'sf' and 'gstat' to enable it.\n")
}

x_hat <- lapply(methods, function(v) as.vector(A %*% v))
resids <- lapply(x_hat, function(v) x_obs - v)

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
row_sums_A <- Matrix::rowSums(A)

metrics <- do.call(rbind, lapply(names(x_hat), function(m) {
  pred <- x_hat[[m]]
  data.frame(
    method = m,
    correlation = cor(x_obs, pred),
    rmse = sqrt(mean((x_obs - pred)^2)),
    mae = mean(abs(x_obs - pred)),
    max_abs_error = max(abs(x_obs - pred)),
    stringsAsFactors = FALSE
  )
}))

cat("Projection diagnostics (lower RMSE/MAE is better)\n")
print(metrics, row.names = FALSE)
cat("\nA row-sum check (should be ~1 for barycentric rows):\n")
cat(sprintf("row_sum_min=%.6f row_sum_max=%.6f\n", min(row_sums_A), max(row_sums_A)))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
out_dir <- "scratch"

# 1) Scatter comparison by method
scatter_df <- do.call(rbind, lapply(names(x_hat), function(m) {
  data.frame(method = m, observed = x_obs, projected = x_hat[[m]])
}))

p_scatter <- ggplot(scatter_df, aes(observed, projected)) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  coord_equal() +
  facet_wrap(~ method, ncol = 3) +
  labs(
    title = "Point-to-Vertex Strategies: Observed vs Projected Back",
    subtitle = "All methods project vertex values back to points with A %*% vertex",
    x = "Observed depth_scaled",
    y = "Projected-back value"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "a-matrix-methods-scatter.png"),
  plot = p_scatter,
  width = 13,
  height = 4.8,
  dpi = 140
)

# 2) Map panels: observed + each projected field
map_pred_df <- rbind(
  data.frame(X = dat$X, Y = dat$Y, value = x_obs, panel = "Observed"),
  do.call(rbind, lapply(names(x_hat), function(m) {
    data.frame(X = dat$X, Y = dat$Y, value = x_hat[[m]], panel = paste0("Projected: ", m))
  }))
)

p_pred_map <- ggplot(map_pred_df, aes(X, Y, color = value)) +
  geom_point(size = 2) +
  coord_equal() +
  facet_wrap(~ panel, ncol = 2) +
  scale_color_viridis_c(option = "C") +
  labs(
    title = "Observed vs Projected-Back Fields",
    color = "value"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "a-matrix-methods-map-predicted.png"),
  plot = p_pred_map,
  width = 10.5,
  height = 8,
  dpi = 140
)

# 3) Residual map panels by method
map_resid_df <- do.call(rbind, lapply(names(resids), function(m) {
  data.frame(X = dat$X, Y = dat$Y, value = resids[[m]], panel = paste0("Residual: ", m))
}))

p_resid_map <- ggplot(map_resid_df, aes(X, Y, color = value)) +
  geom_point(size = 2) +
  coord_equal() +
  facet_wrap(~ panel, ncol = 3) +
  scale_color_gradient2(midpoint = 0) +
  labs(
    title = "Residuals (Observed - Projected Back)",
    color = "residual"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "a-matrix-methods-map-residuals.png"),
  plot = p_resid_map,
  width = 13,
  height = 4.8,
  dpi = 140
)

# 4) Vertex support map for the A^T-normalized method
vertex_df <- data.frame(
  X = vertex_x,
  Y = vertex_y,
  support = support,
  x_vertex = x_vertex_at
)

p_support <- ggplot(vertex_df, aes(X, Y, size = support, color = x_vertex)) +
  geom_point(alpha = 0.85) +
  coord_equal() +
  scale_color_viridis_c(option = "D") +
  labs(
    title = "A^T-Normalized Vertex Values and Support",
    subtitle = "Size = A^T 1 support; Color = projected vertex value",
    color = "vertex value",
    size = "support"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "a-matrix-methods-vertex-support.png"),
  plot = p_support,
  width = 7.2,
  height = 6,
  dpi = 140
)

cat("\nWrote plots to:\n")
cat("- scratch/a-matrix-methods-scatter.png\n")
cat("- scratch/a-matrix-methods-map-predicted.png\n")
cat("- scratch/a-matrix-methods-map-residuals.png\n")
cat("- scratch/a-matrix-methods-vertex-support.png\n")
