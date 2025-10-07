mesh <- fmesher::fm_mesh_2d(st_coordinates(sf_loc)[, 1:2], cutoff = 0.1)

# Compute finite element matrices for the SPDE approach
spde <- fm_fem(mesh)
invM0 <- invsqrtM0 <- spde$c0
diag(invsqrtM0) <- 1 / sqrt(diag(spde$c0))
diag(invM0) <- 1 / diag(spde$c0)
