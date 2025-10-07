
# SIMULATION OF COVARIATE-DEPENDENT SPATIAL CORRELATION
# This script demonstrates covariate-dependent anisotropy in spatial Gaussian random fields
# using the SPDE (Stochastic Partial Differential Equation) approach. The key innovation is
# allowing the spatial scale parameter (kappa) to vary as a function of covariates, creating
# heterogeneous correlation patterns across space. The simulation generates data with a 
# ring-shaped habitat preference pattern and fits a spatial field where correlation changes
# with distance from the origin - areas closer to the preferred ring distance have different
# spatial correlation structure than areas farther away.
# 
# The statistical framework uses a Matérn covariance function implemented via SPDE on a 
# triangular mesh, where the precision matrix Q incorporates both isotropic (M0) and 
# anisotropic (M1, M2) stiffness matrices. Stiffness matrices are fundamental components
# of the SPDE approach that encode spatial relationships between mesh vertices - they
# represent how "stiff" or resistant the spatial field is to changes between neighboring
# locations. The isotropic stiffness matrix M0 captures uniform spatial correlation in
# all directions, while anisotropic matrices M1 and M2 allow correlation to vary by
# direction and be influenced by covariates. The covariate-dependent kappa allows the 
# spatial range to shrink or expand based on environmental gradients, which is particularly 
# useful for modeling species distributions where spatial autocorrelation varies with habitat
# suitability. Model comparison evaluates whether this covariate-dependent approach provides
# better fit than traditional constant correlation models.

# JT: Install development version of tinyVAST
# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force = TRUE )

# JT: Set working directory for this covariate anisotropy SPDE simulation
setwd( R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-07 -- covariate anisotropy SPDE)' )

# Load required libraries
library(tinyVAST)  # For spatial modeling with covariate-dependent anisotropy
library(fmesher)   # For mesh construction and SPDE matrices
library(ggplot2)   # For plotting
library(Matrix)    # For sparse matrix operations
library(sf)        # For spatial data handling

# JT: Euclidean and polar coordinates
# Generate spatial coordinates using Poisson disk sampling for even distribution
euclidean_iz = pracma::poisson2disk( n = 1000 ) - 0.5  # Center around (0,0)
# Convert Cartesian to polar coordinates
polar_iz = cbind( d = sqrt(rowSums(euclidean_iz^2)),              # Distance from origin
                  theta = atan2(euclidean_iz[,2],euclidean_iz[,1]) )  # Angle

# JT: Simulate preference
# Create habitat preference that peaks at distance 0.3 from origin (ring pattern)
data_i = data.frame( "X" = euclidean_iz[,1], "Y" = euclidean_iz[,2],
                     polar_iz, "p" = -1 * abs(polar_iz[,'d'] - 0.3) / 0.2 )

# JT: Make plotting grid
# Create square domain centered at origin
domain = st_polygon( list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0,0)) - 0.5) )
# Create fine-resolution grid for visualization
sf_grid = st_make_grid( domain, cellsize = 0.01 )

# JT: Show preference
# Plot distance from origin
ggplot( data_i ) + geom_point( aes(x=X, y=Y, col=d) )
# Plot preference values (ring pattern)
ggplot( data_i ) + geom_point( aes(x=X, y=Y, col=p) )

# JT: Make mesh
# Create triangular mesh using observation locations
mesh = fm_mesh_2d(
  loc = euclidean_iz,  # Use sampling locations as mesh vertices
  cutoff = 0           # No minimum distance between vertices
)
# Create finite element matrices for SPDE
spde = fm_fem( mesh )
# Create projection matrices from mesh to grid centroids
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(sf_grid)) )$proj$A
# Create projection matrices from mesh to observation locations
A_is = fm_evaluator( mesh, loc=as.matrix(data_i[,c("X","Y")]) )$proj$A

# JT: Visualize mesh
# Plot mesh structure
plot(mesh)
# Overlay observation points
points( euclidean_iz )

# JT: Access internal function and add covariates to mesh
#make_stiffness = tinyVAST:::make_stiffness
#source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\internal.R)')
# Access internal stiffness matrix function
make_stiffness = tinyVAST:::make_stiffness
# Add distance and angle covariates to mesh vertices
mesh = add_vertex_covariates(
  mesh = mesh,
  data = data_i,
  covariates = c("d","theta"),  # Distance and angle from origin
  coords  = c("X", "Y")
)
# Create anisotropic stiffness matrix with covariate dependence
M1 = make_stiffness(
  mesh,
  loc = cbind( mesh$loc[,1:2], as.matrix(mesh$vertex) ),
  H = diag( 10 ^ (c(1, 1, -1, 1)) )  # Anisotropy parameters
)
# Set spatial scale parameter
kappa = 100
# Get isotropic stiffness matrix from SPDE
M0 = spde$c0

# JT: Make indicator to show diffusion
# Find mesh vertex nearest to point (0.3, 0)
nn = RANN::nn2( data = mesh$loc[,1:2], query = matrix(c(0.3,0),nrow=1), k = 1 )
# Create impulse function at that vertex
v0_s = v = rep(0, mesh$n)
v0_s[nn$nn.idx[1]] = 1

# JT: Plot diffusion
# Solve diffusion equation with anisotropic stiffness
invD = Diagonal( n = mesh$n ) + 1/kappa^2 * solve(M0) %*% M1
v1_s = solve(invD, v0_s)
# Project to grid and visualize diffusion pattern
sf_plot = st_sf( sf_grid, v = as.numeric(A_gs %*% v1_s) )
plot(sf_plot, border = NA)

# Simulate spatial random field with covariate-dependent anisotropy
tau = 0.0001  # Precision parameter
# Create second-order stiffness matrix
M2 = M1 %*% solve(M0) %*% M1
# Construct precision matrix for Matérn field with ν=2
Q = tau^2 * ( kappa^4 * M0 + 2*kappa^2*M1 + M2 )
# Sample from multivariate normal with precision Q
x_s = rmvnorm_prec( prec = Q )[,1]
# Project spatial field to observation locations
data_i$x_i = as.numeric(A_is %*% x_s)
# Generate noisy observations
data_i$y_i = rnorm( n = nrow(data_i), mean = data_i$x_i, sd = 0.1 )

# JT: Plot field
# Visualize the simulated spatial field
sf_plot = st_sf( sf_grid, v = as.numeric(A_gs %*% x_s) )
plot(sf_plot, border = NA)

# Plot observed data with spatial field values
ggplot( data_i ) + geom_point( aes( x=X, y=Y, col = y_i) )


# JT: Fit with covariate-anisotropy
# Fit model with distance-dependent spatial correlation
fit = tinyVAST(
  data = data_i,
  formula = y_i ~ 1 + poly(d, 2),  # Quadratic distance effect
  spatial_domain = mesh,
  space_term = "",
  space_columns = c("X","Y"),
  development = list(
    kappa_formula = ~ d  # Distance-dependent kappa parameter
  ),
  control = tinyVASTcontrol( trace = 1, use_anisotropy = TRUE)
)

# JT: Plot diffusion
# Visualize fitted diffusion pattern with estimated parameters
invD = Diagonal( n = mesh$n ) + exp(-2 * fit$opt$par['log_kappa']) * solve(M0) %*% fit$rep$G1
v1_s = solve(invD, v0_s)
sf_plot = st_sf( sf_grid, v = as.numeric(A_gs %*% v1_s) )
plot(sf_plot, border = NA)

# JT: Fit without covariate-anisotropy
# Fit comparison model with constant spatial correlation
fit0 = tinyVAST(
  data = data_i,
  formula = y_i ~ 1 + poly(d, 2),  # Same fixed effects
  spatial_domain = mesh,
  space_term = "",
  space_columns = c("X","Y"),
  development = list(
    #kappa_formula = ~ d  # No distance-dependent kappa
  ),
  control = tinyVASTcontrol( trace = 1, use_anisotropy = TRUE)
)

# Compare model performance
performance = rbind(
  AIC = c(AIC(fit), AIC(fit0)),      # Akaike Information Criterion
  cAIC = c(cAIC(fit), cAIC(fit0)),   # Conditional AIC
  CV = c(cv::cv(fit)[['CV crit']], cv::cv(fit0)[['CV crit']])  # Cross-validation
)
performance
