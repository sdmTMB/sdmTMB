#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_STRING(method); // which method to use (e.g., "diffusion")
  DATA_STRING(dist); // distribution (this example uses "LNP" (lognormal-Poisson))
  DATA_VECTOR(c_i);  // counts for observation i
  DATA_VECTOR(weights_i); // weights for calculating cAIC
  DATA_VECTOR(depth_s);  // counts for observation i

  // SPDE objects (computed from mesh)
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(invsqrtM0);
  DATA_SPARSE_MATRIX(invM0);

  // Projection matrices from mesh nodes to:
  DATA_SPARSE_MATRIX(A_is); // observation level, i
  DATA_SPARSE_MATRIX(A_gs); // prediction grid, g

  // Parameters
  PARAMETER(beta0); // Intercept
  PARAMETER_VECTOR(beta_j); // depth coefficients (linear + quadratic)
  PARAMETER(ln_tau); // scalar of the precision matrix
  PARAMETER(ln_kappa); // decorrelation rate

  // Random effects
  PARAMETER_VECTOR(omega_s); // spatial random effects

  // Objective function
  Type jnll = 0.0; // joint negative log-likelihood

  // Derived quantities
  // Q = precision matrix for the Gaussian random field
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 +
    Type(2.0)*exp(2*ln_kappa)* M1 + M2) * exp(2*ln_tau);
  jnll += GMRF(Q)(omega_s);

  // Probability of random effects (depending on if method is diffusion or not)
  if(method == "diffusion"){
    PARAMETER(ln_kappa2);
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;

    Eigen::SparseMatrix<Type> I_ss( depth_s.size(), depth_s.size() );
    I_ss.setIdentity();
    Eigen::SparseMatrix<Type> invD = I_ss + exp(-2.0 * ln_kappa2) * (invM0 * M1);

    lu.compute(invD);
    REPORT(invD);

    // Solve
    matrix<Type> depth_sz = lu.solve(depth_s.matrix());;

    depth_s = depth_sz.col(0); // turn one column matrix to vector

    // Derived range parameter for diffusion
    Type Range2 = sqrt(8) / exp(ln_kappa2);
    REPORT(Range2);
    REPORT(ln_kappa2);
  }
  
  // Projection to obs / grid
  vector<Type> depth_i = A_is * depth_s;
  vector<Type> depth_g = A_gs * depth_s;
  REPORT(depth_s);
  REPORT(depth_i);
  REPORT(depth_g);
  vector<Type> omega_i = A_is * omega_s;
  vector<Type> omega_g = A_gs * omega_s;
  
  // Polynomial depth effect (linear + quadratic terms)
  vector<Type> pdepth_i = depth_i*beta_j(0) + pow(depth_i,2)*beta_j(1);
  vector<Type> pdepth_g = depth_g*beta_j(0) + pow(depth_g,2)*beta_j(1);

  // Add weights here so that I can set them to 0 when calculating conditional AIC
  // Expected mean counts at observation locations
  vector<Type> mu_i = exp(beta0 + omega_i + pdepth_i);
  PARAMETER(ln_sigma_eta);
  PARAMETER_VECTOR(eta_i);
  for(int i=0; i<c_i.size(); i++){
    if (weights_i(i) > Type(0.0)) {
      jnll -= dnorm(eta_i(i), Type(0.0), exp(ln_sigma_eta), true);
      jnll -= dpois(c_i(i), mu_i(i) * exp(eta_i(i)), true);
    }
  }
  
  // Expected mean counts at prediction grid
  vector<Type> mu_g = exp(beta0 + omega_g + pdepth_g);

  // Reporting (output)
  REPORT(Q);
  REPORT(mu_i); // predicted counts at observations
  REPORT(mu_g); // predicted counts at grid level
  REPORT(pdepth_i);
  REPORT(pdepth_g);

  return jnll;
}

