
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_IVECTOR( options_z );
  // options_z(0)=0: L1 is stationary  |  options_z(0)=1: L1 is innovations
  // options_z(1) NOT USED
  DATA_INTEGER( n_t );
  DATA_VECTOR( a_g );  // Area associated with grid-cell g
  DATA_VECTOR( b_i );  // response for observation i
  DATA_VECTOR( a_i );  // area-swept
  DATA_IVECTOR( t_i );  // Random effect index for observation i
  DATA_MATRIX( Temp_s2t );
  DATA_MATRIX( X_ik );
//  DATA_MATRIX( X_gk );

  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
//  DATA_SPARSE_MATRIX(A_gs);
  DATA_SPARSE_MATRIX(A_is2);

  // SPDE objects
  DATA_SPARSE_MATRIX(M0_ss);
  DATA_SPARSE_MATRIX(M1_ss);
  DATA_SPARSE_MATRIX(M2_ss);
  DATA_SPARSE_MATRIX(M1_s2s2);
  DATA_SPARSE_MATRIX(invM0_s2s2);
  DATA_SPARSE_MATRIX(M_tt);

  // Parameters
  PARAMETER_VECTOR( log_kappaS );
  PARAMETER_VECTOR( kappaT );
  PARAMETER_VECTOR( kappaST );
  PARAMETER_VECTOR( gamma_j );
  PARAMETER_VECTOR( beta_k );
  PARAMETER( ln_tauO );
  PARAMETER( ln_tauE );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi );
  PARAMETER( logit_rhoE );
  PARAMETER( finv_power );
  DATA_INTEGER(sim_gmrf); // simulate GMRFs?

  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_MATRIX( epsilon_st );

  // Objective funcction
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Globals
  int n_i = A_is.rows();
  int n_s2 = A_is2.cols();
  Type rhoE = invlogit( logit_rhoE );
  vector<Type> beta_i = X_ik * beta_k;
  Eigen::SparseMatrix<Type> I_s2s2( n_s2, n_s2 );
  Eigen::SparseMatrix<Type> I_tt( n_t, n_t );
  Eigen::SparseMatrix<Type> I_k2k2( n_s2*n_t, n_s2*n_t );
  I_s2s2.setIdentity();
  I_tt.setIdentity();
  I_k2k2.setIdentity();
  Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;

  // settings
  if( options_z(0) == 1 ){
    M_tt = M_tt - I_tt;
  }

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q_ss = exp(4*ln_kappa)*M0_ss + Type(2.0)*exp(2*ln_kappa)*M1_ss + M2_ss;
  jnll_comp(1) += SCALE( GMRF(Q_ss), 1/exp(ln_tauO) )( omega_s );
  if ( sim_gmrf ) SIMULATE {
    vector<Type> tmp(omega_s.size());
    GMRF(Q_ss).simulate( tmp );
    omega_s = tmp / exp(ln_tauO);
  }
  for( int t=0; t<n_t; t++ ){
    if( t==0 ){
      jnll_comp(2) += SCALE( GMRF(Q_ss), 1 / exp(ln_tauE) / pow( 1.0-pow(rhoE,2), 0.5 ) )( epsilon_st.col(t) );
      if ( sim_gmrf ) SIMULATE {
        vector<Type> tmp(epsilon_st.rows());
        GMRF(Q_ss).simulate( tmp );
        epsilon_st.col(t) = tmp / exp(ln_tauE) / pow( 1.0-pow(rhoE,2), 0.5 );
      }
    }else{
      jnll_comp(2) += SCALE( GMRF(Q_ss), 1 / exp(ln_tauE) )( epsilon_st.col(t) - rhoE*epsilon_st.col(t-1) );
      if ( sim_gmrf ) SIMULATE {
        vector<Type> tmp(epsilon_st.rows());
        GMRF(Q_ss).simulate( tmp );
        tmp = tmp / exp(ln_tauE) + rhoE * vector<Type>(epsilon_st.col(t-1));
        epsilon_st.col(t) = tmp;
      }
    }
  }

  // Spatio-temporal distributed lag
  Eigen::SparseMatrix<Type> P_s2s2( n_s2, n_s2 );
  Eigen::SparseMatrix<Type> P_k2k2( n_s2*n_t, n_s2*n_t );
  P_s2s2 = -1 * invM0_s2s2 * M1_s2s2;
  //P_kk.coeffRef( 0, 0 ) = 0;
  if( log_kappaS.size() > 0 ){
    P_k2k2 = P_k2k2 + exp(-2 * log_kappaS(0) ) * kronecker( I_tt, P_s2s2 );
  }
  if( kappaT.size() > 0 ){
    P_k2k2 = P_k2k2 + kappaT(0) * kronecker( M_tt, I_s2s2 );
  }
  if( kappaST.size() > 0 ){
    P_k2k2 = P_k2k2 + kappaST(0) * exp(-2 * log_kappaS(0) ) * kronecker( M_tt, P_s2s2 );
  }

  // Solve and repack
  Eigen::SparseMatrix<Type> IminusP_k2k2( n_s2*n_t, n_s2*n_t );
  IminusP_k2k2 = I_k2k2 - P_k2k2;
  lu.compute(IminusP_k2k2);
  matrix<Type> Temp_k2 = Temp_s2t.reshaped( n_s2*n_t, 1 );
  matrix<Type> T_k2 = lu.solve(Temp_k2);
  matrix<Type> T_s2t = T_k2.reshaped( n_s2, n_t );
  REPORT( T_s2t );

  // Spatial diffusion
//  Eigen::SparseMatrix<Type> IminusP_ss( n_s, n_s );
//  IminusP_ss = I_ss - P_ss;
//  lu.compute(IminusP_ss);
//  matrix<Type> T_st = lu.solve(Temp_st);;
//  REPORT( T_st )
//  REPORT( P_ss )

  // True density and abundance
//  vector<Type> omega_g( n_g );
//  matrix<Type> epsilon_gt( n_g, n_t );
//  omega_g = A_gs * omega_s;
//  epsilon_gt = A_gs * epsilon_st;
//  array<Type> ln_d_gt( n_g, n_t );
//  for( int t=0; t<n_t; t++){
//  for( int g=0; g<n_g; g++){
//    ln_d_gt(g,t) = beta_t(t) + omega_g(g) + epsilon_gt(g,t);
//  }}

  // Probability of data conditional on random effects
  vector<Type> omega_i( n_i );
  matrix<Type> epsilon_it( n_i, n_t );
  matrix<Type> T_it( n_i, n_t );
  omega_i = A_is * omega_s;
  epsilon_it = A_is * epsilon_st;
  T_it = A_is2 * T_s2t;
  vector<Type> bhat_i( n_i );
  vector<Type> p_i( n_i );
  vector<Type> gamma_i( n_i );
  for( int i=0; i<n_i; i++){
    gamma_i(i) = gamma_j(0) * T_it(i,t_i(i));
    if( gamma_j.size() > 1 ){
      gamma_i(i) += gamma_j(1) * pow( T_it(i,t_i(i)), 2 );
    }
    p_i(i) = beta_i(i) + omega_i(i) + epsilon_it(i,t_i(i)) + gamma_i(i);
    //if( !isNA(b_i(i)) ){
    bhat_i(i) = a_i(i) * exp( p_i(i) );
    jnll_comp(0) -= dtweedie( b_i(i), bhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true );
    SIMULATE{b_i(i) = rtweedie(bhat_i(i), exp(ln_phi), 1.0 + invlogit(finv_power));}
    //}
  }

  // Objective function
  Type jnll = jnll_comp.sum();

//  // Derived quantities
//  vector<Type> b_t( n_t );
//  b_t.setZero();
//  vector<Type> zmean_t( n_t );
//  zmean_t.setZero();
//  for( int t=0; t<n_t; t++){
//    for( int g=0; g<n_g; g++){
//      b_t(t) += a_g(g) * exp(ln_d_gt(g,t));
//    }
//  }

  // Derived quantities
  Type rhoT = 0;
  if( kappaT.size() > 0 ){
    if( options_z(0) == 0 ){
      rhoT = kappaT(0);
    }else{
      rhoT = kappaT(0) / (1.0 + kappaT(0));
    }
    REPORT( rhoT );
    ADREPORT( rhoT );
  }
  if( log_kappaS.size() > 0 ){
    Type MSD;
    if( (options_z(0) == 1) && (kappaT.size() > 0) ){
      MSD = 4 / exp( 2.0 * log_kappaS(0) ) * (1 - rhoT);
    }else{
      MSD = 4 / exp( 2.0 * log_kappaS(0) );
    }
    Type RMSD = pow( MSD, 0.5 );
    REPORT( MSD );
    REPORT( RMSD );
    ADREPORT( MSD );
    ADREPORT( RMSD );
  }

  // Reporting
  REPORT( jnll_comp );
  REPORT( jnll );
//  REPORT( ln_d_gt );
  REPORT( bhat_i );
  REPORT( p_i );
  REPORT( gamma_i );
  REPORT( IminusP_k2k2 );
  ADREPORT( gamma_j );  // Get covariance to sample response curve
  //REPORT( Range );
  //REPORT( SigmaE );
  //REPORT( SigmaO );
//  REPORT( b_t );
//  ADREPORT( b_t );

  SIMULATE {
    REPORT(b_i);
  }

  return jnll;
}
