namespace sdmTMB {
template <class Type>
bool isNA(Type x) {
  return R_IsNA(asDouble(x));
}

// dgengamma
// Written by J. Thorson based on scripts listed below
// using Prentice-1974 parameterization for lambda instead of k, so that lognormal occurs as lambda -> 0
// using mean parameterization to back out theta
// CV is a function of sigma and lambda and NOT mean (i.e., CV is fixed for all values of mean)
// See: C:\Users\James.Thorson\Desktop\Work files\AFSC\2021-10 -- Generalized gamma-lognormal\Explore gengamma.R
template<class Type>
Type dgengamma( Type x,
                Type mean,
                Type sigma,
                Type Q,
                int give_log=0){


  // First convert from mean to mu
  // Rename based on flexdist
  Type lambda = Q;
  // Using https://stats.stackexchange.com/questions/345564/generalized-gamma-log-normal-as-limiting-special-case
  Type k = pow( lambda, -2);
  Type beta = pow( sigma, -1) * lambda;
  // Use wikipedia expression for the mean:   mean = a * gamma((d+1)/p) / gamma(d/p) where d/p = k, theta = a, and beta = p
  // https://en.wikipedia.org/wiki/Generalized_gamma_distribution#Software_implementation
  Type log_theta = log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k );
  // Using https://stats.stackexchange.com/questions/345564/generalized-gamma-log-normal-as-limiting-special-case
  Type mu = log_theta + log(k) / beta;

  // Next evaluate PDF
  // from https://github.com/chjackson/flexsurv-dev/blob/master/src/gengamma.h#L54-L56
  Type y = log(x);
  Type w = (y - mu) / sigma;
  Type qi = pow(Q, -2);
  Type qw = Q * w;                 // 0.5*log(pow(x,2)) as trick for abs(log(x))
  Type logres = -log(sigma*x) + 0.5*log(pow(lambda,2)) * (1 - 2 * qi) + qi * (qw - exp(qw)) - lgamma(qi);

  // return stuff
  if(give_log) return logres; else return exp(logres);
}

// Ordered beta regression (Kubinec 2023) log density.
// eta is the linear predictor, mu = invlogit(eta), psi0 < psi1 are cutpoints
// on the logit scale, phi is the beta precision. Mirrors glmmTMB ordbeta_family.
template<class Type>
Type dordbeta(Type y, Type eta, Type mu, Type phi, Type psi0, Type psi1,
              int give_log = 1) {
  Type logres;
  if (y == Type(0)) {
    // log(1 - invlogit(eta - psi0))
    logres = -logspace_add(Type(0), eta - psi0);
  } else if (y == Type(1)) {
    // log(invlogit(eta - psi1))
    logres = -logspace_add(Type(0), -(eta - psi1));
  } else {
    Type log_F0 = -logspace_add(Type(0), -(eta - psi0)); // log invlogit(eta - psi0)
    Type log_F1 = -logspace_add(Type(0), -(eta - psi1)); // log invlogit(eta - psi1)
    Type log_pmid = logspace_sub(log_F0, log_F1);
    Type s1 = mu * phi;
    Type s2 = (Type(1) - mu) * phi;
    logres = log_pmid + dbeta(y, s1, s2, true);
  }
  if (give_log) return logres; else return exp(logres);
}

// Simulate one ordered beta draw given linear predictor and parameters.
template<class Type>
Type rordbeta(Type eta, Type mu, Type phi, Type psi0, Type psi1) {
  // One uniform selects the component: P(0) = p0, P(1) = p1.
  Type u = runif(Type(0), Type(1));
  if (u < invlogit(psi0 - eta)) return Type(0);
  if (u > Type(1) - invlogit(eta - psi1)) return Type(1);
  Type s1 = mu * phi;
  Type s2 = (Type(1) - mu) * phi;
  return rbeta(s1, s2);
}

// rgengamma
// Written by J. Thorson based on scripts listed below
template<class Type>
Type rgengamma( Type mean,
                Type sigma,
                Type Q){

  // See: C:\Users\James.Thorson\Desktop\Work files\AFSC\2021-10 -- Generalized gamma-lognormal\Explore gengamma.R
  Type lambda = Q;
  Type k = pow( lambda, -2 );
  Type beta = pow( sigma, -1 ) * lambda;
  Type log_theta = log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k );
  Type w = log(rgamma(k, Type(1.0)));
  Type y = w/beta + log_theta;
  return exp(y);
}

extern "C" {
  /* Rmath entry points (Writing R Extensions: 'Numerical analysis subroutines') */
  double Rf_dpois(double x, double lambda, int give_log);
  double Rf_logspace_sub(double logx, double logy);
}

// log P(L <= Y <= U) for Y ~ Poisson(lambda) with fixed integer bounds;
// `U = Inf` is right censoring. Values are computed on the log scale from the
// tail that avoids cancellation. The branches depend on `lambda` but run only
// in the double-valued atomic below, which is evaluated anew at every lambda.
inline double censpois_logprob(double lambda, double L, double U) {
  using atomic::Rmath::Rf_ppois;
  if (!std::isfinite(U)) {
    if (L <= 0) return 0.0;
    return Rf_ppois(L - 1, lambda, 0, 1); // log P(Y > L - 1)
  }
  if (L <= 0) return Rf_ppois(U, lambda, 1, 1);
  double ans = R_NegInf;
  if (U - L > 64) {
    if (lambda < L) { // mass above the interval: difference of upper tails
      ans = Rf_logspace_sub(Rf_ppois(L - 1, lambda, 0, 1),
        Rf_ppois(U, lambda, 0, 1));
    } else if (lambda > U) { // mass below: difference of lower CDFs
      ans = Rf_logspace_sub(Rf_ppois(U, lambda, 1, 1),
        Rf_ppois(L - 1, lambda, 1, 1));
    } else { // lambda inside: 1 minus both tails
      ans = log1p(-(exp(Rf_ppois(L - 1, lambda, 1, 1)) +
        exp(Rf_ppois(U, lambda, 0, 1))));
    }
  }
  if (!std::isfinite(ans)) { // narrow interval (or fallback): sum the PMF
    ans = Rf_dpois(L, lambda, 1);
    for (double k = L + 1; k <= U; k++)
      ans = logspace_add(ans, Rf_dpois(k, lambda, 1));
  }
  return ans;
}

// Atomic with inputs (lambda, L, U); the bounds are data. Using the value y,
// d y / d lambda = [p(L - 1) - p(U)] / P = exp(log p(L - 1) - y) -
// exp(log p(U) - y), with the terms omitted for L = 0 or U = Inf. The rule
// is written with AD types so higher-order derivatives work.
TMB_ATOMIC_VECTOR_FUNCTION(
  censpois_logprob
  ,
  1
  ,
  ty[0] = censpois_logprob(tx[0], tx[1], tx[2]);
  ,
  Type lambda = tx[0];
  double L = asDouble(tx[1]);
  double U = asDouble(tx[2]);
  Type d = Type(0);
  if (L > 0) d += exp(dpois(Type(L - 1), lambda, true) - ty[0]);
  if (std::isfinite(U)) d -= exp(dpois(Type(U), lambda, true) - ty[0]);
  px[0] = d * py[0];
  px[1] = Type(0);
  px[2] = Type(0);
)

template <class Type>
Type censpois_logprob(Type lambda, Type L, Type U) {
  CppAD::vector<Type> tx(3);
  tx[0] = lambda;
  tx[1] = L;
  tx[2] = U;
  return censpois_logprob(tx)[0];
}

// Right-censored (`upr` NA: count >= x), interval-censored
// (x <= count <= upr), or exact (upr == x) Poisson.
template <class Type>
Type dcenspois2(Type x, Type lambda, Type upr, int give_log = 0) {
  Type ll;
  if (!isNA(upr) && upr == x) {
    ll = dpois(x, lambda, true);
  } else {
    ll = censpois_logprob(lambda, x, isNA(upr) ? Type(R_PosInf) : upr);
  }
  if (give_log) {
    return ll;
  } else {
    return exp(ll);
  }
}

template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0) {
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

// List of matrices
template <class Type>
struct LOM_t : vector<matrix<Type> > {
  LOM_t(SEXP x) {  // x = list passed from R
    (*this).resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

// Function to import barrier-SPDE code
// From Olav Nikolai Breivik and Hans Skaug:
template <class Type>
struct spde_barrier_t {
  vector<Type> C0;
  vector<Type> C1;
  Eigen::SparseMatrix<Type> D0;
  Eigen::SparseMatrix<Type> D1;
  Eigen::SparseMatrix<Type> I;
  spde_barrier_t(SEXP x) { /* x = List passed from R */
    C0 = asVector<Type>(getListElement(x, "C0"));
    C1 = asVector<Type>(getListElement(x, "C1"));
    D0 = tmbutils::asSparseMatrix<Type>(getListElement(x, "D0"));
    D1 = tmbutils::asSparseMatrix<Type>(getListElement(x, "D1"));
    I = tmbutils::asSparseMatrix<Type>(getListElement(x, "I"));
  }
};

template<class Type>
Eigen::SparseMatrix<Type> Q_SAR(Type rho, const Eigen::SparseMatrix<Type>& W) {
  int n_s = W.rows();
  Eigen::SparseMatrix<Type> I(n_s, n_s);
  I.setIdentity();
  Eigen::SparseMatrix<Type> B = I - rho * W;
  return B.transpose() * B;
}

template<class Type>
Eigen::SparseMatrix<Type> Q_CAR(Type alpha, const Eigen::SparseMatrix<Type>& W) {
  int n_s = W.rows();
  vector<Type> one(n_s);
  one.setOnes();
  vector<Type> D = W * one;
  Eigen::SparseMatrix<Type> Dmat(n_s, n_s);
  for (int i = 0; i < n_s; ++i) {
    Dmat.coeffRef(i, i) = D(i) > Type(0) ? D(i) : Type(1);
  }
  return Dmat - alpha * W;
}

// Alternative barrier SPDE function matching INLAspacetime implementation
// This version uses INLAspacetime's mathematical formulation but lets sdmTMB
// handle variance scaling through its existing tau mechanism
template <class Type>
Eigen::SparseMatrix<Type> Q_spde_inlaspacetime(spde_barrier_t<Type> spde, 
                                              Type ln_kappa,
                                              vector<Type> barrier_scaling) {
  int n = spde.I.rows();

  // Convert from sdmTMB parameterization 
  Type kappa = exp(ln_kappa);
  Type range = sqrt(Type(8)) / kappa;  // sdmTMB range formula

  // Combine domains as in INLAspacetime:
  // CC = C[[1]] + C[[2]] * range_fraction^2
  // Dmat = D[[1]] + D[[2]] * range_fraction^2
  // Use barrier_scaling(1) as the range_fraction (like sdmTMB's c(1))
  Type range_fraction = barrier_scaling.size() > 1 ? barrier_scaling(1) : Type(0.1);
  vector<Type> CC = spde.C0 + spde.C1 * pow(range_fraction, 2);
  Eigen::SparseMatrix<Type> Dmat = spde.D0 + spde.D1 * pow(range_fraction, 2);

  // Create inverse C diagonal matrix: iC = Diagonal(n, 1 / CC)
  Eigen::SparseMatrix<Type> iC(n, n);
  for (int i = 0; i < n; ++i) {
    iC.coeffRef(i, i) = Type(1) / CC(i);
  }

  // Construct the four combined matrices as in INLAspacetime:
  Eigen::SparseMatrix<Type> ICI = spde.I.transpose() * iC * spde.I;
  Eigen::SparseMatrix<Type> ICD = spde.I.transpose() * iC * Dmat;
  Eigen::SparseMatrix<Type> DCI = Dmat.transpose() * iC * spde.I;
  Eigen::SparseMatrix<Type> DCD = Dmat.transpose() * iC * Dmat;

  // Apply INLAspacetime parameter scaling
  // We use sigma = 1 here since we'll handle the proper scaling later in GMRF calls
  Type range2 = pow(range, 2);
  Type sigma2 = Type(1);  // Unit sigma - actual scaling handled in GMRF
  Type pi2s2 = Type(2) / (M_PI * sigma2);  // 2/(π*σ²) = 2/π with unit sigma

  Type param0 = pi2s2 / range2;              // 2/(π*σ²*r²)
  Type param1 = pi2s2 / Type(8);             // 2/(π*σ²*8)
  Type param2 = param1;                      // 2/(π*σ²*8)  
  Type param3 = range2 * pi2s2 / Type(64);   // 2*r²/(π*σ²*64)

  // Construct precision matrix using INLAspacetime's structure
  Eigen::SparseMatrix<Type> Q = param0 * ICI + param1 * ICD + param2 * DCI + param3 * DCD;

  return Q;
}

template <class Type>
Type barrier_scaling_factor(Type ln_tau, Type ln_kappa) {
  Type tau = exp(ln_tau);
  Type kappa = exp(ln_kappa);
  // marginal std. dev. for SPDE in 2D (alpha = 2)
  return Type(1) / (tau * kappa * sqrt(Type(4.0) * M_PI));
}

template <class Type>
Type minus_one_to_one(Type x) {
  return Type(2) * invlogit(x) - Type(1);
}

template <class Type>
Type log_sum_exp(Type x1, Type x2) {
  Type xmax = x1;
  if (x2 > x1) xmax = x2;
  return xmax + log(exp(x1 - xmax) + exp(x2 - xmax));
}

template <class Type>
matrix<Type> MakeH(vector<Type> x) {
  matrix<Type> H(2, 2);
  H(0, 0) = exp(x(0));
  H(1, 0) = x(1);
  H(0, 1) = x(1);
  H(1, 1) = (1 + x(1) * x(1)) / exp(x(0));
  return H;
}

// https://github.com/hrue/r-inla/blob/devel/r-inla.org/doc/prior/pc.matern.pdf
template <class Type>
Type pc_prior_matern(Type logtau, Type logkappa, Type matern_range,
                     Type matern_SD, Type range_prob, Type SD_prob,
                     int give_log = 0, int share_range = 0, int stan_flag = 0) {
  Type d = 2.;  // dimension
  Type dhalf = d / 2.;
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(SD_prob) / matern_SD;
  Type range = sqrt(8.) / exp(logkappa);
  Type sigma = 1. / sqrt(4. * M_PI * exp(2. * logtau) * exp(2. * logkappa));
  Type range_ll = log(dhalf) + log(lam1) + log(pow(range, -1. - dhalf)) -
                  lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll;
  if (!share_range) penalty += range_ll;

  // Note: these signs are + (and different from inst/jacobian-pcprior-tests)
  // because the jnll is accumulated
  if (stan_flag) {
    penalty += log(sqrt(8.)) - log(pow(range, 2.)); // P(sigma)
    Type C = sqrt(exp(lgamma(1. + dhalf)) * pow(4. * M_PI, dhalf));
    penalty += log(C) + logkappa;
  }
  // std::cout << "PC penalty: " << penalty << "\n";
  if (give_log)
    return penalty;
  else
    return exp(penalty);
}

template <class Type>
Type linear_threshold(Type x, Type slope, Type cutpoint) {
  // linear threshold model. relationship linear up to a point then constant
  // keep all parameters unconstrained - slope and scale can be neg/pos,
  // as can cutpoint if covariate is scaled ~ N(0,1).
  return slope * CppAD::CondExpLt(x, cutpoint, x, cutpoint);
}

template <class Type>
Type logistic_threshold(Type x, Type s50, Type s95, Type scale) {
  // logistic threshold model. similar to length or size based selectivity
  // in fisheries, parameterized by the points at which f(x) = 0.5 or 0.95
  // s50 and scale are unconstrained. s95 has to be > s50 though, so modelled as
  // s95 = s50 + exp(b(1))
  // Type s95 = s50 + exp(soffset); // this done outside function
  Type pred = (scale)*Type(1.0) /
              (Type(1.0) + exp(-log(Type(19.0)) * (x - s50) / (s95 - s50)));
  return pred;
}

// from glmmTMB distrib.h
// alpha = size (dispersion param), k = truncation point, mu = mean
double rtruncated_nbinom(double alpha, int k, double mu) {
  int m;
  double mdoub;
  double p = alpha / (mu + alpha);
  double q = mu / (mu + alpha);

  if (alpha <= 0.0)
    throw std::range_error("non-positive size in k-truncated-neg-bin simulator\n");
  if (mu <= 0.0)
    throw std::range_error("non-positive mu in k-truncated-neg-bin simulator\n");
  if (k < 0)
    throw std::range_error("negative k in k-truncated-neg-bin simulator\n");

  mdoub = (k + 1.0) * p - alpha * q;
  if (mdoub < 0.0)
    mdoub = 0.0;
  m = mdoub;
  if (m < mdoub)
    m = m + 1;
  /* since p < 1.0 and q > 0.0 we have 0.0 <= mdoub < k + 1
     hence 0 <= m <= k + 1 */

  for (;;) {
    double x = rnbinom(alpha + m, p) + m;
    if (m > 0) {
      double a = 1.0;
      int j;
      double u = unif_rand();
      for (j = 0; j < m; ++j)
        a *= (k + 1 - j) / (x - j);
      if (u < a && x > k)
        return x;
    } else {
      if (x > k)
        return x;
    }
  }
}

template<class Type>
Type calc_rf_sigma(Type ln_tau, Type ln_kappa) {
  Type sigma = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_tau + Type(2.) * ln_kappa));
  return sigma;
}

// from glmmTMB distrib.h
extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  double Rf_logspace_sub (double logx, double logy);
  void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
}

// from glmmTMB distrib.h
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  logit_invcloglog
  ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
,
  // ATOMIC_REVERSE
  px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
)

template<class Type>
Type logit_invcloglog(Type x) {
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return logit_invcloglog(tx)[0];
}

// get sign of double, only for REPORT use, from tinyVAST
template<class Type>
Type sign(Type x){
  if (x == Type(0.0)) return Type(0.0);
  return (x > Type(0.0)) ? Type(1.0) : Type(-1.0);
}

// Deviance for the Tweedie
// https://en.wikipedia.org/wiki/Tweedie_distribution#Properties
// From tinyVAST:
template<class Type>
Type devresid_tweedie( Type y,
                       Type mu,
                       Type p ){
  Type c1 = pow( y, 2.0-p ) / (1.0-p) / (2.0-p);
  Type c2 = y * pow( mu, 1.0-p ) / (1.0-p);
  Type c3 = pow( mu, 2.0-p ) / (2.0-p);
  Type deviance = 2 * (c1 - c2 + c3 );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
}

// Binomial deviance residual for y successes out of `size` trials.
template<class Type>
Type devresid_binomial(Type y, Type size, Type logit_p) {
  Type log_p = -logspace_add(Type(0), -logit_p);
  Type log_one_minus_p = -logspace_add(Type(0), logit_p);
  Type deviance = 0;
  if (y > 0) deviance += y * (log(y / size) - log_p);
  if (y < size) deviance += (size - y) * (log((size - y) / size) - log_one_minus_p);
  return sign(y - size * exp(log_p)) * pow(Type(2) * deviance, 0.5);
}

// From tinyVAST:
template<class Type>
Type devresid_nbinom2( Type y,
                       Type logmu,
                       Type logtheta ){
  Type logp1 = dnbinom_robust( y, log(y + Type(1e-10)), Type(2.0) * log(y + Type(1e-10)) - logtheta, true );
  Type logp2 = dnbinom_robust( y, logmu, Type(2.0) * logmu - logtheta, true );
  Type deviance = 2 * (logp1 - logp2);
  Type devresid = sign( y - exp(logmu) ) * pow( deviance, 0.5 );
  return devresid;
}

// Beta-binomial distribution
// Modified from glmmTMB
template<class Type>
Type dbetabinom_robust(Type y, Type loga, Type logb, Type n, int give_log=0)
{
  Type logres =
    lgamma(n + Type(1)) - lgamma(y + Type(1)) - lgamma(n - y + Type(1)) +
    lgamma(exp(loga) + exp(logb)) + lgamma(y + exp(loga)) + lgamma(n - y + exp(logb)) -
    lgamma(n + exp(loga) + exp(logb)) - lgamma(exp(loga)) - lgamma(exp(logb));
  if(!give_log) return exp(logres);
  else return logres;
}

// Helper for low-rank sparse Hessian bias correction
// Used by derived quantities (index, cog, weighted_avg, eao)
// See Kristensen et al. (2016) doi:10.18637/jss.v070.i05 and
// TMB's newton::Tag for low-rank sparse Hessian approximation
template <class Type>
Type add_lowrank_bias_correction(
    const vector<Type>& values,
    const vector<Type>& eps_index,
    const vector<int>& include,
    int offset,
    Type jnll
) {
  if (eps_index.size() > 0) {
    for (int t = 0; t < values.size(); t++) {
      if (include.size() == 0 || include(t) != 0) {
        Type S = newton::Tag(values(t));
        jnll += eps_index(t + offset) * S;
      }
    }
  }
  return jnll;
}

}  // namespace sdmTMB
