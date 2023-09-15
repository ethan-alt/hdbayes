functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>                        K; // total number of historical datasets (excluding the current data)
  int<lower=0>                        n; // number of observations in the current data
  int<lower=0>                        p;
  int<lower=0,upper=1>                ind_disp; // indicator for whether dispersion is included
  vector[n]                           y; // response for the current data
  matrix[n,p]                         X; // design matrix for the current data
  matrix[p + ind_disp, K]             theta_hats; // the kth column consists of the MLE for beta and LOG DISPERSION for the kth dataset (the first dataset is the current data)
  array[K] cov_matrix[p + ind_disp]   theta_covars; // the kth element is the inverse fisher info for (beta, LOG PHI) for the kth dataset (the first dataset is the current data)
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[n]                           offs; // offset for current data


}
parameters {
  vector[p] beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
  vector<lower=0,upper=1>[K] a0s;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // prior of a0s is beta (uniform if 1's)
  if ( a0_shape1 != 1 || a0_shape2 != 1)
    a0s ~ beta(a0_shape1, a0_shape2);

  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offs); // current data likelihood

    // asympt. power prior
    for ( k in 1:K ) {
      target += multi_normal_lpdf(beta | theta_hats[, k], inv(a0s[k]) * theta_covars[k]);
    }
  }
  else {
    // change of variables:
    //   eta = (beta, phi) --> theta = (beta, log(phi))
    //   Jacobian       = blkdiag ( I, 1/phi )
    //   log |det(Jacobian)| = -log(phi)
    real log_disp      = log(dispersion[1]);
    vector[p+1] theta  = append_row( beta, log_disp );  // theta = ( beta, log(phi) )'
    for ( k in 1:K ) {
      target          += multi_normal_lpdf(theta | theta_hats[, k], inv(a0s[k]) * theta_covars[k]);
    }
    target            += -log_disp;   // jacobian

    // likelihood of current data
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offs);
  }
}

