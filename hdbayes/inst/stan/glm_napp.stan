functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>                        n;
  int<lower=0>                        p;
  int<lower=0,upper=1>                ind_disp;     // indicator for whether dispersion is included
  vector[n]                           y;
  matrix[n,p]                         X;
  vector[p + ind_disp]                theta_hat;   // MLE for beta and LOG DISPERSION
  matrix[p + ind_disp, p + ind_disp]  theta_cov;   // inverse fisher info for (beta, LOG PHI)
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[n]                           offset;
}
parameters {
  vector[p] beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
  real<lower=0,upper=1> a0;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // prior of a0 is beta (uniform if 1's)
  if ( a0_shape1 != 1 || a0_shape2 != 1 )
    a0 ~ beta(a0_shape1, a0_shape2);

  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);   // current data likelihood
    beta    ~ multi_normal(theta_hat, inv(a0) * theta_cov);     // asympt. power prior
    target += multi_normal_lpdf(beta | theta_hat, inv(a0) * theta_cov);
  }
  else {
    // change of variables:
    //   eta = (beta, phi) --> theta = (beta, log(phi))
    //   Jacobian       = blkdiag ( I, 1/phi )
    //   log |det(Jacobian)| = -log(phi)
    real log_disp      = log(dispersion[1]);
    vector[p+1] theta  = append_row( beta, log_disp );  // theta = ( beta, log(phi) )'
    target            += multi_normal_lpdf(theta | theta_hat, inv(a0) * theta_cov);
    target            += -log_disp;   // jacobian

    // likelihood of current data
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);
  }
}

