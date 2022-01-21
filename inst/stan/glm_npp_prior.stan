functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>          n0;
  int<lower=0>          p;
  vector[n0]            y0;
  matrix[n0,p]          X0;
  vector[p]             beta_mean;
  matrix[p,p]           beta_cov;
  real<lower=0>         disp_shape;
  real<lower=0>         disp_scale;
  int<lower=1,upper=5>  dist;
  int<lower=1,upper=9>  link;
  vector[n0]            offset0;
  real<lower=0,upper=1> a0;
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  beta ~ multi_normal(beta_mean, beta_cov);
  if ( dist <= 2 ) {
    target += a0 * glm_lp(y0, beta, 1.0, X0, dist, link, offset0);   // historical data likelihood
  }
  else {
    dispersion  ~ inv_gamma(disp_shape,  disp_scale);                        // prior for dispersion
    target += a0 * glm_lp(y0, beta, dispersion[1], X0, dist, link, offset0); // historical data likelihood
  }
}

