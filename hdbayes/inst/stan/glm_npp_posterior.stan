functions{
#include include/expfam_loglik.stan
#include include/finda0closest.stan
}
data {
  int<lower=0>         n;
  int<lower=0>         n0;
  int<lower=0>         p;
  int<lower=0>         k;          // number of a0s for which we have log nc
  vector[n]            y;
  matrix[n,p]          X;
  vector[n0]           y0;
  matrix[n0,p]         X0;
  vector[p]            beta_mean;
  matrix[p,p]          beta_cov;
  real<lower=0>        disp_shape;
  real<lower=0>        disp_scale;
  vector<lower=0,upper=1>[k] a0_vec;
  vector[k]            lognca0_vec;
  real<lower=0>        a0_shape1;
  real<lower=0>        a0_shape2;
  int<lower=1,upper=5> dist;
  int<lower=1,upper=9> link;
  vector[n]            offset;
  vector[n0]           offset0;
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
  real<lower=0,upper=1> a0;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  if ( a0_shape1 != 1 || a0_shape2 != 1)
    a0 ~ beta(a0_shape1, a0_shape2);
  beta ~ multi_normal(beta_mean, beta_cov);                                 // prior for beta
  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);         // current data likelihood
    target += a0 * glm_lp(y0, beta, 1.0, X0, dist, link, offset0);   // historical data likelihood
  }
  else {
    dispersion  ~ inv_gamma(disp_shape,  disp_scale);                         // prior for dispersion
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);      // current data likelihood
    target += a0 * glm_lp(y0, beta, dispersion[1], X0, dist, link, offset0);  // historical data likelihood
  }
  ## Subtract log nc from power prior
  target += -pp_lognc(a0, a0_vec, lognca0_vec);
}

