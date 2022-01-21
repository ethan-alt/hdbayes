functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>         n;
  int<lower=0>         n0;
  int<lower=0>         p;
  vector[n]            y;
  matrix[n,p]          X;
  vector[n0]           y0;
  matrix[n0,p]         X0;
  vector[p]            hp_norm_mean;
  matrix[p,p]          hp_norm_cov;
  int<lower=0>         hp_iw_df;
  matrix[p,p]          hp_iw_scale;
  real<lower=0>        disp0_shape;
  real<lower=0>        disp0_scale;
  real<lower=0>        disp_shape;
  real<lower=0>        disp_scale;
  int<lower=1,upper=5> dist;
  int<lower=1,upper=9> link;
  vector[n]            offset;
  vector[n0]           offset0;
}
parameters {
  vector[p] beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
  vector[p] beta0;
  real<lower=0> dispersion0[(dist > 2) ? 1 :  0];
  vector[p]       beta_mean;
  cov_matrix[p]   beta_cov;
}
model {
  beta      ~ multi_normal(beta_mean, beta_cov);
  beta0     ~ multi_normal(beta_mean, beta_cov);
  beta_mean ~ multi_normal(hp_norm_mean, hp_norm_cov);
  beta_cov  ~ inv_wishart(hp_iw_df, hp_iw_scale);
  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);    // current data likelihood
    target += glm_lp(y0, beta0, 1.0, X0, dist, link, offset0);   // historical data likelihood
  }
  else {
    dispersion0 ~ inv_gamma(disp0_shape, disp0_scale);                    // prior for historical dispersion
    dispersion  ~ inv_gamma(disp_shape,  disp_scale);                     // prior for dispersion
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);  // current data likelihood
    target += glm_lp(y0, beta0, dispersion0[1], X0, dist, link, offset0); // historical data likelihood
  }
}

