functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>          n;                 // sample size of current data
  int<lower=0>          n0;                // sample size of historical data
  int<lower=0>          p;                 // number of covariates (incl intercept)
  vector[n]             y;                 // response vector for current data
  matrix[n,p]           X;                 // design mtx for current data
  vector[n0]            y0;                // response vector for historical data
  matrix[n0,p]          X0;                // design mtx for historical data
  vector[p]             coef_hp_mean;      // mean for normal hyperprior
  matrix[p,p]           coef_hp_cov;       // covariance for normal hyperprior
  vector[p]             coef_vague_mean;   // mean for vague prior on beta
  matrix[p,p]           coef_vague_cov;    // cov for vague prior on beta
  int<lower=0>          iw_hp_df;          // df for inv-wishart hyperprior
  matrix[p,p]           iw_hp_scale;       // scale mtx for inv-wishart hyperprior
  real<lower=0>         disp_shape;        // shape parameter for inv-gamma prior on dispersion
  real<lower=0>         disp_scale;        // scale parameter for inv-gamma prior on dispersion
  real<lower=0>         disp0_shape;       // shape parameter for inv-gamma prior on historical dispersion
  real<lower=0>         disp0_scale;       // scale parameter for inv-gamma prior on historical dispersion
  real<lower=0.000001,upper=0.999999> w;   // mixing parameter
  int<lower=1,upper=5>  dist;              // index for distribution
  int<lower=1,upper=9>  link;              // index for link
  vector[n]             offset;            // current data offset
  vector[n0]            offset0;           // historical data offset
}
transformed data{
  real logw   = log(w);
  real log1mw = log1m(w);
}
parameters {
  vector[p]         beta;
  vector[p]         beta0;
  vector[p]         coef_mean;
  cov_matrix[p]     coef_cov;
  real<lower=0>     dispersion[(dist > 2) ? 1 :  0];
  real<lower=0>     dispersion0[(dist > 2) ? 1 :  0];
}
model {
  coef_mean ~ multi_normal(coef_hp_mean, coef_hp_cov);         // hyperprior for mean of coefficients
  coef_cov  ~ inv_wishart(iw_hp_df, iw_hp_scale);              // hyperprior for cov of coefficients
  target += multi_normal_lpdf(beta0 | coef_mean, coef_cov);               // prior on beta0
  // robust MAP prior = w * hierarchical prior + (1 - w) * vague prior
  target   += log_sum_exp(
      logw   + multi_normal_lpdf(beta | coef_mean,       coef_cov)
    , log1mw + multi_normal_lpdf(beta | coef_vague_mean, coef_vague_cov)
  );
  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);    // current data likelihood
    target += glm_lp(y0, beta0, 1.0, X0, dist, link, offset0);   // historical data likelihood
  }
  else {
    dispersion  ~ inv_gamma(disp_shape,  disp_scale);
    dispersion0 ~ inv_gamma(disp0_shape, disp0_scale);
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);    // current data likelihood
    target += glm_lp(y0, beta0, dispersion0[1], X0, dist, link, offset0);   // historical data likelihood
  }
}

