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
  vector[p]             mean_beta;         // mean for normal initial prior on coefficients
  matrix[p,p]           cov_beta;          // covariance mtx for normal initial prior on coefficients
  real<lower=0>         disp_shape;        // shape parameter for inv-gamma initial prior on dispersion
  real<lower=0>         disp_scale;        // scale parameter for inv-gamma initial prior on dispersion
  real<lower=0,upper=1> a0;                // power prior parameter (fixed)
  int<lower=1,upper=5>  dist;
  int<lower=1,upper=9>  link;
  vector[n]             offset;
  vector[n0]            offset0;
}
parameters {
  vector[p]     beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
}
model {
  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);       // current data likelihood
    target += a0 * glm_lp(y0, beta,  1.0, X0, dist, link, offset0); // power prior
    beta    ~ multi_normal(mean_beta, cov_beta);                    // initial prior beta ~ N(mu0, cov0)
  }
  else {
    target     += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);       // current data likelihood
    target     += a0 * glm_lp(y0, beta,  dispersion[1],  X0, dist, link, offset0); // historical data likelihood
    beta        ~ multi_normal(mean_beta, dispersion[1] * cov_beta);               // initial prior beta | disp ~ N(mu0, disp * cov0)
    dispersion  ~ inv_gamma(disp_shape, disp_scale);                                // initial prior disp ~ inv_gamma(shape0, scale0)
  }
}

