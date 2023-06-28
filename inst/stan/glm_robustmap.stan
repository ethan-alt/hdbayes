functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>           n1;                 // sample size of current data
  int<lower=0>           p;                  // number of covariates (including the intercept)
  vector[n1]             y1;                 // response vector for current data
  matrix[n1,p]           X1;                 // design matrix for current data
  int<lower=0>           G;                  // number of mixture components in the mixture approximation to the prior induced by the BHM
  simplex[G]             probs;              // mixing proportions for the robust MAP prior
  matrix[p, G]           means;              // mean matrix for the mixture approximation to the prior induced by the BHM
  array[G] cov_matrix[p] covars;             // the gth element is the covariance matrix for the gth component of the mixture approximation
  vector[p]              norm_vague_mean;    // mean for vague prior on beta
  vector<lower=0>[p]     norm_vague_sd;      // sd for vague prior on beta
  real                   disp_mean1;         // mean for the half-normal prior for dispersion of current data
  real<lower=0>          disp_sd1;           // sd for the half-normal prior for dispersion of current data
  real<lower=0,upper=1>  w;                  // mixing parameter
  int<lower=1,upper=5>   dist;               // index for distribution
  int<lower=1,upper=9>   link;               // index for link
  vector[n1]             offs1;              // current data offset
}
transformed data{
  array[G] cov_matrix[p] precisions;
  for ( i in 1:G )
    precisions[i] = inverse_spd(covars[i]);
}
parameters {
  vector[p]                            beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
}
model {
  // robust MAP prior = w * mixture approximation to the prior induced by the BHM + (1 - w) * vague prior
  target += log_mix(
    w
    , multi_normal_mix_lpdf(beta | probs, means, precisions)
    , normal_lpdf(beta | norm_vague_mean, norm_vague_sd)
  );

  if ( dist <= 2 ) {
    target += glm_lp(y1, beta, 1.0, X1, dist, link, offs1);    // current data likelihood
  }
  else {
    dispersion ~ normal(disp_mean1, disp_sd1); // half-normal prior for dispersion
    target += glm_lp(y1, beta, dispersion[1], X1, dist, link, offs1);    // current data likelihood
  }
}

