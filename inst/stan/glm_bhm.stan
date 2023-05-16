functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>                       K; // total number of datasets (including the current data)
  int<lower=0>                       N; // total number of observations (including the current data)
  array[K] int<lower=0, upper=N>     start_idx; // starting index of each data in the stacked version
  array[K] int<lower=0, upper=N>     end_idx; // ending index of each data in the stacked version
  int<lower=0>                       p;
  vector[N]                          y; // response for the stacked data
  matrix[N,p]                        X; // design matrix for the stacked data
  vector[p]                          meta_mean_mean; // mean for the hyperprior for mean of each coefficient
  array[p] real<lower=0>             meta_mean_sd; // sd for the hyperprior for mean of each coefficient
  vector[p]                          meta_sd_mean; // mean for the hyperprior for sd of each coefficient
  array[p] real<lower=0>             meta_sd_sd; // sd for the hyperprior for sd of each coefficient
  vector[K]                          disp_mean; // mean for the half-normal prior for dispersion (same for all datasets)
  array[K] real<lower=0>             disp_sd; // sd for the half-normal prior for dispersion (same for all datasets)
  int<lower=1,upper=5>               dist;
  int<lower=1,upper=9>               link;
  vector[N]                          offs; // offset
}
parameters {
  matrix[p, K]           beta; // kth column is the vector of coefficients for the kth dataset
  vector[p]              beta_mean;
  array[p] real<lower=0> beta_sd;
  array[(dist > 2) ? K :  0] real<lower=0> dispersion;
}
model {
  for ( j in 1:p ) {
    beta_mean[j] ~ normal(meta_mean_mean[j], meta_mean_sd[j]); // hyperprior for mean of coefficients
    beta_sd[j]   ~ normal(meta_sd_mean[j], meta_sd_sd[j]); // hyperprior for sd of coefficients (half-normal)
    beta[j, ]    ~ normal(beta_mean[j], beta_sd[j]);
  }

  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
  else {
    for ( k in 1:K ) {
      dispersion[k] ~ normal(disp_mean[k], disp_sd[k]); // half-normal prior for dispersion
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], dispersion[k], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
}

