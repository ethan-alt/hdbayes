functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>                        K; // total number of datasets (including the current data)
  int<lower=0>                        N; // total number of observations (including the current data)
  array[K] int<lower=0, upper=N>      start_idx; // starting index of each data in the stacked version
  array[K] int<lower=0, upper=N>      end_idx; // ending index of each data in the stacked version
  int<lower=0>                        p;
  vector[N]                           y; // response for the stacked data
  matrix[N,p]                         X; // design matrix for the stacked data
  vector[p]                           meta_mean_mean; // mean for the hyperprior for mean of each coefficient
  vector<lower=0>[p]                  meta_mean_sd; // sd for the hyperprior for mean of each coefficient
  vector[p]                           meta_sd_mean; // mean for the hyperprior for sd of each coefficient
  vector<lower=0>[p]                  meta_sd_sd; // sd for the hyperprior for sd of each coefficient
  vector[K]                           disp_mean; // mean for the half-normal prior for dispersion of each dataset
  vector<lower=0>[K]                  disp_sd; // sd for the half-normal prior for dispersion of each dataset
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
parameters {
  matrix[p, K]                        beta; // kth column is the vector of coefficients for the kth dataset
  vector[p]                           beta_mean;
  vector<lower=0>[p]                  beta_sd;
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;
}
model {
  beta_mean       ~ normal(meta_mean_mean, meta_mean_sd); // hyperprior for mean of coefficients
  beta_sd         ~ normal(meta_sd_mean, meta_sd_sd); // hyperprior for sd of coefficients (half-normal)
  for ( j in 1:p ) {
    beta[j, ]     ~ normal(beta_mean[j], beta_sd[j]);
  }

  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
  else {
    dispersion ~ normal(disp_mean, disp_sd); // half-normal prior for dispersion
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], dispersion[k], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
}

