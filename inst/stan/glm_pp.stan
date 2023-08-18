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
  vector[p]                           mean_beta; // mean for normal initial prior on coefficients
  vector<lower=0>[p]                  sd_beta; //sd for normal initial prior on coefficients
  vector<lower=0,upper=1>[K]          a0_vals; // power prior parameter (fixed) for each dataset (with first element = 1 for current data)
  real                                disp_mean; // mean for the half-normal prior on dispersion parameter
  real<lower=0>                       disp_sd; // sd for the half-normal prior on dispersion parameter
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
parameters {
  vector[p]     beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
}
model {
  beta ~ normal(mean_beta, sd_beta); // initial prior beta ~ N(mu0, sd0)
  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      target += a0_vals[k] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // current data likelihood and power prior
    }
  }
  else {
    dispersion ~ normal(disp_mean, disp_sd); // half-normal prior for dispersion
    for ( k in 1:K ) {
      target += a0_vals[k] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, dispersion[1], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);  // current data likelihood and power prior
    }
  }
}

