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
  vector[p]                           beta0_mean;
  vector<lower=0>[p]                  beta0_sd;
  vector[K]                           disp_mean; // mean for the half-normal prior for dispersion of each dataset
  vector<lower=0>[K]                  disp_sd; // sd for the half-normal prior for dispersion of each dataset
  vector<lower=0>[p]                  tau; // commensurate prior parameter
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
transformed data{
  vector[p] comm_sd = inv_sqrt(tau);   // sd hyperparameter for commensurate prior
}
parameters {
  vector[p] beta;  // regression coefficients for current data
  vector[p] beta0; // regression coefficients for historical datasets
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;
}
model {
  beta0   ~ normal(beta0_mean, beta0_sd); // initial prior for beta0
  beta    ~ normal(beta0, comm_sd);       // commensurate prior for beta
  if ( dist <= 2 ) {
    target += glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, 1.0, X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]);    // current data likelihood

    for ( k in 2:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta0, 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // historical data likelihood
    }
  }
  else {
    dispersion ~ normal(disp_mean, disp_sd); // half-normal prior for dispersion
    target += glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, dispersion[1], X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]);    // current data likelihood

    for ( k in 2:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta0, dispersion[k], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // historical data likelihood
    }
  }
}

