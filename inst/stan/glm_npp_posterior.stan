functions{
#include include/expfam_loglik.stan
#include include/finda0closest.stan
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
  real                                disp_mean; // mean for the half-normal prior on dispersion parameter
  real<lower=0>                       disp_sd; // sd for the half-normal prior on dispersion parameter
  int<lower=0>                        s; // number of a0s for which we have log nc
  vector<lower=0,upper=1>[s]          a0_lognc;
  matrix[s,K-1]                       lognc; // the j-th column is the log nc for a0_lognc using the j-th historical datasets
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
  vector<lower=0,upper=1>[K-1] a0_vals;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  if ( a0_shape1 != 1 || a0_shape2 != 1)
    a0_vals ~ beta(a0_shape1, a0_shape2);
  beta ~ normal(mean_beta, sd_beta); // prior for beta
  if ( dist <= 2 ) {
    target += glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, 1.0, X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]); // current data likelihood

    for ( k in 2:K ) {
      target += a0_vals[k-1] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // power prior
    }
  }
  else {
    dispersion ~ normal(disp_mean, disp_sd); // half-normal prior for dispersion
    target +=  glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, dispersion[1], X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]); // current data likelihood

    for ( k in 2:K ) {
      target += a0_vals[k-1] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, dispersion[1], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);  // power prior
    }
  }
  // Subtract log nc from power prior
  for ( k in 2:K ) {
      target += -pp_lognc(a0_vals[k-1], a0_lognc, lognc[, k-1]);
    }
}

