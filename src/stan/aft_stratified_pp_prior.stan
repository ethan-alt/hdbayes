// fit AFT models (Log-Normal, Weibull, or Log-Logistic)
functions{
  // compute log-likelihood for an uncensored observation
  real aft_model_obs_lpdf(real y_obs, real eta_obs, real scale, int dist) {
    real loglik;

    // Compute likelihood
    if ( dist == 1 ) { // log-normal
      loglik = normal_lpdf(y_obs | eta_obs, scale);    // uncensored data
    }
    else if ( dist == 2 ) { // log-logistic
      loglik = logistic_lpdf(y_obs | eta_obs, scale);  // uncensored data
    }
    else if ( dist == 3 ) { // weibull
      loglik = gumbel_lpdf(-y_obs | -eta_obs, scale);  // uncensored data
    }
    else reject("Distribution not supported.");

    return(loglik);
  }

  // compute log-likelihood for a censored observation
  real aft_model_cen_lpdf(real y_cen, real eta_cen, real scale, int dist) {
    real loglik;

    // Compute likelihood
    if ( dist == 1 ) { // log-normal
      loglik = normal_lccdf(y_cen | eta_cen, scale);   // censored data
    }
    else if ( dist == 2 ) { // log-logistic
      loglik = logistic_lccdf(y_cen | eta_cen, scale); // censored data
    }
    else if ( dist == 3 ) { // weibull
      loglik = gumbel_lcdf(-y_cen | -eta_cen, scale);  // censored data
    }
    else reject("Distribution not supported.");

    return(loglik);
  }
}
data {
  int<lower=1,upper=3>                 dist;
  int<lower=0>                         n0_obs;          // number of uncensored observations in historical data
  int<lower=0>                         n0_cen;          // number of censored observations in historical data
  int<lower=0>                         p;               // number of regression coefficients (including intercept)
  vector[n0_obs]                       y0_obs;          // log of observed event time (uncensored) in historical data
  vector[n0_cen]                       y0_cen;          // log of censored time in current data
  matrix[n0_obs,p]                     X0_obs;          // design mtx for historical data (uncensored)
  matrix[n0_cen,p]                     X0_cen;          // design mtx for historical data (censored)
  int<lower = 0>                       K;               // number of strata
  array[n0_obs] int<lower=1,upper=K>   stratumID0_obs;  // strata assignment for historical data (uncensored)
  array[n0_cen] int<lower=1,upper=K>   stratumID0_cen;  // strata assignment for historical data (censored)
  vector<lower=0, upper=1>[K]          a0s;             // power prior parameter for each stratum
  vector[p]                            beta_mean;       // mean for normal initial prior on coefficients
  vector<lower=0>[p]                   beta_sd;         // sd for normal initial prior on coefficients
  real                                 scale_mean;      // location parameter for half-normal prior on scale
  real<lower=0>                        scale_sd;        // scale parameter for half-normal prior on scale
}
transformed data {
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);
}
parameters {
  matrix[p, K]           betaMat; // each column is a stratum-specific vector of coefficients
  row_vector<lower=0>[K] scaleVec;
}
model {
  // Temporary calculations
  matrix[n0_obs, K] Eta0_obs     = X0_obs * betaMat;     // linear predictors for historical data uncensored
  matrix[n0_cen, K] Eta0_cen     = X0_cen * betaMat;     // linear predictors for historical data censored

  for ( k in 1:K ) {
    // half-normal prior on scale
    target += normal_lpdf( scaleVec[k] | scale_mean, scale_sd ) - scale_prior_lognc;
    // normal initial prior on beta
    target += normal_lpdf( betaMat[, k] | beta_mean, beta_sd );
  }

  // power prior
  for( i in 1:n0_obs ){
    target += a0s[ stratumID0_obs[i] ] * aft_model_obs_lpdf(y0_obs[i] | Eta0_obs[i, stratumID0_obs[i]], scaleVec[ stratumID0_obs[i] ], dist);
  }
  for( i in 1:n0_cen ){
    target += a0s[ stratumID0_cen[i] ] * aft_model_cen_lpdf(y0_cen[i] | Eta0_cen[i, stratumID0_cen[i]], scaleVec[ stratumID0_cen[i] ], dist);
  }
}

