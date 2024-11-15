// fit AFT models (Log-Normal, Weibull, or Log-Logistic)
// sample from power prior (Not posterior)
functions{
  real aft_model_lp(vector y_obs, vector y_cen, vector eta_obs, vector eta_cen, real scale, int dist) {
    real loglik = 0;

    // Compute likelihood
    if ( dist == 1 ) { // log-normal
      loglik += normal_lpdf(y_obs | eta_obs, scale);    // uncensored data
      loglik += normal_lccdf(y_cen | eta_cen, scale);   // censored data
    }
    else if ( dist == 2 ) { // log-logistic
      loglik += logistic_lpdf(y_obs | eta_obs, scale);  // uncensored data
      loglik += logistic_lccdf(y_cen | eta_cen, scale); // censored data
    }
    else if ( dist == 3 ) { // weibull
      loglik += gumbel_lpdf(-y_obs | -eta_obs, scale);  // uncensored data
      loglik += gumbel_lcdf(-y_cen | -eta_cen, scale);  // censored data
    }
    else reject("Distribution not supported.");

    return(loglik);
  }
}
data {
  int<lower=1,upper=3>       dist;
  int<lower=0>               n0_obs;           // number of uncensored observations in historical data
  int<lower=0>               n0_cen;           // number of censored observations in historical data
  int<lower=0>               p;                // number of regression coefficients (including intercept)
  vector[n0_obs]             y0_obs;           // log of observed event time (uncensored) in historical data
  vector[n0_cen]             y0_cen;           // log of censored time in current data
  matrix[n0_obs,p]           X0_obs;           // design mtx for historical data (uncensored)
  matrix[n0_cen,p]           X0_cen;           // design mtx for historical data (censored)
  real<lower=0,upper=1>      a0;               // power prior parameter
  vector[p]                  beta_mean;        // mean for normal initial prior on coefficients
  vector<lower=0>[p]         beta_sd;          // sd for normal initial prior on coefficients
  real                       scale_mean;       // location parameter for half-normal prior on scale
  real<lower=0>              scale_sd;         // scale parameter for half-normal prior on scale
}
transformed data {
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);
}
parameters {
  vector[p] beta;
  real<lower=0> scale;
}
model {
  // Temporary calculations
  vector[n0_obs] eta0_obs     = X0_obs * beta;     // linear predictor for historical data uncensored
  vector[n0_cen] eta0_cen     = X0_cen * beta;     // linear predictor for historical data censored

  // half-normal prior on scale
  target += normal_lpdf( scale | scale_mean, scale_sd ) - scale_prior_lognc;

  // normal initial prior on beta
  target += normal_lpdf( beta | beta_mean, beta_sd );

  // power prior
  target += a0 * aft_model_lp(y0_obs, y0_cen, eta0_obs, eta0_cen, scale, dist);
}

