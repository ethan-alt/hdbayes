// fit AFT models (Log-Normal, Weibull, or Log-Logistic)
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
  int<lower=1,upper=3>       dist;
  int<lower=0>               n0_obs;           // number of uncensored observations in historical data
  int<lower=0>               n0_cen;           // number of censored observations in historical data
  int<lower=0>               p;                // number of regression coefficients (including intercept)
  vector[n0_obs]             y0_obs;           // log of observed event time (uncensored) in historical data
  vector[n0_cen]             y0_cen;           // log of censored time in current data
  matrix[n0_obs,p]           X0_obs;           // design mtx for historical data (uncensored)
  matrix[n0_cen,p]           X0_cen;           // design mtx for historical data (censored)
  vector[p]                  beta0_mean;       // mean for normal initial prior on coefficients in historical data model
  vector<lower=0>[p]         beta0_sd;         // sd for normal initial prior on coefficients in historical data model
  real<lower=0,upper=1>      p_spike;          // prior probability of spike
  real                       mu_spike;         // location parameter for half-normal prior (spike)
  real<lower=0>              sigma_spike;      // scale parameter for half-normal prior (spike)
  real                       mu_slab;          // location parameter for half-normal prior (slab)
  real<lower=0>              sigma_slab;       // scale parameter for half-normal prior (slab)
  real                       scale_mean;       // location parameter for half-normal prior on scale
  real<lower=0>              scale_sd;         // scale parameter for half-normal prior on scale
}
transformed data {
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);
  real lognc_spike       = normal_lccdf(0 | mu_spike, sigma_spike); // \Phi(mu_spike / sigma_spike)
  real lognc_slab        = normal_lccdf(0 | mu_slab, sigma_slab);   // \Phi(mu_slab / sigma_slab)
}
parameters {
  vector[p]          beta;  // regression coefficients for current data model
  vector[p]          beta0; // regression coefficients for historical data model
  real<lower=0>      scale0;
  vector<lower=0>[p] comm_prec;  // commensurability parameter
}
transformed parameters {
  vector[p] comm_sd = inv_sqrt(comm_prec);   // sd for commensurate prior
}
model {
  // Temporary calculations
  vector[n0_obs] eta0_obs     = X0_obs * beta0;    // linear predictor for historical data uncensored
  vector[n0_cen] eta0_cen     = X0_cen * beta0;    // linear predictor for historical data censored

  // half-normal prior on scale
  target += normal_lpdf( scale0 | scale_mean, scale_sd ) - scale_prior_lognc;

  // spike and slab prior for commensurability:
  for ( i in 1:p ) {
    target += log_mix(
      p_spike
      , normal_lpdf(comm_prec[i] | mu_spike, sigma_spike) - lognc_spike
      , normal_lpdf(comm_prec[i] | mu_slab, sigma_slab) - lognc_slab
    );
  }
  // initial prior for beta0
  target += normal_lpdf(beta0 | beta0_mean, beta0_sd);
  // commensurate prior for beta
  target += normal_lpdf(beta  | beta0, comm_sd);

  // historical data likelihood
  target += aft_model_lp(y0_obs, y0_cen, eta0_obs, eta0_cen, scale0, dist);
}

