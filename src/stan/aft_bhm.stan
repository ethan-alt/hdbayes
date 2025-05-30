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
  int<lower=0>               n;                // total number of observations in current data
  int<lower=0>               n_obs;            // number of uncensored observations in current data
  int<lower=0>               n_cen;            // number of censored observations in current data
  int<lower=0>               n0_obs;           // number of uncensored observations in historical data
  int<lower=0>               n0_cen;           // number of censored observations in historical data
  int<lower=0>               p;                // number of regression coefficients (including intercept)
  vector[n_obs]              y_obs;            // log of observed event time (uncensored) in current data
  vector[n_cen]              y_cen;            // log of censored time in current data
  matrix[n_obs,p]            X_obs;            // design mtx for current data (uncensored)
  matrix[n_cen,p]            X_cen;            // design mtx for current data (censored)
  vector[n0_obs]             y0_obs;           // log of observed event time (uncensored) in historical data
  vector[n0_cen]             y0_cen;           // log of censored time in current data
  matrix[n0_obs,p]           X0_obs;           // design mtx for historical data (uncensored)
  matrix[n0_cen,p]           X0_cen;           // design mtx for historical data (censored)
  vector[p]                  meta_mean_mean;   // mean for the hyperprior for mean of each coefficient
  vector<lower=0>[p]         meta_mean_sd;     // sd for the hyperprior for mean of each coefficient
  vector[p]                  meta_sd_mean;     // location parameter for half-normal prior on sd of each coefficient
  vector<lower=0>[p]         meta_sd_sd;       // scale parameter for half-normal prior on sd of each coefficient
  real                       scale_mean;       // location parameter for half-normal prior on scale
  real<lower=0>              scale_sd;         // scale parameter for half-normal prior on scale
  int<lower=0,upper=1>       get_loglik;       // whether to generate log-likelihood matrix
}
transformed data {
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);
  real lognc_beta_sd     = normal_lccdf(0 | meta_sd_mean, meta_sd_sd);
}
parameters {
  vector[p]             beta_raw;
  vector[p]             beta0_raw;
  vector[p]             beta_mean;
  vector<lower=0>[p]    beta_sd;
  real<lower=0>         scale;
  real<lower=0>         scale0;
}
transformed parameters{
  vector[p] beta              = beta_mean + beta_sd .* beta_raw;   // coefficients for current data model
  vector[p] beta0             = beta_mean + beta_sd .* beta0_raw;  // coefficients for historical data model
}
model {
  // Temporary calculations
  vector[n_obs]  eta_obs      = X_obs * beta;      // linear predictor for current data uncensored
  vector[n_cen]  eta_cen      = X_cen * beta;      // linear predictor for current data censored
  vector[n0_obs] eta0_obs     = X0_obs * beta0;    // linear predictor for historical data uncensored
  vector[n0_cen] eta0_cen     = X0_cen * beta0;    // linear predictor for historical data censored

  // half-normal prior on scale
  target += normal_lpdf( scale | scale_mean, scale_sd ) - scale_prior_lognc;
  target += normal_lpdf( scale0 | scale_mean, scale_sd ) - scale_prior_lognc;

  // hyperprior for mean of coefficients
  target += normal_lpdf(beta_mean | meta_mean_mean, meta_mean_sd);
  // hyperprior for sd of coefficients (half-normal)
  target += normal_lpdf(beta_sd | meta_sd_mean, meta_sd_sd) - lognc_beta_sd;

  // prior on coefficients
  target += std_normal_lpdf(beta_raw);  // implies beta ~ normal(beta_mean, beta_sd);
  target += std_normal_lpdf(beta0_raw); // implies beta0 ~ normal(beta_mean, beta_sd);

  // historical data likelihood
  target += aft_model_lp(y0_obs, y0_cen, eta0_obs, eta0_cen, scale0, dist);

  // current data likelihood
  target += aft_model_lp(y_obs, y_cen, eta_obs, eta_cen, scale, dist);
}
generated quantities{
  vector[(get_loglik == 1) ? n : 0] log_lik;                  // log likelihood matrix
  vector[n_obs]                     eta_obs   = X_obs * beta; // linear predictor for current data uncensored
  vector[n_cen]                     eta_cen   = X_cen * beta; // linear predictor for current data censored

  if ( get_loglik == 1 ) {
    for( i in 1:n_obs ){
      log_lik[i] = aft_model_obs_lpdf(y_obs[i] | eta_obs[i], scale, dist);
    }
    for( i in (n_obs+1):n ){
      log_lik[i] = aft_model_cen_lpdf(y_cen[i-n_obs] | eta_cen[i-n_obs], scale, dist);
    }
  }
}

