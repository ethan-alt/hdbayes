// fit AFT models (Log-Normal, Weibull, or Log-Logistic)
// sample from the LEAP (Not posterior)
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

  real aft_model_obs_mixture_lp(real y_obs, row_vector eta_obs, row_vector scale, vector log_probs, int dist) {
    int       K = rows(log_probs);
    vector[K] contribs;              // log likelihood contribution of each component

    for (k in 1:K) {
      contribs[k] = log_probs[k] + aft_model_obs_lpdf(y_obs | eta_obs[k], scale[k], dist);
    }

    return(log_sum_exp(contribs));
  }

  real aft_model_cen_mixture_lp(real y_cen, row_vector eta_cen, row_vector scale, vector log_probs, int dist) {
    int       K = rows(log_probs);
    vector[K] contribs;              // log likelihood contribution of each component

    for (k in 1:K) {
      contribs[k] = log_probs[k] + aft_model_cen_lpdf(y_cen | eta_cen[k], scale[k], dist);
    }

    return(log_sum_exp(contribs));
  }

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return
      -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x);
  }
}
data {
  int<lower=1,upper=3>       dist;
  int<lower=0>               K;                // number of components in mixture
  int<lower=0>               n0_obs;           // number of uncensored observations in historical data
  int<lower=0>               n0_cen;           // number of censored observations in historical data
  int<lower=0>               p;                // number of regression coefficients (including intercept)
  vector[n0_obs]             y0_obs;           // log of observed event time (uncensored) in historical data
  vector[n0_cen]             y0_cen;           // log of censored time in current data
  matrix[n0_obs,p]           X0_obs;           // design mtx for historical data (uncensored)
  matrix[n0_cen,p]           X0_cen;           // design mtx for historical data (censored)
  vector[p]                  beta_mean;        // mean for normal initial prior on coefficients
  vector<lower=0>[p]         beta_sd;          // sd for normal initial prior on coefficients
  real                       scale_mean;       // location parameter for half-normal prior on scale
  real<lower=0>              scale_sd;         // scale parameter for half-normal prior on scale
  vector<lower=0>[K]         prob_conc;        // dirichlet parameters for probs
  real<lower=0,upper=1> gamma_lower;           // lower bound for probability of being exchangeable
  real<lower=gamma_lower,upper=1> gamma_upper; // upper bound for probability of being exchangeable
}
transformed data {
  real gamma_shape1      = prob_conc[1];
  real gamma_shape2      = sum(prob_conc[2:K]);
  real logit_gamma_lower = logit(gamma_lower);
  real logit_gamma_upper = logit(gamma_upper);
  real lognc_logit_gamma = 0;
  int  K_gt_2            = (K > 2) ? (1) : (0);
  vector[K-1] conc_delta = prob_conc[2:K];
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);

  if( gamma_upper != 1 || gamma_lower != 0 )
    lognc_logit_gamma = lognc_logit_gamma + log_diff_exp( beta_lcdf(gamma_upper | gamma_shape1, gamma_shape2), beta_lcdf(gamma_lower | gamma_shape1, gamma_shape2) );
}
parameters {
  matrix[p, K] betaMat;
  row_vector<lower=0>[K] scaleVec;
  simplex[K-1]           delta_raw;
  real<lower=logit_gamma_lower,upper=logit_gamma_upper> logit_gamma; // logit of probability of being exchangeable
}
model {
  // Temporary calculations
  matrix[n0_obs, K] Eta0_obs = (X0_obs * betaMat);        // linear predictors for historical data uncensored
  matrix[n0_cen, K] Eta0_cen = (X0_cen * betaMat);        // linear predictors for historical data censored
  vector[K]         logprobs = append_row(logit_gamma, log(delta_raw)) - log1p_exp(logit_gamma);

  // Priors
  // If two components, get beta prior on gamma;
  // If >2 components, get a dirichlet prior on raw delta
  target += logit_beta_lpdf(logit_gamma | gamma_shape1, gamma_shape2) - lognc_logit_gamma;

  if (K_gt_2)
    target += dirichlet_lpdf(delta_raw | conc_delta);

  for ( k in 1:K ) {
    // half-normal prior on scale
    target += normal_lpdf( scaleVec[k] | scale_mean, scale_sd ) - scale_prior_lognc;
    // normal initial prior on beta
    target += normal_lpdf( betaMat[, k] | beta_mean, beta_sd );
  }

  // historical data likelihood
  for ( i in 1:n0_obs ){
    target += aft_model_obs_mixture_lp(y0_obs[i], Eta0_obs[i, ], scaleVec, logprobs, dist);
  }
  for ( i in 1:n0_cen ){
    target += aft_model_cen_mixture_lp(y0_cen[i], Eta0_cen[i, ], scaleVec, logprobs, dist);
  }
}

