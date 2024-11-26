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

  // find index of x, j, such that x0 is closest to x[j] without
  // going over. Uses binary search algorithm
  int findClosestIndex(real x0, vector x) {
    int K = rows(x);
    int i = 1;
    int j = K;
    int mid;
    // check corner cases
    if ( x0 < x[2] )
      return 1;
    if ( x0 == x[K] ) {
      return K;
    }
    // conduct binary search
    while ( i <= j ) {
      mid = (i + j) %/% 2;
      // if x0 < x[mid], index must lie in left half
      if ( x0 < x[mid] ) {
        // if x0 is larger than x[mid-1], return mid-1; else update j
        if ( mid > 2 &&  x0 > x[mid - 1] )
          return mid - 1;
        j = mid;
      }
      // otherwise, index must lie in right half
      else {
        // if x0 is less than x[mid + 1], return mid; else update i
        if ( mid < K && x0 < x[mid + 1] )
          return mid;
        i = mid + 1;
      }
    }
    reject("Error in finding midpoint");
    return(0); // never reached
  }

  // approximate lognc of power prior
  //
  // * @param a0       power prior param to obtain lognc
  // * @param a0vec    fine grid of power prior parameters for which we have estimates
  // * @param lognca0  estimate lognc pertaining to fine grid a0vec
  //
  // * @return linearly interpolated log normalizing constant.
  real pp_lognc(real a0, vector a0vec, vector lognca0) {
    // find index of a0vec closest to a0
    int i = findClosestIndex(a0, a0vec);
    // if not exact match, use linear interpolation to get estimated lognc
    if ( a0 != a0vec[i] ) {
      real x1 = a0vec[i];
      real x2 = a0vec[i + 1];
      real y1 = lognca0[i];
      real y2 = lognca0[i + 1];
      return y1 + (y2 - y1) * (a0 - x1) / (x2 - x1);
    }
    return lognca0[i];
  }

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return
      -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x);
  }
}
data {
  int<lower=1,upper=3>                dist;
  int<lower=0>                        n;                // total number of observations in current data
  int<lower=0>                        n_obs;            // number of uncensored observations in current data
  int<lower=0>                        n_cen;            // number of censored observations in current data
  int<lower=0>                        n0_obs;           // number of uncensored observations in historical data
  int<lower=0>                        n0_cen;           // number of censored observations in historical data
  int<lower=0>                        p;                // number of regression coefficients (including intercept)
  vector[n_obs]                       y_obs;            // log of observed event time (uncensored) in current data
  vector[n_cen]                       y_cen;            // log of censored time in current data
  matrix[n_obs,p]                     X_obs;            // design mtx for current data (uncensored)
  matrix[n_cen,p]                     X_cen;            // design mtx for current data (censored)
  vector[n0_obs]                      y0_obs;           // log of observed event time (uncensored) in historical data
  vector[n0_cen]                      y0_cen;           // log of censored time in current data
  matrix[n0_obs,p]                    X0_obs;           // design mtx for historical data (uncensored)
  matrix[n0_cen,p]                    X0_cen;           // design mtx for historical data (censored)
  vector[p]                           beta_mean;        // mean for normal initial prior on coefficients
  vector<lower=0>[p]                  beta_sd;          // sd for normal initial prior on coefficients
  real                                scale_mean;       // location parameter for half-normal prior on scale
  real<lower=0>                       scale_sd;         // scale parameter for half-normal prior on scale
  int<lower=0>                        s;                // number of a0s for which we have log nc
  vector<lower=0,upper=1>[s]          a0_lognc;
  vector[s]                           lognc;            // log normalizing constant for each a0 value specified in a0_lognc using the historical data
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  real<lower=0,upper=1>               a0_lower;         // lower bounds for a0s
  real<lower=a0_lower,upper=1>        a0_upper;         // upper bounds for a0s
  int<lower=0,upper=1>                get_loglik;       // whether to generate log-likelihood matrix
}
transformed data{
  real scale_prior_lognc = normal_lccdf(0 | scale_mean, scale_sd);
  real logit_a0_lower    = logit(a0_lower);
  real logit_a0_upper    = logit(a0_upper);
  real lognc_logit_a0    = 0;

  if( a0_upper != 1 || a0_lower != 0 )
    lognc_logit_a0 = lognc_logit_a0 + log_diff_exp( beta_lcdf(a0_upper | a0_shape1, a0_shape2), beta_lcdf(a0_lower | a0_shape1, a0_shape2) );
}
parameters {
  vector[p] beta;
  real<lower=0> scale;
  real<lower=logit_a0_lower,upper=logit_a0_upper> logit_a0;
}
transformed parameters {
  real<lower=a0_lower,upper=a0_upper> a0 = inv_logit(logit_a0);
}
model {
  // Temporary calculations
  vector[n_obs]  eta_obs      = X_obs * beta;      // linear predictor for current data uncensored
  vector[n_cen]  eta_cen      = X_cen * beta;      // linear predictor for current data censored
  vector[n0_obs] eta0_obs     = X0_obs * beta;     // linear predictor for historical data uncensored
  vector[n0_cen] eta0_cen     = X0_cen * beta;     // linear predictor for historical data censored

  // half-normal prior on scale
  target += normal_lpdf( scale | scale_mean, scale_sd ) - scale_prior_lognc;

  // normal initial prior on beta
  target += normal_lpdf( beta | beta_mean, beta_sd );

  // prior on logit(a0)
  target += logit_beta_lpdf(logit_a0 | a0_shape1, a0_shape2) - lognc_logit_a0;

  // power prior
  target += a0 * aft_model_lp(y0_obs, y0_cen, eta0_obs, eta0_cen, scale, dist);

  // current data likelihood
  target += aft_model_lp(y_obs, y_cen, eta_obs, eta_cen, scale, dist);

  // Subtract log nc from power prior
  target += -pp_lognc(a0, a0_lognc, lognc);
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

