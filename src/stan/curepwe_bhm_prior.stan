// fit standard cure rate model: S(t) = pcure + (1 - pcure) * survival function for PWE model
functions {
  //' PWE log likelihood
  //' @param y               failure / censoring time
  //' @param eta             linear predictor
  //' @param lambda          vector of baseline hazards
  //' @param log_lambda      vector of log(lambda)
  //' @param breaks          (J+1)-dim vector giving intervals
  //' @param j               index of the interval for which the individual failed / was censored: j \in {1, ..., J}
  //' @param death_ind       integer giving whether individual died (death_ind == 1) or was censored (death_ind == 0)
  //' @param lambda_d_breaks (J-1)-dim vector giving lambda[j] * (s[j] - s[j-1]), j = 1, ..., J-1
  real pwe_lpdf(real y, real eta, vector lambda, vector log_lambda, vector breaks, int j, int death_ind, vector cumblhaz) {
    real loglik;

    // Initialize cumhaz to lambda[j] ( y - s[j] )
    real cumhaz = lambda[j] * ( y - breaks[j] );

    // add on (sum_{g=1}^{j-1} lambda[j] * ( s[j] - s[j-1] )
    if ( j > 1 )
      cumhaz += cumblhaz[j-1];

    // Initialize loglik = - cumhaz * exp(x'beta) = log(Survival)
    loglik = -cumhaz * exp(eta);

    // If individual experienced event, add log hazard: log(lambda[j]) + x'beta
    if ( death_ind == 1 )
      loglik += log_lambda[j] + eta;

    return(loglik);
  }
}
data {
  int<lower=0>                    n0;             // historical data sample size
  int<lower=0>                    J;              // number of time intervals
  int<lower=0>                    p;              // number of regression coefficients
  vector[n0]                      y0;             // event times in historical data
  matrix[n0,p]                    X0;             // historical data design matrix (EXCLUDING intercept term)
  array[n0] int<lower=1,upper=J>  intindx0;       // index giving interval into which obs i failed / was censored for historical data
  array[n0] int<lower=0,upper=1>  death_ind0;     // event indicator (1 = event; 0 = censored) for historical data
  vector[J+1]                     breaks;         // J+1-dim interval of breaks
  vector[p]                       meta_mean_mean; // mean for the hyperprior for mean of each coefficient
  vector<lower=0>[p]              meta_mean_sd;   // sd for the hyperprior for mean of each coefficient
  vector[p]                       meta_sd_mean;   // mean for the hyperprior for sd of each coefficient
  vector<lower=0>[p]              meta_sd_sd;     // sd for the hyperprior for sd of each coefficient
  vector[J]                       hazard_mean;    // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;      // scale parameter for half-normal prior on baseline hazards
  real                            logit_p_cured_mean; // mean for normal prior on logit_p_cured
  real<lower=0>                   logit_p_cured_sd;   // sd for normal prior on logit_p_cured
}
transformed data {
  vector[n0] logcens0;
  real lognc_beta_sd = normal_lccdf(0 | meta_sd_mean, meta_sd_sd);
  real lognc_hazard  = normal_lccdf(0 | hazard_mean, hazard_sd);

  // compute censoring indicator
  for ( i in 1:n0 )
    logcens0[i] = log(1 - death_ind0[i]);
}
parameters {
  real                  logit_p_cured0; // logit of proportion of cured individuals for historical data model
  vector<lower=0>[J]    lambda0;        // the J hazards for each interval for historical data model
  vector[p]             beta0_raw;
  vector[p]             beta_mean;
  vector<lower=0>[p]    beta_sd;
}
transformed parameters{
  real<lower=0,upper=1> p_cured0 = inv_logit(logit_p_cured0); // proportion of cured individuals for historical data model
  vector[p] beta0 = beta_mean + beta_sd .* beta0_raw;  // regression coefficients for historical data model
}
model {
  vector[n0]  eta0       = X0 * beta0;
  vector[J]   log_lambda0 = log(lambda0);
  vector[J-1] cumblhaz0;

  // compute lambda[j] * (s[j] * s[j-1])
  cumblhaz0 = cumulative_sum( lambda0[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  // half-normal prior on baseline hazards
  target += normal_lpdf(lambda0 | hazard_mean, hazard_sd) - lognc_hazard;

  // prior on logit_p_cured
  target += normal_lpdf(logit_p_cured0 | logit_p_cured_mean, logit_p_cured_sd);

  target += normal_lpdf(beta_mean | meta_mean_mean, meta_mean_sd); // hyperprior for mean of coefficients
  target += normal_lpdf(beta_sd | meta_sd_mean, meta_sd_sd) - lognc_beta_sd; // hyperprior for sd of coefficients (half-normal)

  // prior on regression coefficients
  target += std_normal_lpdf(beta0_raw); // implies beta0 ~ normal(beta_mean, beta_sd);

  // historical data likelihood
  for ( i in 1:n0 ){
    target += log_mix(
      p_cured0,
      logcens0[i],
      pwe_lpdf(y0[i] | eta0[i], lambda0, log_lambda0, breaks, intindx0[i], death_ind0[i], cumblhaz0)
    );
  }
}

