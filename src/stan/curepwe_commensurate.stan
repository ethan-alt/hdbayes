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
  int<lower=0>                    n1;             // current data sample size
  int<lower=0>                    n0;             // historical data sample size
  int<lower=0>                    J;              // number of time intervals
  int<lower=0>                    p;              // number of regression coefficients
  vector[n1]                      y1;             // event times in current data
  vector[n0]                      y0;             // event times in historical data
  matrix[n1,p]                    X1;             // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;             // historical data design matrix (EXCLUDING intercept term)
  array[n1] int<lower=1,upper=J>  intindx;        // index giving interval into which obs i failed / was censored for current data
  array[n0] int<lower=1,upper=J>  intindx0;       // index giving interval into which obs i failed / was censored for historical data
  array[n1] int<lower=0,upper=1>  death_ind;      // event indicator (1 = event; 0 = censored) for current data
  array[n0] int<lower=0,upper=1>  death_ind0;     // event indicator (1 = event; 0 = censored) for historical data
  vector[J+1]                     breaks;         // J+1-dim interval of breaks
  vector[p]                       beta0_mean;     // mean for normal initial prior on coefficients in historical data model
  vector<lower=0>[p]              beta0_sd;       // sd for normal initial prior on coefficients in historical data model
  real<lower=0,upper=1>           p_spike;        // prior probability of spike
  real                            mu_spike;       // location parameter for half-normal prior (spike)
  real<lower=0>                   sigma_spike;    // scale parameter for half-normal prior (spike)
  real                            mu_slab;        // location parameter for half-normal prior (slab)
  real<lower=0>                   sigma_slab;     // scale parameter for half-normal prior (slab)
  vector[J]                       hazard_mean;    // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;      // scale parameter for half-normal prior on baseline hazards
  real                            logit_p_cured_mean; // mean for normal prior on logit_p_cured
  real<lower=0>                   logit_p_cured_sd;   // sd for normal prior on logit_p_cured
  int<lower=0,upper=1>            get_loglik;     // whether to generate log-likelihood matrix
}
transformed data {
  vector[n1] logcens;
  vector[n0] logcens0;
  real lognc_spike   = normal_lccdf(0 | mu_spike, sigma_spike); // \Phi(mu_spike / sigma_spike)
  real lognc_slab    = normal_lccdf(0 | mu_slab, sigma_slab);   // \Phi(mu_slab / sigma_slab)
  real lognc_hazard  = normal_lccdf(0 | hazard_mean, hazard_sd);

  // compute censoring indicators
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);

  for ( i in 1:n0 )
    logcens0[i] = log(1 - death_ind0[i]);
}
parameters {
  real                  logit_p_cured;  // logit of proportion of cured individuals for current data model
  real                  logit_p_cured0; // logit of proportion of cured individuals for historical data model
  vector<lower=0>[J]    lambda;       // the J hazards for each interval for current data model
  vector<lower=0>[J]    lambda0;      // the J hazards for each interval for historical data model
  vector[p]             beta;
  vector[p]             beta0;
  vector<lower=0>[p]    comm_prec;    // commensurability parameter
}
transformed parameters{
  real<lower=0,upper=1> p_cured  = inv_logit(logit_p_cured); // proportion of cured individuals for current data model
  real<lower=0,upper=1> p_cured0 = inv_logit(logit_p_cured0); // proportion of cured individuals for historical data model
  vector[p] comm_sd = inv_sqrt(comm_prec);   // sd for commensurate prior
}
model {
  vector[n1]  eta         = X1 * beta;
  vector[n0]  eta0        = X0 * beta0;
  vector[J]   log_lambda  = log(lambda);
  vector[J]   log_lambda0 = log(lambda0);
  vector[J-1] cumblhaz;
  vector[J-1] cumblhaz0;

  // compute lambda[j] * (s[j] * s[j-1])
  cumblhaz = cumulative_sum( lambda[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );
  cumblhaz0 = cumulative_sum( lambda0[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  // half-normal prior on baseline hazards
  target += normal_lpdf(lambda | hazard_mean, hazard_sd) - lognc_hazard;
  target += normal_lpdf(lambda0 | hazard_mean, hazard_sd) - lognc_hazard;

  // prior on logit_p_cured
  target += normal_lpdf(logit_p_cured | logit_p_cured_mean, logit_p_cured_sd);
  target += normal_lpdf(logit_p_cured0 | logit_p_cured_mean, logit_p_cured_sd);

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
  for ( i in 1:n0 ){
    target += log_mix(
      p_cured0,
      logcens0[i],
      pwe_lpdf(y0[i] | eta0[i], lambda0, log_lambda0, breaks, intindx0[i], death_ind0[i], cumblhaz0)
    );
  }

  // current data likelihood
  for ( i in 1:n1 ){
    target += log_mix(
      p_cured,
      logcens[i],
      pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz)
    );
  }
}
generated quantities{
  vector[(get_loglik == 1) ? n1 : 0] log_lik; //to generate log likelihood matrices
  vector[n1]  eta        = X1 * beta;
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz   = cumulative_sum( lambda[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  if ( get_loglik == 1 ) {
    for( i in 1:n1 ){
      log_lik[i] = log_mix(
        p_cured,
        logcens[i],
        pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz)
      );
    }
  }
}

