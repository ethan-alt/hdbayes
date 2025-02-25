// fit PWE model
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
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of regression coefficients
  vector[n1]                      y1;            // event times in current data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  array[n1] int<lower=1,upper=J>  intindx;       // index giving interval into which obs i failed / was censored for current data
  array[n1] int<lower=0,upper=1>  death_ind;     // event indicator (1 = event; 0 = censored) for current data
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  vector[p]                       beta_mean;     // mean for normal initial prior on coefficients
  vector<lower=0>[p]              beta_sd;       //sd for normal initial prior on coefficients
  vector[J]                       hazard_mean;   // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;     // scale parameter for half-normal prior on baseline hazards
  int<lower=0,upper=1>            get_loglik;    // whether to generate log-likelihood matrix
}
transformed data{
  real lognc_hazard = normal_lccdf(0 | hazard_mean, hazard_sd);
}
parameters {
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  vector[p]             beta;         // regression coefficients
}
model {
  vector[n1]  eta        = X1 * beta;
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;

  // compute lambda[j] * (s[j] * s[j-1])
  cumblhaz = cumulative_sum( lambda[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  // half-normal prior on baseline hazards
  target += normal_lpdf(lambda | hazard_mean, hazard_sd) - lognc_hazard;

  // prior on beta
  target += normal_lpdf(beta | beta_mean, beta_sd);

  // current data likelihood
  for ( i in 1:n1 )
    target += pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);
}
generated quantities{
  vector[(get_loglik == 1) ? n1 : 0] log_lik; //to generate log likelihood matrices
  vector[n1]  eta        = X1 * beta;
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz   = cumulative_sum( lambda[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  if ( get_loglik == 1 ) {
    for( i in 1:n1 ){
      log_lik[i] = pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);
    }
  }
}

