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
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of regression coefficients
  vector[n0]                      y0;            // event times in historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  array[n0] int<lower=1,upper=J>  intindx0;      // index giving interval into which obs i failed / was censored for historical data
  array[n0] int<lower=0,upper=1>  death_ind0;    // event indicator (1 = event; 0 = censored) for historical data
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  int<lower = 0>                  K;             // number of strata
  array[n0] int<lower=1,upper=K>  stratumID0;    // strata assignment for historical data
  vector<lower=0, upper=1>[K]     a0s;           // power prior parameter for each stratum
  vector[p]                       beta_mean;     // mean for normal initial prior on coefficients
  vector<lower=0>[p]              beta_sd;       // sd for normal initial prior on coefficients
  vector[J]                       hazard_mean;   // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;     // scale parameter for half-normal prior on baseline hazards
  real                            logit_p_cured_mean; // mean for normal prior on logit_p_cured
  real<lower=0>                   logit_p_cured_sd;   // sd for normal prior on logit_p_cured
}
transformed data{
  vector[n0] logcens0;
  real lognc_hazard = normal_lccdf(0 | hazard_mean, hazard_sd);

  // compute censoring indicator
  for ( i in 1:n0 )
    logcens0[i] = log(1 - death_ind0[i]);
}
parameters {
  vector[K]             logit_p_curedVec; // Kx1 vector of logit of proportions of cured individuals
  matrix<lower=0>[J,K]  lambdaMat;        // JxK matrix of hazards for each interval
  matrix[p,K]           betaMat;          // pxK matrix of regression coefficients; p = number of covars, K = number of strata
}
transformed parameters{
  vector<lower=0,upper=1>[K] p_curedVec = inv_logit(logit_p_curedVec); // Kx1 vector of proportions of cured individuals
}
model {
  matrix[n0,K]  eta0Mat       = X0 * betaMat;
  matrix[J,K]   log_lambdaMat = log(lambdaMat);
  matrix[J-1,K] cumblhazMat;

  // noninformative initial prior
  for ( k in 1:K ) {
    target += normal_lpdf(betaMat[, k] | beta_mean, beta_sd);
    // half-normal prior on baseline hazards
    target += normal_lpdf(lambdaMat[, k] | hazard_mean, hazard_sd) - lognc_hazard;

    // compute lambda[j] * (s[j] * s[j-1])
    cumblhazMat[, k] = cumulative_sum( lambdaMat[1:(J-1), k] .* ( breaks[2:J] - breaks[1:(J-1)] ) );
  }

  // prior on logit_p_curedVec
  target += normal_lpdf(logit_p_curedVec | logit_p_cured_mean, logit_p_cured_sd);

  // power prior
  for ( i in 1:n0 ){
    target += a0s[ stratumID0[i] ] * log_mix(
      p_curedVec[ stratumID0[i] ],
      logcens0[i],
      pwe_lpdf(y0[i] | eta0Mat[i, stratumID0[i]], lambdaMat[, stratumID0[i]], log_lambdaMat[, stratumID0[i]], breaks, intindx0[i], death_ind0[i], cumblhazMat[, stratumID0[i]])
    );
  }
}

