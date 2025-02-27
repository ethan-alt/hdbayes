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

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return(
      -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x)
    );
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
}
data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of regression coefficients
  vector[n1]                      y1;            // event times in current data
  vector[n0]                      y0;            // event times in historical data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  array[n1] int<lower=1,upper=J>  intindx;       // index giving interval into which obs i failed / was censored for current data
  array[n0] int<lower=1,upper=J>  intindx0;      // index giving interval into which obs i failed / was censored for historical data
  array[n1] int<lower=0,upper=1>  death_ind;     // event indicator (1 = event; 0 = censored) for current data
  array[n0] int<lower=0,upper=1>  death_ind0;    // event indicator (1 = event; 0 = censored) for historical data
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  vector[p]                       beta_mean;     // mean for normal initial prior on coefficients
  vector<lower=0>[p]              beta_sd;       //sd for normal initial prior on coefficients
  vector[J]                       hazard_mean;   // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;     // scale parameter for half-normal prior on baseline hazards
  int<lower=0>                    s;             // number of a0s for which we have log nc
  vector<lower=0,upper=1>[s]      a0_lognc;
  vector[s]                       lognc;         // log normalizing constant for each a0 value specified in a0_lognc using the historical data
  real<lower=0>                   a0_shape1;
  real<lower=0>                   a0_shape2;
  real<lower=0,upper=1>           a0_lower;      // lower bounds for a0s
  real<lower=a0_lower,upper=1>    a0_upper;      // upper bounds for a0s
  int<lower=0,upper=1>            get_loglik;    // whether to generate log-likelihood matrix
}
transformed data{
  real lognc_hazard = normal_lccdf(0 | hazard_mean, hazard_sd);
  real logit_a0_lower    = logit(a0_lower);
  real logit_a0_upper    = logit(a0_upper);
  real lognc_logit_a0    = 0;

  if( a0_upper != 1 || a0_lower != 0 )
    lognc_logit_a0 = lognc_logit_a0 + log_diff_exp( beta_lcdf(a0_upper | a0_shape1, a0_shape2), beta_lcdf(a0_lower | a0_shape1, a0_shape2) );
}
parameters {
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  vector[p]             beta;         // regression coefficients
  real<lower=logit_a0_lower,upper=logit_a0_upper> logit_a0;
}
transformed parameters {
  real<lower=a0_lower,upper=a0_upper> a0 = inv_logit(logit_a0);
}
model {
  vector[n1]  eta        = X1 * beta;
  vector[n0]  eta0       = X0 * beta;
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;

  // compute lambda[j] * (s[j] * s[j-1])
  cumblhaz = cumulative_sum( lambda[1:(J-1)] .* ( breaks[2:J] - breaks[1:(J-1)] ) );

  // half-normal prior on baseline hazards
  target += normal_lpdf(lambda | hazard_mean, hazard_sd) - lognc_hazard;

  // noninformative initial prior
  target += normal_lpdf(beta | beta_mean, beta_sd);

  // prior on logit(a0)
  target += logit_beta_lpdf(logit_a0 | a0_shape1, a0_shape2) - lognc_logit_a0;

  // power prior
  for ( i in 1:n0 )
    target += a0 * pwe_lpdf(y0[i] | eta0[i], lambda, log_lambda, breaks, intindx0[i], death_ind0[i], cumblhaz);

  // current data likelihood
  for ( i in 1:n1 )
    target += pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);

  // Subtract log nc from power prior
  target += -pp_lognc(a0, a0_lognc, lognc);
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

