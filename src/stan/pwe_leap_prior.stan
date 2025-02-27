// fit PWE model
// sample from the LEAP (Not posterior)
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

  //' compute log likelihood for a mixture of PWE models
  //'
  real pwe_mixture_lpdf(
    real y, row_vector eta, matrix lambdaMat, matrix log_lambdaMat, vector breaks, int j, int death_ind, matrix cumblhazMat, vector log_probs
  ) {
    int K = rows(log_probs);
    vector[K] contribs; // log likelihood contribution of each component

    for (k in 1:K) {
      contribs[k] = log_probs[k] + pwe_lpdf(y | eta[k], lambdaMat[, k], log_lambdaMat[, k], breaks, j, death_ind, cumblhazMat[, k]);
    }

    return log_sum_exp(contribs);
  }

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return(
      -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x)
    );
  }
}
data {
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of regression coefficients
  int<lower=0>                    K;             // number of components in mixture
  vector[n0]                      y0;            // event times in historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  array[n0] int<lower=1,upper=J>  intindx0;      // index giving interval into which obs i failed / was censored for historical data
  array[n0] int<lower=0,upper=1>  death_ind0;    // event indicator (1 = event; 0 = censored) for historical data
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  vector<lower=0>[K]              conc;          // concentration parameter for exchangeability
  real<lower=0,upper=1>           gamma_lower;   // lower bound for probability of being exchangeable
  real<lower=gamma_lower,upper=1> gamma_upper;   // upper bound for probability of being exchangeable
  vector[p]                       beta_mean;     // mean for normal initial prior on coefficients
  vector<lower=0>[p]              beta_sd;       // sd for normal initial prior on coefficients
  vector[J]                       hazard_mean;   // location parameter for half-normal prior on baseline hazards
  vector<lower=0>[J]              hazard_sd;     // scale parameter for half-normal prior on baseline hazards
}
transformed data {
  real gamma_shape1      = conc[1];
  real gamma_shape2      = sum(conc[2:K]);
  real logit_gamma_lower = logit(gamma_lower);
  real logit_gamma_upper = logit(gamma_upper);
  int K_gt_2             = (K > 2) ? (1) : (0);
  vector[K-1] conc_delta = conc[2:K];
  real lognc_hazard      = normal_lccdf(0 | hazard_mean, hazard_sd);
  real lognc_logit_gamma = 0;

  if( gamma_upper != 1 || gamma_lower != 0 )
    lognc_logit_gamma = lognc_logit_gamma + log_diff_exp( beta_lcdf(gamma_upper | gamma_shape1, gamma_shape2), beta_lcdf(gamma_lower | gamma_shape1, gamma_shape2) );
}
parameters {
  matrix<lower=0>[J,K]  lambdaMat;    // JxK matrix of hazards for each interval
  matrix[p,K]           betaMat;      // pxK matrix of regression coefficients; p = number of covars, K = number of components
  simplex[K-1]          delta_raw;
  real<lower=logit_gamma_lower,upper=logit_gamma_upper> logit_gamma; // logit of probability of being exchangeable
}
model {
  matrix[n0, K] eta0Mat       = X0 * betaMat;
  matrix[J,K]   log_lambdaMat = log(lambdaMat);
  vector[K]     logprobs      = append_row(logit_gamma, log(delta_raw)) - log1p_exp(logit_gamma);
  matrix[J-1,K] cumblhazMat;

  // noninformative initial prior
  for ( k in 1:K ) {
    // normal initial prior on betaMat
    target += normal_lpdf(betaMat[, k] | beta_mean, beta_sd);
    // half-normal prior on baseline hazards
    target += normal_lpdf(lambdaMat[, k] | hazard_mean, hazard_sd) - lognc_hazard;

    // compute lambda[j] * (s[j] * s[j-1])
    cumblhazMat[, k] = cumulative_sum( lambdaMat[1:(J-1), k] .* ( breaks[2:J] - breaks[1:(J-1)] ) );
  }

  // historical data likelihood
  for ( i in 1:n0 )
    target += pwe_mixture_lpdf(y0[i] | eta0Mat[i, ], lambdaMat, log_lambdaMat, breaks, intindx0[i], death_ind0[i], cumblhazMat, logprobs);

  // If two components, get beta prior on gamma;
  // If >2 components, get a dirichlet prior on raw delta
  target += logit_beta_lpdf(logit_gamma | gamma_shape1, gamma_shape2) - lognc_logit_gamma;

  if (K_gt_2)
    target += dirichlet_lpdf(delta_raw | conc_delta);
}

