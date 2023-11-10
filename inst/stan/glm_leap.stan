functions{
#include include/expfam_loglik.stan
#include include/leap_helpers.stan
}
data {
  int<lower=0>          n;                      // sample size of current data
  int<lower=0>          n0;                     // sample size of historical data
  int<lower=0>          p;                      // number of covariates (incl intercept)
  int<lower=0>          K;                      // number of components in mixture
  vector[n]             y;                      // response vector for current data
  matrix[n,p]           X;                      // design mtx for current data
  vector[n0]            y0;                     // response vector for historical data
  matrix[n0,p]          X0;                     // design mtx for historical data
  matrix[p,K]           mean_beta;              // mean for normal initial prior on coefficients
  matrix[p*K,p]         cov_beta;               // covariance mtx for normal initial prior on coefficients
  real<lower=0>         disp_shape;             // shape parameter for inv-gamma initial prior on dispersion
  real<lower=0>         disp_scale;             // scale parameter for inv-gamma initial prior on dispersion
  vector<lower=0>[K]    conc;                   // Concentration parameter for exchangeability
  real<lower=0,upper=1> gamma_lower;            // Lower bound for probability of being exchangeable.
  real<lower=gamma_lower,upper=1> gamma_upper;  // Upper bound for probability of being exchangeable.
  int<lower=1,upper=5>  dist;
  int<lower=1,upper=9>  link;
  matrix[n,K]           offs;
  matrix[n0,K]          offs0;
}
transformed data {
  real gamma_shape1 = conc[1];
  real gamma_shape2 = sum(conc[2:K]);
  int K_gt_2 = (K > 2) ? (1) : (0);
  vector[K-1] conc_delta = conc[2:K];
}
parameters {
  matrix[p,K] betaMat;       // pxK matrix of regression coefficients; p = number of covars, K = number of components
  vector<lower=0>[K] dispersion;  // K-dim vector of dispersion params
  real<lower=gamma_lower,upper=gamma_upper> gamma;  // probability of being exchangeable
  simplex[K-1] delta_raw;
}
transformed parameters {
  vector[n0] lp01;                           // log probability for first component
  vector[n0] contrib;                        // log probability summing over all components
  vector[K] probs;
  vector[K] logProbs;
  matrix[n0, K] contribs;

  probs[1]   = gamma;
  probs[2:K] = (1 - gamma) * delta_raw;
  logProbs = log(probs);
  // Compute probability of being in first component and marginalized log probability
  contribs = glm_mixture_contrib(y0, betaMat, dispersion, probs, X0, dist, link, offs0);

  for (i in 1:n0) contrib[i] = log_sum_exp(contribs[i,]);
}
model {
  if ( dist <= 2 ) {
    target += glm_lp(y, betaMat[, 1], 1.0, X, dist, link, offs[,1]); // current data likelihood
    target += contrib; // historical data likelihood
    // initial priors
    for ( k in 1:K ) {
      betaMat[, k] ~ multi_normal(mean_beta[,k], cov_beta[(1+(k-1)*p):(k*p),]);
    }
  }
  else {
    target += glm_lp(y, betaMat[, 1], dispersion[1], X, dist, link, offs[,1]); // current data likelihood
    target += contrib; // historical data likelihood
    // initial priors
    for ( k in 1:K ) {
      betaMat[, k] ~ multi_normal(mean_beta[,k], dispersion[k] * cov_beta[(1+(k-1)*p):(k*p),]);
      dispersion[k] ~ normal(disp_shape, disp_scale);
    }
  }
  // If two components, get beta prior on gamma;
  // If >2 components, get a dirichlet prior on raw delta
  if ( gamma_shape1 != 1 || gamma_shape2 != 1)
    gamma ~ beta(gamma_shape1, gamma_shape2);
  if (K_gt_2)
    delta_raw ~ dirichlet(conc_delta);
}
generated quantities {
  vector[p] beta = betaMat[, 1];
}
