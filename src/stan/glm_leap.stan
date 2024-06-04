functions{
  //' Compute mean from linear predictor in a GLM
  //'
  //' @param eta linear predictor
  //' @param link integer giving link function
  vector lp2mean(vector eta, int link) {
    if (link == 1) return(eta);                        // identity link
    else if (link == 2) return exp(eta);               // log link
    else if (link == 3) return inv_logit(eta);         // logit link
    else if (link == 4) return inv(eta);               // inverse link
    else if (link == 5) return Phi_approx(eta);        // probit link
    else if (link == 6) return atan(eta) / pi() + 0.5; // cauchit link
    else if (link == 7) return inv_cloglog(eta);       // complementary log-log link
    else if (link == 8) return square(eta);            // sqrt link
    else if (link == 9) return inv_sqrt(eta);          // 1/mu^2 link
    else reject("Link not supported");
    return eta; // never reached
  }

  real normal_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 1 )
      theta = lp2mean(theta, link);
    return normal_lpdf(y | theta, sqrt(phi) );
  }

  real bernoulli_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 3 )
      theta = logit( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( log1p_exp(theta) );
  }

  real poisson_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 2 )
      theta = log( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( exp(theta) + lgamma(y + 1) );
  }

  real gamma_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    real tau        = inv(phi); // shape parameter
    vector[n] theta = X * beta + offs;
    if ( link != 4 )
      theta = inv( lp2mean(theta, link) );
    return gamma_lpdf(y | tau, tau * theta );
  }

  real invgauss_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n                 = rows(y);
    real tau              = inv(phi); // shape parameter
    real log_2pi          = 1.837877066409345483560659;  // log(2*pi)
    vector[n] theta       = X * beta + offs;
    if ( link != 9 )
      theta = inv_square( lp2mean(theta, link) );
    return 0.5 * (
              n * (log(tau) - log_2pi) - 3 * sum(log(y))
            - tau * dot_self( (y .* sqrt(theta) - 1) .* inv_sqrt(y) )
          );
  }

  real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offs) {
    // Compute likelihood
    if (dist == 1) {     // Bernoulli
      return bernoulli_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 2) {  // Poisson
      return poisson_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 3) {  // Normal
      return normal_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 4) { // Gamma
      return gamma_glm_lp(y, beta, phi, X, link, offs);
    }
    else if (dist == 5) { // Inverse-Gaussian
      return invgauss_glm_lp(y, beta, phi, X, link, offs);
    }
    else reject("Distribution not supported");
    return 0; // never reached;
  }

  //' Compute mean from linear predictor in a GLM
  //'
  //' @param eta linear predictor
  //' @param link integer giving link function
  matrix lp2mean(matrix eta, int link) {
    if (link == 1) return(eta);                        // identity link
    else if (link == 2) return exp(eta);               // log link
    else if (link == 3) return inv_logit(eta);         // logit link
    else if (link == 4) return inv(eta);               // inverse link
    else if (link == 5) return Phi_approx(eta);        // probit link
    else if (link == 6) return atan(eta) / pi() + 0.5; // cauchit link
    else if (link == 7) return inv_cloglog(eta);       // complementary log-log link
    else if (link == 8) return square(eta);            // sqrt link
    else if (link == 9) return inv_sqrt(eta);          // 1/mu^2 link
    else reject("Link not supported");
    return eta; // never reached
  }

  //' Log likelihood contribution to mixture of normal GLMs
  //'
  //' Computes the log likelihood contribution of each individual to each component
  //' of a mixture of normal GLMs
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //' @param offs offset
  //'
  //' @return n x K matrix of log likelihood contributions
  matrix normal_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
  ) {
    // compute logarithm of normalizing constant
    real log_2pi = 1.837877066409345483560659;  // log(2*pi)
    int n = rows(y);
    int p = cols(X);
    int K = rows(probs);
    matrix [n,K] theta;
    vector[K] log_probs = log(probs);
    vector[K] inv_disp = inv(disp);
    matrix[n,K] contrib;
    vector[n] y_sq = square(y);

    // Compute canonical parameter: theta = mu
    theta = X * beta + offs;
    if ( link != 1 )
      theta = lp2mean(theta, link);

    // Likelihood contribution:
    //  log f(y | theta, phi) = y * theta - b(theta) + c(y, phi)
    //     b(theta) = 0.5 * theta^2
    //     c(y, phi) = -0.5 * (y^2 / phi + log(2 * pi * phi))
    for ( k in 1:K ) {
      contrib[, k] =
          log_probs[k]
        + inv_disp[k] * (y .* theta[, k] - 0.5 * square(theta[, k]))
        - 0.5 * ( y_sq * inv_disp[k] + log_2pi + log(disp[k]) );
    }

    // Return n x K matrix of contributions to mixture likelihood
    return contrib;
  }
  //' Log likelihood contribution to mixture of Bernoulli GLMs
  //'
  //' Computes the log likelihood contribution of each individual to each component
  //' of a mixture of Bernoulli GLMs
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters (ignored for the likelihood calculation)
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //' @param offs offset
  //'
  //' @return n x K matrix of log likelihood contributions
  matrix bernoulli_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
  ) {
    int n = rows(y);
    int p = cols(X);
    int K = rows(probs);
    matrix [n,K] theta;
    vector[K] log_probs = log(probs);
    matrix[n,K] contrib;

    // Compute canonical parameter: theta = logit(mu)
    theta = X * beta + offs;
    if ( link != 3)
      theta = logit( lp2mean(theta, link) );

    // Likelihood contribution:
    //  log f(y | theta) = y * theta - b(theta) + c(y, phi)
    //    b(theta) = log(1 + exp(theta))
    //    c(y, phi) = 0
    for ( k in 1:K )
      contrib[, k] = log_probs[k] + y .* theta[, k] - log1p_exp(theta[, k]);

    // Return n x K matrix of contributions to mixture likelihood
    return contrib;
  }
  //' Log likelihood contribution to mixture of Poisson GLMs
  //'
  //' Computes the log likelihood contribution of each individual to each component
  //' of a mixture of Poisson GLMs
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters (ignored for the likelihood calculation)
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //'
  //' @return n x K matrix of log likelihood contributions
  matrix poisson_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
  ) {
    int n = rows(y);
    int p = cols(X);
    int K = rows(probs);
    matrix [n,K] theta;
    vector[K] log_probs = log(probs);
    matrix[n,K] contrib;
    vector[n] log_y_factorial = lgamma(y+1);

    // Compute canonical parameter: theta = log(mu)
    theta = X * beta + offs;
    if ( link != 2)
      theta = log( lp2mean(theta, link) );

    // Likelihood contribution:
    //  log f(y | theta) = y * theta - b(theta) + c(y, phi)
    //    b(theta)  = exp(theta)
    //    c(y, phi) = -log(y!)
    for ( k in 1:K )
      contrib[, k] = log_probs[k] + y .* theta[, k] - exp(theta[, k]) - log_y_factorial;

    // Return n x K matrix of contributions to mixture likelihood
    return contrib;
  }

  //' Log likelihood contribution to mixture of Gamma GLMs
  //'
  //' Computes the log likelihood contribution of each individual to each component
  //' of a mixture of Poisson GLMs
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //' @param offs offset
  //'
  //' @return n x K matrix of log likelihood contributions
  matrix gamma_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
    ) {
    int n = rows(y);
    int p = cols(X);
    int K = rows(probs);
    matrix [n,K] theta;
    vector[K] log_probs = log(probs);
    matrix[n,K] contrib;
    vector[n] log_y = log(y);
    vector[K] inv_disp = inv(disp);

    // Compute canonical parameter: theta = 1 / mu
    theta = X * beta + offs;
    if ( link != 4 )
      theta = inv( lp2mean(theta, link) );

    // Likelihood contribution:
    //  log f(y | theta) = 1 / phi * [ -y * theta - b(theta) ] + c(y, phi)
    //    b(theta)  = -log(theta)
    //    c(y, phi) = (1/phi - 1) * log(y) - (1 / phi) * log(phi) - lgamma(1 / phi)
    for ( k in 1:K )
      contrib[, k] = log_probs[k]
         + inv_disp[k] * ( -y .* theta[, k] + log(theta[, k]) )
         + (inv_disp[k] - 1) * log_y - inv_disp[k] * log(disp[k]) - lgamma(inv_disp[k]);
      ;

    // Return n x K matrix of contributions to mixture likelihood
    return contrib;
  }


  //' Log likelihood contribution to mixture of inverse-Gaussian GLMs
  //'
  //' Computes the log likelihood contribution of each individual to each
  //' component of a mixture of inverse-Gaussian GLMs
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //' @param offs offset
  //'
  //' @return n x K matrix of log likelihood contributions
  matrix invgauss_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
  ) {
    real log_2pi = 1.837877066409345483560659;  // log(2*pi)
    int n = rows(y);
    int p = cols(X);
    int K = rows(probs);
    matrix [n,K] theta;
    vector[K] log_probs = log(probs);
    matrix[n,K] contrib;
    vector[n] log_y_cubed = 3 * log(y);
    vector[n] inv_y = inv(y);
    vector[K] inv_disp = inv(disp);

    // Compute canonical parameter: theta = 1 / mu^2
    //  technically it is -1 / (2*mu^2)
    theta = X * beta + offs;
    if ( link != 9 )
      theta = inv_square( lp2mean(theta, link) );

    // Likelihood contribution:
    //  log f(y | theta) = 1 / phi * [ -0.5 * y * theta - b(theta) ] + c(y, phi)
    //    b(theta)  = -sqrt(theta)
    //    c(y, phi) = -0.5 * ( 1 / (y * phi) + log(phi) + log(2 * pi) + log(y^3) )
    for ( k in 1:K )
      contrib[, k] = log_probs[k]
         + inv_disp[k] * ( -0.5 * y .* theta[, k] + sqrt(theta[, k]) )
         - 0.5 * ( inv_disp[k] * inv_y + log(disp[k]) + log_2pi + log_y_cubed )
      ;

    // Return n x K matrix of contributions to mixture likelihood
    return contrib;
  }
  //' Get log likelihood contribution for each individual for each component
  //'
  //' Computes a n x K matrix of log likelihood contributions for a mixture of
  //' generalized linear models
  //'
  //' @param y vector of responses
  //' @param beta p x K matrix of regression coefficients
  //' @param phi K-dimensional vector of dispersion parameters
  //' @param probs K-dimensional vector of mixture component probabilities
  //' @param matrix X n x p design matrix
  //' @param dist index pertaining to distribution
  //' @param link index pertaining to link function
  //' @param offs n-dimensional vector of offsets
  //'
  //' @return n x K matrix of log likelihood contributions:
  //'   res[i, k] = log density for subject i in class k
  matrix glm_mixture_contrib(
    vector y, matrix beta, vector disp, vector probs
    , matrix X, int dist, int link, matrix offs
  ) {
    if (dist == 1)
      return bernoulli_glm_mixture_contrib(y, X, beta, disp, probs, link, offs);
    else if (dist == 2)
      return poisson_glm_mixture_contrib(y, X, beta, disp, probs, link, offs);
    else if (dist == 3)
      return normal_glm_mixture_contrib(y, X, beta, disp, probs, link, offs);
    else if (dist == 4)
      return gamma_glm_mixture_contrib(y, X, beta, disp, probs, link, offs);
    else if (dist == 5)
      return invgauss_glm_mixture_contrib(y, X, beta, disp, probs, link, offs);
    else reject("Distribution not supported");
    return beta; // never reached
  }


  //' Compute log likelihood for a mixture of generalized linear models
  //'
  //' Sum of log likelihood contributions from each individual
  //'
  //' @param y vector of responses
  //' @param X design matrix (incl. intercept if applicable)
  //' @param beta matrix of regression coefficients
  //' @param disp vector of dispersion parameters
  //' @param probs vector of component probabilities
  //' @param link index of link function
  //' @param offs offset
  //'
  //' @return sum of log likelihood (a real number)
  real glm_mixture_lp(
    vector y, matrix beta, vector disp, vector probs
    , matrix X, int dist, int link, matrix offs
  ) {
    int n = rows(y);
    int K = rows(probs);
    vector[n] contrib; // log probability summing over all components
    matrix[n, K] contribs = glm_mixture_contrib(y, beta, disp, probs, X, dist, link, offs);

    for (i in 1:n) {
      contrib[i] = log_sum_exp(contribs[i,]);
    }
    return sum(contrib);
  }
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
  matrix[p,K]           sd_beta;                // sd for normal initial prior on coefficients
  vector[K]             disp_mean;              // mean for the half-normal prior for dispersion of each dataset
  vector<lower=0>[K]    disp_sd;                // sd for the half-normal prior for dispersion of each dataset
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
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;  // K-dim vector of dispersion params
  real<lower=gamma_lower,upper=gamma_upper> gamma;  // probability of being exchangeable
  simplex[K-1] delta_raw;
}
transformed parameters {
  simplex[K] probs;
  // Compute probability of being in first component and marginalized log probability
  probs[1]   = gamma;
  probs[2:K] = (1 - gamma) * delta_raw;
}
model {
  // initial priors
  for ( k in 1:K ) {
    target += normal_lpdf(betaMat[, k] | mean_beta[,k], sd_beta[,k]);
  }
  if ( dist <= 2 ) {
    // historical data likelihood
    target += glm_mixture_lp(y0, betaMat, ones_vector(K), probs, X0, dist, link, offs0);
    // current data likelihood
    target += glm_lp(y, betaMat[, 1], 1.0, X, dist, link, offs[,1]);
  }
  else {
    // historical data likelihood
    target += glm_mixture_lp(y0, betaMat, dispersion, probs, X0, dist, link, offs0);
    // half-normal prior for dispersion
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - normal_lccdf(0 | disp_mean, disp_sd);
    // current data likelihood
    target += glm_lp(y, betaMat[, 1], dispersion[1], X, dist, link, offs[,1]);
  }

  // If two components, get beta prior on gamma;
  // If >2 components, get a dirichlet prior on raw delta
  if ( gamma_shape1 != 1 || gamma_shape2 != 1)
    target += beta_lpdf(gamma | gamma_shape1, gamma_shape2);
  if (K_gt_2)
    target += dirichlet_lpdf(delta_raw | conc_delta);
}
generated quantities {
  vector[p] beta = betaMat[, 1];
}

