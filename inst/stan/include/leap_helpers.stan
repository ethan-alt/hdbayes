
//' Compute mean from linear predictor in a GLM
//'
//' @param eta linear predictor
//' @param link integer giving link function
matrix lp2mean(matrix eta, int link) {
  if (link == 1) return(eta);                        // identity link
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
//'
//' @return n x K matrix of log likelihood contributions
matrix normal_glm_mixture_contrib(
  vector y, matrix X, matrix beta, vector disp, vector probs, int link, matrix offs
) {
  // compute logarithm of normalizing constant
  real log_2pi = 1.837877066409345483560659;  // log(2*pi)
  int n = rows(X);
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
  int n = rows(X);
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
  int n = rows(X);
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
  int K = rows(probs);
  int n = rows(y);
  int p = cols(X);
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
  //  log f(y | theta) = 1 / phi * [ y * theta - b(theta) ] + c(y, phi)
  //    b(theta)  = -log(theta) [technically b(theta) = -log(-theta) but absorb negative in coefficients]
  //    c(y, phi) = (1/phi - 1) * log(y) - (1 / phi) * log(phi) - lgamma(1 / phi)
  for ( k in 1:K )
    contrib[, k] = log_probs[k]
       + inv_disp[k] * ( y .* theta[, k] + log(theta[, k]) )
       + (inv_disp[k] - 1) * log_y - inv_disp * log(disp[k]) - lgamma(inv_disp[k]);
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
  int K = rows(probs);
  int n = rows(y);
  int p = cols(X);
  matrix [n,K] theta;
  vector[K] log_probs = log(probs);
  matrix[n,K] contrib;
  vector[n] log_y_cubed = 3 * log(y);
  vector[n] inv_y = inv(y);
  vector[K] inv_disp = inv(disp);
  real log_2pi = 1.837877066409345483560659;  // log(2*pi)

  // Compute canonical parameter: theta = 1 / mu
  //  technically it is -1 / (2*mu^2), but negative can be absorbed in reg coefs
  theta = X * beta + offs;
  if ( link != 4 )
    theta = inv( lp2mean(theta, link) );

  // Likelihood contribution:
  //  log f(y | theta) = 1 / phi * [ 0.5 * y * theta - b(theta) ] + c(y, phi)
  //    b(theta)  = -sqrt(theta)
  //    c(y, phi) = -0.5 * ( 1 / (y * phi) + log(phi) + log(2 * pi) + log(y^3) )
  for ( k in 1:K )
    contrib[, k] = log_probs[k]
       + inv_disp[k] * ( y .* theta[, k] + sqrt(theta[, k]) )
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
  else reject("Distribution not supported")
  return beta; // never reached
}
