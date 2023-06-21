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

//' Compute the density of a mixture of multivariate normal distributions
//'
//' @param x vector to evaluate density
//' @param probs vector giving mixing proportions
//' @param means matrix, the gth column giving means for the gth component of the mixture density
//' @param precisions array of precision matrices
real multi_normal_mix_lpdf(vector x, vector probs, matrix means, array[] matrix precisions) {
  int G = size(x);
  vector[G] log_contrib;
  for ( i in 1:G) {
    log_contrib[i] = log(probs[i]) + multi_normal_prec_lpdf(x | means[, i], precisions[i]);
  }
  return log_sum_exp(log_contrib);
}
