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
}
data {
  int<lower=0>                        K; // total number of historical datasets (excluding the current data)
  int<lower=0>                        n; // number of observations in the current data
  int<lower=0>                        p;
  int<lower=0,upper=1>                ind_disp; // indicator for whether dispersion is included
  vector[n]                           y; // response for the current data
  matrix[n,p]                         X; // design matrix for the current data
  matrix[p + ind_disp, K]             theta_hats; // the kth column consists of the MLE for beta and LOG DISPERSION for the kth historical dataset
  array[K] cov_matrix[p + ind_disp]   theta_covars; // the kth element is the inverse fisher info for (beta, LOG PHI) for the kth historical dataset
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[n]                           offs; // offset for current data


}
parameters {
  vector[p] beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
  vector<lower=0,upper=1>[K] a0s;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // prior of a0s is beta (uniform if 1's)
  if ( a0_shape1 != 1 || a0_shape2 != 1)
    target += beta_lpdf(a0s | a0_shape1, a0_shape2);

  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offs); // current data likelihood

    // asympt. power prior
    for ( k in 1:K ) {
      target += multi_normal_lpdf(beta | theta_hats[, k], inv(a0s[k]) * theta_covars[k]);
    }
  }
  else {
    // change of variables:
    //   eta = (beta, phi) --> theta = (beta, log(phi))
    //   Jacobian       = blkdiag ( I, 1/phi )
    //   log |det(Jacobian)| = -log(phi)
    real log_disp      = log(dispersion[1]);
    vector[p+1] theta  = append_row( beta, log_disp );  // theta = ( beta, log(phi) )'
    for ( k in 1:K ) {
      target          += multi_normal_lpdf(theta | theta_hats[, k], inv(a0s[k]) * theta_covars[k]);
    }
    target            += -log_disp;   // jacobian

    // likelihood of current data
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offs);
  }
}

