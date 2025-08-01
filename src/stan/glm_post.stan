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
  real lp2mean(real eta, int link) {
    if (link == 1) return eta;                        // identity link
    else if (link == 2) return exp(eta);              // log link
    else if (link == 3) return inv_logit(eta);        // logit link
    else if (link == 4) return inv(eta);              // inverse link
    else if (link == 5) return Phi_approx(eta);       // probit link
    else if (link == 6) return atan(eta) / pi() + 0.5; // cauchit link
    else if (link == 7) return inv_cloglog(eta);      // complementary log-log link
    else if (link == 8) return square(eta);           // sqrt link
    else if (link == 9) return inv_sqrt(eta);         // 1 / mu^2 link
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
  real normal_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if ( link != 1 )
      theta = lp2mean(theta, link);
    return normal_lpdf(y | theta, sqrt(phi));
  }

  real bernoulli_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 3 )
      theta = logit( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( log1p_exp(theta) );
  }
  real bernoulli_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if ( link != 3 )
      theta = logit( lp2mean(theta, link) );
    return y * theta - log1p_exp(theta);
  }

  real poisson_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    vector[n] theta = X * beta + offs;
    if ( link != 2 )
      theta = log( lp2mean(theta, link) );
    return dot_product(y, theta) - sum( exp(theta) + lgamma(y + 1) );
  }
  real poisson_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if ( link != 2 )
      theta = log( lp2mean(theta, link) );
    return y * theta - exp(theta) - lgamma(y + 1);
  }

  real gamma_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n           = rows(y);
    real tau        = inv(phi); // shape parameter
    vector[n] theta = X * beta + offs;
    if ( link != 4 )
      theta = inv( lp2mean(theta, link) );
    return gamma_lpdf(y | tau, tau * theta);
  }
  real gamma_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    real tau = inv(phi); // shape parameter
    if ( link != 4 )
      theta = inv( lp2mean(theta, link) );
    return gamma_lpdf(y | tau, tau * theta);
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
  real invgauss_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta   = dot_product(x, beta) + offs;
    real tau     = inv(phi); // shape parameter
    real log_2pi = 1.837877066409345483560659;  // log(2*pi)
    if ( link != 9 )
      theta = inv_square( lp2mean(theta, link) );
    return 0.5 * (
            log(tau) - log_2pi - 3 * log(y)
            - tau * square( (y * sqrt(theta) - 1) ) * inv(y)
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
  real glm_lpdf(real y, vector beta, real phi, row_vector x, int dist, int link, real offs) {
    // Compute likelihood
    if (dist == 1) {     // Bernoulli
      return bernoulli_glm_lpdf(y | beta, phi, x, link, offs);
    }
    else if (dist == 2) {  // Poisson
      return poisson_glm_lpdf(y | beta, phi, x, link, offs);
    }
    else if (dist == 3) {  // Normal
      return normal_glm_lpdf(y | beta, phi, x, link, offs);
    }
    else if (dist == 4) { // Gamma
      return gamma_glm_lpdf(y | beta, phi, x, link, offs);
    }
    else if (dist == 5) { // Inverse-Gaussian
      return invgauss_glm_lpdf(y | beta, phi, x, link, offs);
    }
    else reject("Distribution not supported");
    return 0; // never reached;
  }
}
data {
  int<lower=0>                        n; // number of observations in current data
  int<lower=0>                        p;
  vector[n]                           y; // response for current data
  matrix[n,p]                         X; // design matrix for current data
  vector[p]                           mean_beta; // mean for normal prior on coefficients
  vector<lower=0>[p]                  sd_beta; //sd for normal prior on coefficients
  real                                disp_mean; // mean for the half-normal prior on dispersion parameter
  real<lower=0>                       disp_sd; // sd for the half-normal prior on dispersion parameter
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[n]                           offs; // offset
  int<lower=0,upper=1>                get_loglik; // whether to generate log-likelihood matrix
}
transformed data{
  real lognc_disp = normal_lccdf(0 | disp_mean, disp_sd);
}
parameters {
  vector[p]     beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
}
model {
  target+= normal_lpdf(beta | mean_beta, sd_beta); // prior on beta ~ N(mu0, sd0)
  if ( dist <= 2 ) {
    target += glm_lp(y, beta, 1.0, X, dist, link, offs); // current data likelihood
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;  // half-normal prior for dispersion
    target += glm_lp(y, beta, dispersion[1], X, dist, link, offs); // current data likelihood
  }
}
generated quantities{
  vector[(get_loglik == 1) ? n : 0] log_lik; // to generate log likelihood matrices

  if ( get_loglik == 1 ) {
    real phi_val = (dist <= 2) ? 1.0 : dispersion[1];
    for( i in 1:n ){
      log_lik[i] = glm_lpdf(y[i] | beta, phi_val, X[i, ], dist, link, offs[i]);
    }
  }
}

