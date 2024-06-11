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
  int<lower=0>                        N; // total number of observations (excluding the current data)
  array[K] int<lower=0, upper=N>      start_idx; // starting index of each data in the stacked version
  array[K] int<lower=0, upper=N>      end_idx; // ending index of each data in the stacked version
  int<lower=0>                        p;
  vector[N]                           y; // response for the stacked data
  matrix[N,p]                         X; // design matrix for the stacked data
  vector[p]                           beta0_mean;
  vector<lower=0>[p]                  beta0_sd;
  vector[K]                           disp_mean;    // location parameter for the half-normal prior on dispersion of each historical dataset
  vector<lower=0>[K]                  disp_sd;      // scale parameter for the half-normal prior on dispersion of each historical dataset
  real<lower=0,upper=1>               p_spike;      // prior probability of spike
  real                                mu_spike;     // location parameter for half-normal prior (spike)
  real<lower=0>                       sigma_spike;  // scale parameter for half-normal prior (spike)
  real                                mu_slab;      // location parameter for half-normal prior (slab)
  real<lower=0>                       sigma_slab;   // scale parameter for half-normal prior (slab)
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
transformed data{
  real lognc_spike = normal_lccdf(0 | mu_spike, sigma_spike); // \Phi(mu_spike / sigma_spike)
  real lognc_slab  = normal_lccdf(0 | mu_slab, sigma_slab); // \Phi(mu_slab / sigma_slab)
  real lognc_disp  = normal_lccdf(0 | disp_mean, disp_sd);
}
parameters {
  vector[p] beta;  // regression coefficients for current data
  vector[p] beta0; // regression coefficients for historical datasets
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;
  vector<lower=0>[p] comm_prec;  // commensurability parameter
}
transformed parameters {
  vector[p] comm_sd = inv_sqrt(comm_prec);   // sd for commensurate prior
}
model {
  // spike and slab prior for commensurability:
  for ( i in 1:p ) {
    target += log_mix(
      p_spike
      , normal_lpdf(comm_prec[i] | mu_spike, sigma_spike) - lognc_spike
      , normal_lpdf(comm_prec[i] | mu_slab, sigma_slab) - lognc_slab
    );
  }
  // initial prior for beta0
  target += normal_lpdf(beta0 | beta0_mean, beta0_sd);
  // commensurate prior for beta
  target += normal_lpdf(beta  | beta0, comm_sd);

  if ( dist <= 2 ) {
    target += glm_lp(y, beta0, 1.0, X, dist, link, offs); // historical data likelihood
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;  // half-normal prior for dispersion
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta0, dispersion[k], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // historical data likelihood
    }
  }
}

