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
  int<lower=0>                        K; // total number of datasets (including the current data)
  int<lower=0>                        N; // total number of observations (including the current data)
  array[K] int<lower=0, upper=N>      start_idx; // starting index of each data in the stacked version
  array[K] int<lower=0, upper=N>      end_idx; // ending index of each data in the stacked version
  int<lower=0>                        p;
  vector[N]                           y; // response for the stacked data
  matrix[N,p]                         X; // design matrix for the stacked data
  vector[p]                           meta_mean_mean; // mean for the hyperprior for mean of each coefficient
  vector<lower=0>[p]                  meta_mean_sd; // sd for the hyperprior for mean of each coefficient
  vector[p]                           meta_sd_mean; // mean for the hyperprior for sd of each coefficient
  vector<lower=0>[p]                  meta_sd_sd; // sd for the hyperprior for sd of each coefficient
  vector[K]                           disp_mean; // mean for the half-normal prior for dispersion of each dataset
  vector<lower=0>[K]                  disp_sd; // sd for the half-normal prior for dispersion of each dataset
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
parameters {
  //matrix[p, K]                        beta; // kth column is the vector of coefficients for the kth dataset
  matrix[p, K]                        beta_raw; // for non-centered parameterization
  vector[p]                           beta_mean;
  vector<lower=0>[p]                  beta_sd;
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;
}
transformed parameters{
  matrix[p, K]                        beta; // kth column is the vector of coefficients for the kth dataset
  for ( j in 1:p ) {
    beta[j, ] = beta_mean[j] + beta_sd[j] * beta_raw[j, ];
  }
}
model {
  beta_mean       ~ normal(meta_mean_mean, meta_mean_sd); // hyperprior for mean of coefficients
  beta_sd         ~ normal(meta_sd_mean, meta_sd_sd); // hyperprior for sd of coefficients (half-normal)
  for ( j in 1:p ) {
    beta_raw[j, ]   ~ std_normal(); // implies beta[j, ] ~ normal(beta_mean[j], beta_sd[j]);
  }

  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
  else {
    dispersion ~ normal(disp_mean, disp_sd); // half-normal prior for dispersion
    for ( k in 1:K ) {
      target += glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta[, k], dispersion[k], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);   // data likelihood
      }
  }
}

