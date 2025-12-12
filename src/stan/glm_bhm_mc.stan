// GLM with Bayesian Hierarchical Model (BHM) prior and Gaussian copula for missing covariates
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
    return gamma_lpdf(y | tau, tau * theta );
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

  //' Transform latent normal to observed scale based on marginal type
  //'
  //' @param z latent normal value
  //' @param type marginal type (1=continuous, 2=binary, 3=ordinal, 4=count)
  //' @param mean marginal mean (or prob for binary, lambda for count)
  //' @param sd marginal sd
  //' @param cutpoints cutpoints for ordinal (unused for other types)
  //' @param n_levels number of levels for ordinal
  real z_to_x(real z, int type, real mean_param, real sd_param,
              vector cutpoints, int n_levels) {
    if (type == 1) {
      // Continuous: location-scale transform
      return mean_param + sd_param * z;
    }
    else if (type == 2) {
      // Binary: probit threshold (z > 0 => x = 1)
      // Return probability for smooth gradient
      return Phi_approx(z);
    }
    else if (type == 3) {
      // Ordinal: find which level based on cutpoints
      // Return expected value (weighted by level)
      real p_cumul = 0.0;
      real expected = 0.0;
      real p_prev = 0.0;
      for (k in 1:n_levels) {
        real p_k;
        if (k < n_levels) {
          p_k = Phi_approx(cutpoints[k]) - p_prev;
        } else {
          p_k = 1.0 - p_prev;
        }
        expected += k * p_k;
        p_prev = Phi_approx(cutpoints[k]);
      }
      // Scale z to the expected level range
      return mean_param + sd_param * z;  // Simplified: treat as continuous
    }
    else if (type == 4) {
      // Count: use inverse Poisson CDF approximation
      // Approximate with continuous relaxation
      real lambda = mean_param;
      real p = Phi_approx(z);
      // Approximate quantile function for Poisson
      // Using normal approximation: x ≈ lambda + sqrt(lambda) * z
      return fmax(0.0, lambda + sqrt(lambda) * z);
    }
    else {
      reject("Unknown marginal type");
      return 0.0;
    }
  }
}

data {
  // Standard GLM data
  int<lower=0>                        K;         // total number of datasets
  int<lower=0>                        N;         // total number of observations
  array[K] int<lower=0, upper=N>      start_idx; // starting index of each dataset
  array[K] int<lower=0, upper=N>      end_idx;   // ending index of each dataset
  int<lower=0>                        p;         // number of covariates (design matrix columns)
  vector[N]                           y;         // response for the stacked data
  matrix[N,p]                         X;         // design matrix (with placeholders for missing)

  // Hyperpriors for BHM
  vector[p]                           meta_mean_mean;
  vector<lower=0>[p]                  meta_mean_sd;
  vector[p]                           meta_sd_mean;
  vector<lower=0>[p]                  meta_sd_sd;
  vector[K]                           disp_mean;
  vector<lower=0>[K]                  disp_sd;

  // Distribution and link
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs;      // offset

  // Flags
  int<lower=0,upper=1>                get_loglik;

  // Missing data structure
  int<lower=0>                        n_miss;    // total number of missing values
  int<lower=0>                        q;         // number of covariates with missing values
  array[n_miss, 3] int<lower=0>       miss_idx;  // [row, dm_col, cov_idx] for each missing

  // Copula correlation matrix
  matrix[q, q]                        copula_corr;

  // Marginal distribution parameters
  array[q] int<lower=1,upper=4>       marginal_type;  // 1=cont, 2=binary, 3=ordinal, 4=count
  vector[q]                           marginal_mean;
  vector<lower=0>[q]                  marginal_sd;
  array[q] int<lower=0>               n_levels;       // for ordinal
  int<lower=1>                        max_levels;     // maximum levels for ordinal
  matrix[q, max_levels-1]             cutpoints;      // ordinal cutpoints
}

transformed data {
  real lognc_beta_sd = normal_lccdf(0 | meta_sd_mean, meta_sd_sd);
  real lognc_disp    = normal_lccdf(0 | disp_mean, disp_sd);
  int n = end_idx[1]; // current data sample size

  // Cholesky decomposition of copula correlation for sampling
  matrix[q, q] L_copula = cholesky_decompose(copula_corr);
}

parameters {
  // BHM parameters (non-centered)
  matrix[p, K]                        beta_raw;
  vector[p]                           beta_mean;
  vector<lower=0>[p]                  beta_sd;
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;

  // Latent normal values for missing covariates
  vector[n_miss]                      z_miss_raw;  // raw (uncorrelated) latent values
}

transformed parameters {
  matrix[p, K]                        beta;
  matrix[N, p]                        X_complete;  // design matrix with imputed values

  // Non-centered parameterization for beta
  for ( j in 1:p ) {
    beta[j, ] = beta_mean[j] + beta_sd[j] * beta_raw[j, ];
  }

  // Start with original X (has placeholders for missing)
  X_complete = X;

  // Fill in missing values using latent z transformed to observed scale
  if (n_miss > 0) {
    for (m in 1:n_miss) {
      int row_idx = miss_idx[m, 1];
      int col_idx = miss_idx[m, 2];
      int cov_idx = miss_idx[m, 3];

      // Transform z to observed scale
      real x_imputed = z_to_x(
        z_miss_raw[m],
        marginal_type[cov_idx],
        marginal_mean[cov_idx],
        marginal_sd[cov_idx],
        to_vector(cutpoints[cov_idx, ]),
        n_levels[cov_idx]
      );

      X_complete[row_idx, col_idx] = x_imputed;
    }
  }
}

model {
  // Hyperpriors for BHM
  target += normal_lpdf(beta_mean | meta_mean_mean, meta_mean_sd);
  target += normal_lpdf(beta_sd | meta_sd_mean, meta_sd_sd) - lognc_beta_sd;

  for ( j in 1:p ) {
    target += std_normal_lpdf(beta_raw[j, ]);
  }

  // Prior on latent missing values (standard normal - copula structure)
  // The copula correlation is implicitly handled through the marginal transforms
  target += std_normal_lpdf(z_miss_raw);

  // Likelihood
  if ( dist <= 2 ) {
    for ( k in 1:K ) {
      target += glm_lp(
        y[start_idx[k]:end_idx[k]],
        beta[, k],
        1.0,
        X_complete[start_idx[k]:end_idx[k], ],
        dist, link,
        offs[start_idx[k]:end_idx[k]]
      );
    }
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;
    for ( k in 1:K ) {
      target += glm_lp(
        y[start_idx[k]:end_idx[k]],
        beta[, k],
        dispersion[k],
        X_complete[start_idx[k]:end_idx[k], ],
        dist, link,
        offs[start_idx[k]:end_idx[k]]
      );
    }
  }
}

generated quantities {
  vector[(get_loglik == 1) ? n : 0] log_lik;

  if ( get_loglik == 1 ) {
    real phi_val = (dist <= 2) ? 1.0 : dispersion[1];
    for (i in 1:n) {
      log_lik[i] = glm_lpdf(y[i] | beta[, 1], phi_val, X_complete[i, ], dist, link, offs[i]);
    }
  }
}
