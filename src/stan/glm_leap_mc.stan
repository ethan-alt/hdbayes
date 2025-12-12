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

  //' Transform latent normal z to observed covariate scale
  //'
  //' @param z latent normal value (standard normal scale)
  //' @param type_id 1=continuous, 2=binary, 3=ordinal, 4=count
  //' @param mu marginal mean (for continuous) or prob (for binary) or lambda (for count)
  //' @param sigma marginal sd (for continuous)
  //' @param cutpoints vector of cutpoints for ordinal
  //' @param n_levels number of ordinal levels
  real z_to_x(real z, int type_id, real mu, real sigma,
              vector cutpoints, int n_levels) {
    if (type_id == 1) {
      // Continuous: location-scale transform
      return mu + sigma * z;
    } else if (type_id == 2) {
      // Binary: probit model threshold at 0
      return z > 0 ? 1.0 : 0.0;
    } else if (type_id == 3) {
      // Ordinal: find which interval z falls into
      for (k in 1:(n_levels - 1)) {
        if (z <= cutpoints[k])
          return k - 1.0;
      }
      return n_levels - 1.0;
    } else if (type_id == 4) {
      // Count: inverse CDF via Poisson
      // Use normal CDF to get probability, then find quantile
      real p = Phi(z);
      real lambda = mu;
      // Simple quantile approximation for Poisson
      if (p < 1e-6) return 0.0;
      if (p > 1 - 1e-6) return lambda + 5 * sqrt(lambda);
      // Newton-style search for quantile
      int k = 0;
      real cum_p = exp(-lambda);
      while (cum_p < p && k < 100) {
        k += 1;
        cum_p += exp(-lambda + k * log(lambda) - lgamma(k + 1));
      }
      return k;
    }
    return mu + sigma * z; // fallback
  }

  real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offs) {
    int n = rows(y);
    vector[n] theta = X * beta + offs;

    // Compute mean from linear predictor
    if (dist == 1) { // Bernoulli
      if (link != 3) theta = logit(lp2mean(theta, link));
      return dot_product(y, theta) - sum(log1p_exp(theta));
    }
    else if (dist == 2) { // Poisson
      if (link != 2) theta = log(lp2mean(theta, link));
      return dot_product(y, theta) - sum(exp(theta) + lgamma(y + 1));
    }
    else if (dist == 3) { // Normal
      if (link != 1) theta = lp2mean(theta, link);
      return normal_lpdf(y | theta, sqrt(phi));
    }
    else if (dist == 4) { // Gamma
      real tau = inv(phi);
      if (link != 4) theta = inv(lp2mean(theta, link));
      return gamma_lpdf(y | tau, tau * theta);
    }
    else if (dist == 5) { // Inverse-Gaussian
      real tau = inv(phi);
      real log_2pi = 1.837877066409345483560659;
      if (link != 9) theta = inv_square(lp2mean(theta, link));
      return 0.5 * (n * (log(tau) - log_2pi) - 3 * sum(log(y))
                    - tau * dot_self((y .* sqrt(theta) - 1) .* inv_sqrt(y)));
    }
    else reject("Distribution not supported");
    return 0;
  }

  real glm_lpdf(real y, vector beta, real phi, row_vector x, int dist, int link, real offs) {
    real theta = dot_product(x, beta) + offs;

    if (dist == 1) { // Bernoulli
      if (link != 3) theta = logit(lp2mean(theta, link));
      return y * theta - log1p_exp(theta);
    }
    else if (dist == 2) { // Poisson
      if (link != 2) theta = log(lp2mean(theta, link));
      return y * theta - exp(theta) - lgamma(y + 1);
    }
    else if (dist == 3) { // Normal
      if (link != 1) theta = lp2mean(theta, link);
      return normal_lpdf(y | theta, sqrt(phi));
    }
    else if (dist == 4) { // Gamma
      real tau = inv(phi);
      if (link != 4) theta = inv(lp2mean(theta, link));
      return gamma_lpdf(y | tau, tau * theta);
    }
    else if (dist == 5) { // Inverse-Gaussian
      real tau = inv(phi);
      real log_2pi = 1.837877066409345483560659;
      if (link != 9) theta = inv_square(lp2mean(theta, link));
      return 0.5 * (log(tau) - log_2pi - 3 * log(y)
                    - tau * square((y * sqrt(theta) - 1)) * inv(y));
    }
    else reject("Distribution not supported");
    return 0;
  }

  //' Compute mean from linear predictor in a GLM (matrix version)
  matrix lp2mean(matrix eta, int link) {
    if (link == 1) return(eta);
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv_logit(eta);
    else if (link == 4) return inv(eta);
    else if (link == 5) return Phi_approx(eta);
    else if (link == 6) return atan(eta) / pi() + 0.5;
    else if (link == 7) return inv_cloglog(eta);
    else if (link == 8) return square(eta);
    else if (link == 9) return inv_sqrt(eta);
    else reject("Link not supported");
    return eta;
  }

  //' Log likelihood contribution to mixture of normal GLMs
  matrix normal_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector log_probs, int link, matrix offs
  ) {
    real log_2pi = 1.837877066409345483560659;
    int n = rows(y);
    int K = rows(log_probs);
    matrix [n,K] theta;
    vector[K] inv_disp = inv(disp);
    matrix[n,K] contrib;
    vector[n] y_sq = square(y);

    theta = X * beta + offs;
    if (link != 1) theta = lp2mean(theta, link);

    for (k in 1:K) {
      contrib[, k] = log_probs[k]
        + inv_disp[k] * (y .* theta[, k] - 0.5 * square(theta[, k]))
        - 0.5 * (y_sq * inv_disp[k] + log_2pi + log(disp[k]));
    }
    return contrib;
  }

  //' Log likelihood contribution to mixture of Bernoulli GLMs
  matrix bernoulli_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector log_probs, int link, matrix offs
  ) {
    int n = rows(y);
    int K = rows(log_probs);
    matrix [n,K] theta;
    matrix[n,K] contrib;

    theta = X * beta + offs;
    if (link != 3) theta = logit(lp2mean(theta, link));

    for (k in 1:K)
      contrib[, k] = log_probs[k] + y .* theta[, k] - log1p_exp(theta[, k]);
    return contrib;
  }

  //' Log likelihood contribution to mixture of Poisson GLMs
  matrix poisson_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector log_probs, int link, matrix offs
  ) {
    int n = rows(y);
    int K = rows(log_probs);
    matrix [n,K] theta;
    matrix[n,K] contrib;
    vector[n] log_y_factorial = lgamma(y+1);

    theta = X * beta + offs;
    if (link != 2) theta = log(lp2mean(theta, link));

    for (k in 1:K)
      contrib[, k] = log_probs[k] + y .* theta[, k] - exp(theta[, k]) - log_y_factorial;
    return contrib;
  }

  //' Log likelihood contribution to mixture of Gamma GLMs
  matrix gamma_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector log_probs, int link, matrix offs
  ) {
    int n = rows(y);
    int K = rows(log_probs);
    matrix [n,K] theta;
    matrix[n,K] contrib;
    vector[n] log_y = log(y);
    vector[K] inv_disp = inv(disp);

    theta = X * beta + offs;
    if (link != 4) theta = inv(lp2mean(theta, link));

    for (k in 1:K)
      contrib[, k] = log_probs[k]
         + inv_disp[k] * (-y .* theta[, k] + log(theta[, k]))
         + (inv_disp[k] - 1) * log_y - inv_disp[k] * log(disp[k]) - lgamma(inv_disp[k]);
    return contrib;
  }

  //' Log likelihood contribution to mixture of inverse-Gaussian GLMs
  matrix invgauss_glm_mixture_contrib(
    vector y, matrix X, matrix beta, vector disp, vector log_probs, int link, matrix offs
  ) {
    real log_2pi = 1.837877066409345483560659;
    int n = rows(y);
    int K = rows(log_probs);
    matrix [n,K] theta;
    matrix[n,K] contrib;
    vector[n] log_y_cubed = 3 * log(y);
    vector[n] inv_y = inv(y);
    vector[K] inv_disp = inv(disp);

    theta = X * beta + offs;
    if (link != 9) theta = inv_square(lp2mean(theta, link));

    for (k in 1:K)
      contrib[, k] = log_probs[k]
         + inv_disp[k] * (-0.5 * y .* theta[, k] + sqrt(theta[, k]))
         - 0.5 * (inv_disp[k] * inv_y + log(disp[k]) + log_2pi + log_y_cubed);
    return contrib;
  }

  //' Get log likelihood contribution for mixture of GLMs
  matrix glm_mixture_contrib(
    vector y, matrix beta, vector disp, vector log_probs,
    matrix X, int dist, int link, matrix offs
  ) {
    if (dist == 1)
      return bernoulli_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs);
    else if (dist == 2)
      return poisson_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs);
    else if (dist == 3)
      return normal_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs);
    else if (dist == 4)
      return gamma_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs);
    else if (dist == 5)
      return invgauss_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs);
    else reject("Distribution not supported");
    return beta;
  }

  //' Compute log likelihood for a mixture of GLMs
  real glm_mixture_lp(
    vector y, matrix beta, vector disp, vector log_probs,
    matrix X, int dist, int link, matrix offs
  ) {
    int n = rows(y);
    int K = rows(log_probs);
    vector[n] contrib;
    matrix[n, K] contribs = glm_mixture_contrib(y, beta, disp, log_probs, X, dist, link, offs);

    for (i in 1:n) {
      contrib[i] = log_sum_exp(contribs[i,]);
    }
    return sum(contrib);
  }

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x);
  }
}

data {
  int<lower=0>          n;                      // sample size of current data
  int<lower=0>          n0;                     // sample size of historical data
  int<lower=0>          p;                      // number of covariates (incl intercept)
  int<lower=0>          K;                      // number of components in mixture
  vector[n]             y;                      // response vector for current data
  matrix[n,p]           X;                      // design mtx for current data (with placeholders for missing)
  vector[n0]            y0;                     // response vector for historical data
  matrix[n0,p]          X0;                     // design mtx for historical data (with placeholders for missing)
  matrix[p,K]           mean_beta;              // mean for normal initial prior on coefficients
  matrix<lower=0>[p,K]  sd_beta;                // sd for normal initial prior on coefficients
  vector[K]             disp_mean;              // mean for the half-normal prior for dispersion
  vector<lower=0>[K]    disp_sd;                // sd for the half-normal prior for dispersion
  vector<lower=0>[K]    conc;                   // Concentration parameter for exchangeability
  real<lower=0,upper=1> gamma_lower;            // Lower bound for probability of being exchangeable
  real<lower=gamma_lower,upper=1> gamma_upper;  // Upper bound for probability of being exchangeable
  int<lower=1,upper=5>  dist;
  int<lower=1,upper=9>  link;
  matrix[n,K]           offs;                   // offset matrix for current data
  matrix[n0,K]          offs0;                  // offset matrix for historical data
  int<lower=0,upper=1>  get_loglik;             // whether to generate log-likelihood matrix

  // Missing covariate data (Gaussian copula parameters)
  int<lower=0>          n_miss;                 // total number of missing values in current data
  int<lower=0>          n_miss0;                // total number of missing values in historical data
  int<lower=0>          q;                      // number of covariates that can have missing values
  array[n_miss > 0 ? n_miss : 1, 3] int miss_idx;   // [row, dm_col, cov_idx] for current data
  array[n_miss0 > 0 ? n_miss0 : 1, 3] int miss_idx0; // [row, dm_col, cov_idx] for historical data
  matrix[q > 0 ? q : 1, q > 0 ? q : 1] copula_corr; // correlation matrix for copula
  array[q > 0 ? q : 1] int marginal_type;       // 1=continuous, 2=binary, 3=ordinal, 4=count
  vector[q > 0 ? q : 1] marginal_mean;          // mean for continuous, prob for binary, lambda for count
  vector<lower=0>[q > 0 ? q : 1] marginal_sd;   // sd for continuous
  array[q > 0 ? q : 1] int n_levels;            // number of levels for ordinal
  matrix[q > 0 ? q : 1, 9] cutpoints;           // cutpoints for ordinal (max 10 levels)
  int<lower=1>          max_levels;
}

transformed data {
  real gamma_shape1      = conc[1];
  real gamma_shape2      = sum(conc[2:K]);
  real logit_gamma_lower = logit(gamma_lower);
  real logit_gamma_upper = logit(gamma_upper);
  int K_gt_2             = (K > 2) ? (1) : (0);
  vector[K-1] conc_delta = conc[2:K];
  real lognc_disp        = normal_lccdf(0 | disp_mean, disp_sd);
  real lognc_logit_gamma = 0;

  // Cholesky factor of copula correlation for sampling
  matrix[q > 0 ? q : 1, q > 0 ? q : 1] L_corr = cholesky_decompose(copula_corr);

  if (gamma_upper != 1 || gamma_lower != 0)
    lognc_logit_gamma = lognc_logit_gamma + log_diff_exp(
      beta_lcdf(gamma_upper | gamma_shape1, gamma_shape2),
      beta_lcdf(gamma_lower | gamma_shape1, gamma_shape2)
    );
}

parameters {
  matrix[p,K] betaMat;                                          // regression coefficients
  vector<lower=0>[(dist > 2) ? K : 0] dispersion;               // dispersion params
  real<lower=logit_gamma_lower,upper=logit_gamma_upper> logit_gamma;
  simplex[K-1] delta_raw;

  // Latent normal values for missing covariates (non-centered parameterization)
  vector[n_miss > 0 ? n_miss : 0] z_miss_raw;    // current data
  vector[n_miss0 > 0 ? n_miss0 : 0] z_miss0_raw; // historical data
}

transformed parameters {
  // Build complete design matrices with imputed values
  matrix[n, p] X_complete = X;
  matrix[n0, p] X0_complete = X0;

  // Impute missing values in current data
  if (n_miss > 0) {
    for (m in 1:n_miss) {
      int i = miss_idx[m, 1];       // row index
      int j = miss_idx[m, 2];       // design matrix column
      int cov_j = miss_idx[m, 3];   // covariate index in q

      real x_imp = z_to_x(z_miss_raw[m], marginal_type[cov_j],
                          marginal_mean[cov_j], marginal_sd[cov_j],
                          cutpoints[cov_j]', n_levels[cov_j]);
      X_complete[i, j] = x_imp;
    }
  }

  // Impute missing values in historical data
  if (n_miss0 > 0) {
    for (m in 1:n_miss0) {
      int i = miss_idx0[m, 1];      // row index
      int j = miss_idx0[m, 2];      // design matrix column
      int cov_j = miss_idx0[m, 3];  // covariate index in q

      real x_imp = z_to_x(z_miss0_raw[m], marginal_type[cov_j],
                          marginal_mean[cov_j], marginal_sd[cov_j],
                          cutpoints[cov_j]', n_levels[cov_j]);
      X0_complete[i, j] = x_imp;
    }
  }
}

model {
  vector[K] logprobs = append_row(logit_gamma, log(delta_raw)) - log1p_exp(logit_gamma);

  // Prior on latent missing values (standard normal - copula structure)
  if (n_miss > 0)
    z_miss_raw ~ std_normal();
  if (n_miss0 > 0)
    z_miss0_raw ~ std_normal();

  // Initial priors on regression coefficients
  for (k in 1:K) {
    target += normal_lpdf(betaMat[, k] | mean_beta[,k], sd_beta[,k]);
  }

  if (dist <= 2) {
    // Historical data likelihood (mixture)
    target += glm_mixture_lp(y0, betaMat, ones_vector(K), logprobs, X0_complete, dist, link, offs0);
    // Current data likelihood
    target += glm_lp(y, betaMat[, 1], 1.0, X_complete, dist, link, offs[,1]);
  }
  else {
    // Historical data likelihood (mixture)
    target += glm_mixture_lp(y0, betaMat, dispersion, logprobs, X0_complete, dist, link, offs0);
    // Half-normal prior for dispersion
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;
    // Current data likelihood
    target += glm_lp(y, betaMat[, 1], dispersion[1], X_complete, dist, link, offs[,1]);
  }

  // Beta prior on gamma (probability of exchangeability)
  target += logit_beta_lpdf(logit_gamma | gamma_shape1, gamma_shape2) - lognc_logit_gamma;

  // Dirichlet prior on delta (if K > 2)
  if (K_gt_2)
    target += dirichlet_lpdf(delta_raw | conc_delta);
}

generated quantities {
  vector[p] beta = betaMat[, 1];
  real<lower=gamma_lower,upper=gamma_upper> gamma = inv_logit(logit_gamma);
  simplex[K] probs;
  vector[(get_loglik == 1) ? n : 0] log_lik;

  // Compute mixture probabilities
  probs[1]   = gamma;
  probs[2:K] = (1 - gamma) * delta_raw;

  // Log-likelihood for model comparison
  if (get_loglik == 1) {
    real phi_val = (dist <= 2) ? 1.0 : dispersion[1];
    for (i in 1:n) {
      log_lik[i] = glm_lpdf(y[i] | beta, phi_val, X_complete[i, ], dist, link, offs[i, 1]);
    }
  }
}
