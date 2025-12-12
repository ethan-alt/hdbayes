// GLM with Power Prior and Gaussian copula for missing covariates
functions{
  vector lp2mean(vector eta, int link) {
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
  real lp2mean(real eta, int link) {
    if (link == 1) return eta;
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

  real normal_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n = rows(y);
    vector[n] theta = X * beta + offs;
    if (link != 1) theta = lp2mean(theta, link);
    return normal_lpdf(y | theta, sqrt(phi));
  }
  real normal_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if (link != 1) theta = lp2mean(theta, link);
    return normal_lpdf(y | theta, sqrt(phi));
  }

  real bernoulli_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n = rows(y);
    vector[n] theta = X * beta + offs;
    if (link != 3) theta = logit(lp2mean(theta, link));
    return dot_product(y, theta) - sum(log1p_exp(theta));
  }
  real bernoulli_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if (link != 3) theta = logit(lp2mean(theta, link));
    return y * theta - log1p_exp(theta);
  }

  real poisson_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n = rows(y);
    vector[n] theta = X * beta + offs;
    if (link != 2) theta = log(lp2mean(theta, link));
    return dot_product(y, theta) - sum(exp(theta) + lgamma(y + 1));
  }
  real poisson_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    if (link != 2) theta = log(lp2mean(theta, link));
    return y * theta - exp(theta) - lgamma(y + 1);
  }

  real gamma_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n = rows(y);
    real tau = inv(phi);
    vector[n] theta = X * beta + offs;
    if (link != 4) theta = inv(lp2mean(theta, link));
    return gamma_lpdf(y | tau, tau * theta);
  }
  real gamma_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    real tau = inv(phi);
    if (link != 4) theta = inv(lp2mean(theta, link));
    return gamma_lpdf(y | tau, tau * theta);
  }

  real invgauss_glm_lp(vector y, vector beta, real phi, matrix X, int link, vector offs) {
    int n = rows(y);
    real tau = inv(phi);
    real log_2pi = 1.837877066409345483560659;
    vector[n] theta = X * beta + offs;
    if (link != 9) theta = inv_square(lp2mean(theta, link));
    return 0.5 * (n * (log(tau) - log_2pi) - 3 * sum(log(y))
                  - tau * dot_self((y .* sqrt(theta) - 1) .* inv_sqrt(y)));
  }
  real invgauss_glm_lpdf(real y, vector beta, real phi, row_vector x, int link, real offs) {
    real theta = dot_product(x, beta) + offs;
    real tau = inv(phi);
    real log_2pi = 1.837877066409345483560659;
    if (link != 9) theta = inv_square(lp2mean(theta, link));
    return 0.5 * (log(tau) - log_2pi - 3 * log(y)
                  - tau * square((y * sqrt(theta) - 1)) * inv(y));
  }

  real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offs) {
    if (dist == 1) return bernoulli_glm_lp(y, beta, phi, X, link, offs);
    else if (dist == 2) return poisson_glm_lp(y, beta, phi, X, link, offs);
    else if (dist == 3) return normal_glm_lp(y, beta, phi, X, link, offs);
    else if (dist == 4) return gamma_glm_lp(y, beta, phi, X, link, offs);
    else if (dist == 5) return invgauss_glm_lp(y, beta, phi, X, link, offs);
    else reject("Distribution not supported");
    return 0;
  }
  real glm_lpdf(real y, vector beta, real phi, row_vector x, int dist, int link, real offs) {
    if (dist == 1) return bernoulli_glm_lpdf(y | beta, phi, x, link, offs);
    else if (dist == 2) return poisson_glm_lpdf(y | beta, phi, x, link, offs);
    else if (dist == 3) return normal_glm_lpdf(y | beta, phi, x, link, offs);
    else if (dist == 4) return gamma_glm_lpdf(y | beta, phi, x, link, offs);
    else if (dist == 5) return invgauss_glm_lpdf(y | beta, phi, x, link, offs);
    else reject("Distribution not supported");
    return 0;
  }

  real z_to_x(real z, int type, real mean_param, real sd_param,
              vector cutpoints, int n_levels) {
    if (type == 1) {
      return mean_param + sd_param * z;
    }
    else if (type == 2) {
      return Phi_approx(z);
    }
    else if (type == 3) {
      return mean_param + sd_param * z;
    }
    else if (type == 4) {
      real lambda = mean_param;
      return fmax(0.0, lambda + sqrt(lambda) * z);
    }
    else {
      reject("Unknown marginal type");
      return 0.0;
    }
  }
}

data {
  int<lower=0>                        K;
  int<lower=0>                        N;
  array[K] int<lower=0, upper=N>      start_idx;
  array[K] int<lower=0, upper=N>      end_idx;
  int<lower=0>                        p;
  vector[N]                           y;
  matrix[N,p]                         X;
  vector[p]                           mean_beta;
  vector<lower=0>[p]                  sd_beta;
  vector[K]                           a0_vals;
  real                                disp_mean;
  real<lower=0>                       disp_sd;
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs;
  int<lower=0,upper=1>                get_loglik;

  // Missing data structure
  int<lower=0>                        n_miss;
  int<lower=0>                        q;
  array[n_miss, 3] int<lower=0>       miss_idx;
  matrix[q, q]                        copula_corr;
  array[q] int<lower=1,upper=4>       marginal_type;
  vector[q]                           marginal_mean;
  vector<lower=0>[q]                  marginal_sd;
  array[q] int<lower=0>               n_levels;
  int<lower=1>                        max_levels;
  matrix[q, max_levels-1]             cutpoints;
}

transformed data {
  real lognc_disp = normal_lccdf(0 | disp_mean, disp_sd);
  int n = end_idx[1];
}

parameters {
  vector[p]                           beta;
  vector<lower=0>[(dist > 2) ? 1 : 0] dispersion;
  vector[n_miss]                      z_miss_raw;
}

transformed parameters {
  matrix[N, p] X_complete = X;

  if (n_miss > 0) {
    for (m in 1:n_miss) {
      int row_idx = miss_idx[m, 1];
      int col_idx = miss_idx[m, 2];
      int cov_idx = miss_idx[m, 3];

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
  target += normal_lpdf(beta | mean_beta, sd_beta);
  target += std_normal_lpdf(z_miss_raw);

  if (dist <= 2) {
    for (k in 1:K) {
      target += a0_vals[k] * glm_lp(
        y[start_idx[k]:end_idx[k]],
        beta, 1.0,
        X_complete[start_idx[k]:end_idx[k], ],
        dist, link, offs[start_idx[k]:end_idx[k]]
      );
    }
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;
    for (k in 1:K) {
      target += a0_vals[k] * glm_lp(
        y[start_idx[k]:end_idx[k]],
        beta, dispersion[1],
        X_complete[start_idx[k]:end_idx[k], ],
        dist, link, offs[start_idx[k]:end_idx[k]]
      );
    }
  }
}

generated quantities {
  vector[(get_loglik == 1) ? n : 0] log_lik;

  if (get_loglik == 1) {
    real phi_val = (dist <= 2) ? 1.0 : dispersion[1];
    for (i in 1:n) {
      log_lik[i] = glm_lpdf(y[i] | beta, phi_val, X_complete[i, ], dist, link, offs[i]);
    }
  }
}
