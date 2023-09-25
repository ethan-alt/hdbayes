data {
  int<lower=0>                        K; // total number of historical datasets (excluding the current data)
  int<lower=0>                        n; // number of observations in the current data
  vector[K]                           n0s; // sample size for each historical dataset
  int<lower=0>                        p;
  vector[n]                           y; // response for the current data after adjusting for offsets
  matrix[n,p]                         X; // design matrix for the current data
  array[K] matrix[p,p]                hist_prec; // each element is X0'X0 for a historical dataset
  matrix[p,K]                         hist_mle; // each column is betahat0 for a historical dataset
  vector[K]                           sumy0sq; // each element is sum(y0^2) for a historical dataset
  vector[p]                           mean_beta; // initital prior mean for beta
  cov_matrix[p]                       cov_beta;  // covariance for initial prior for beta is cov_beta * sigmasq
  real<lower=0>                       sigmasq_shape;  // shape parameter for inv gamma prior for sigmasq
  real<lower=0>                       sigmasq_scale;  // scale parameter for inv gamma prior for sigmasq
  real<lower=0>                       a0_shape1;      // shape parameter for a0s
  real<lower=0>                       a0_shape2;      // shape2 parameter for a0s
}
transformed data{
  matrix[p,K]              hist_prec_mle; // each column is X0'X0 * betahat0 for a historical dataset
  matrix[p,p]              prior_prec      = inverse_spd(cov_beta); // cov_beta^(-1) = Omega0
  vector[p]                prior_prec_mean = prior_prec * mean_beta; // Omega0 * mean_beta
  real                     quad_b0         = quad_form(prior_prec, mean_beta);  // mean_beta' Omega0 mean_beta

  for ( k in 1:K ) {
    hist_prec_mle[, k] = hist_prec[k] * hist_mle[, k];  // X0'X0 * betahat0
  }
}
parameters {
  vector[p] beta;
  real<lower=0> sigmasq;
  vector<lower=0,upper=1>[K] a0s;
}
model {
  // Compute parameters for implied normal-inverse-Gamma prior on beta, sigma^2
  vector[p]   mean_a0; // mu* = (sum of a0 X0'X0 + Omega0)^(-1) (sum of a_0 X0'X0 betahat0 + Omega0 mean_beta) = mean | sigmasq, a0
  real        scale_a0;
  matrix[p,p] prec_a0          = prior_prec; // Omega* =  Omega0 + (sum of a0 X0'X0 over all historical data) = sigmasq * precision matrix | sigmasq, a0
  vector[p]   a0_hist_prec_mle = rep_vector(0.0, p); // sum of a_0 X0'X0 betahat0 over all historical data

  // shape for sigmasq = delta0 + 0.5 * (sum of a0 * n0 over all historical data)
  real        shape_a0         = sigmasq_shape + 0.5 * dot_product(a0s, n0s);
  // scale for sigmasq = gamma0 + 0.5 * [ (sum of a0 * y0' y0 over all historical data) + mean_beta' Omega0 mean_beta - mu*' Omega* mu*]

  for ( i in 1:K ) {
    prec_a0          += a0s[i] * hist_prec[i];
    a0_hist_prec_mle += a0s[i] * hist_prec_mle[, i];
  }

  mean_a0   = mdivide_left_spd(prec_a0, a0_hist_prec_mle + prior_prec_mean);
  scale_a0  = sigmasq_scale + 0.5 * ( dot_product(a0s, sumy0sq) + quad_b0 - quad_form(prec_a0, mean_a0) );

  // iid beta prior for a0s
  if ( a0_shape1 != 1 || a0_shape2 != 1)
    a0s ~ beta(a0_shape1, a0_shape2);

  // power prior is normal-inverse Gamma
  target += inv_gamma_lpdf(sigmasq | shape_a0, scale_a0);
  target += multi_normal_prec_lpdf(beta | mean_a0, prec_a0 * inv(sigmasq) );

  // likelihood of current data
  target += normal_id_glm_lpdf(y | X, 0.0, beta, sqrt(sigmasq));
}

