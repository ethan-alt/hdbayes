functions{
#include include/expfam_loglik.stan
}
data {
  int<lower=0>          n0;
  int<lower=0>          p;
  vector[n0]            y0;
  matrix[n0,p]          X0;
  vector[p]             beta_mean;
  vector[p]             beta_sd;
  real                  disp_mean;
  real<lower=0>         disp_sd;
  int<lower=1,upper=5>  dist;
  int<lower=1,upper=9>  link;
  vector[n0]            offs0;
  real<lower=0,upper=1> a0;
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += normal_lpdf(beta|beta_mean, beta_sd);
  if ( dist <= 2 ) {
    target += a0 * glm_lp(y0, beta, 1.0, X0, dist, link, offs0);   // historical data likelihood
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - normal_lccdf(0 | disp_mean, disp_sd); // prior for dispersion
    target += a0 * glm_lp(y0, beta, dispersion[1], X0, dist, link, offs0); // historical data likelihood
  }
}

