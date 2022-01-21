real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offset) {
  int n    = rows(y);
  real tau = inv(phi);
  vector[n] eta = X * beta + offset;
  vector[n] mu;
  if      (link == 1) mu = eta;                      // identity link
  else if (link == 2) mu = exp(eta);                 // log link
  else if (link == 3) mu = inv_logit(eta);           // logit link
  else if (link == 4) mu = inv(eta);                 // inverse link
  else if (link == 5) mu = Phi_approx(eta);          // probit link
  else if (link == 6) mu = atan(eta) / pi() + 0.5;   // cauchit link
  else if (link == 7) mu = inv_cloglog(eta);         // complementary log-log link
  else if (link == 8) mu = square(eta);              // sqrt link
  else if (link == 9) mu = inv(sqrt(eta));           // 1/mu^2 link
  else reject("Link not supported");
  // Compute likelihood
  if (dist == 1) {     // Bernoulli
    return dot_product(y, log(mu)) + dot_product(1 - y, log(1 - mu));
  }
  else if (dist == 2) {  // Poisson
    return dot_product(y, log(mu)) - sum(mu + lgamma(y + 1));
  }
  else if (dist == 3) {  // Normal
    return normal_lpdf(y | mu, sqrt(phi) );
  }
  else if (dist == 4) { // Gamma
    return gamma_lpdf(y | tau, tau * mu );
  }
  else if (dist == 5) { // Inverse-Gaussian
    return 
        0.5 * (
            n * (log(tau) - 1.8378770664) - 3 * sum(log(y))
          - tau * dot_self( (y .* inv(mu) - 1) .* inv_sqrt(y) )
        );
  }
  else reject("Distribution not supported");
  return 0; // never reached;
}
