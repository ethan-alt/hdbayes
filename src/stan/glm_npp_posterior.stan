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

  // find index of x, j, such that x0 is closest to x[j] without
  // going over. Uses binary search algorithm
  int findClosestIndex(real x0, vector x) {
    int K = rows(x);
    int i = 1;
    int j = K;
    int mid;
    // check corner cases
    if ( x0 < x[2] )
      return 1;
    if ( x0 == x[K] ) {
      return K;
    }
    // conduct binary search
    while ( i <= j ) {
      mid = (i + j) %/% 2;
      // if x0 < x[mid], index must lie in left half
      if ( x0 < x[mid] ) {
        // if x0 is larger than x[mid-1], return mid-1; else update j
        if ( mid > 2 &&  x0 > x[mid - 1] )
          return mid - 1;
        j = mid;
      }
      // otherwise, index must lie in right half
      else {
        // if x0 is less than x[mid + 1], return mid; else update i
        if ( mid < K && x0 < x[mid + 1] )
          return mid;
        i = mid + 1;
      }
    }
    reject("Error in finding midpoint");
    return(0); // never reached
  }

  // approximate lognc of power prior
  //
  // * @param a0       power prior param to obtain lognc
  // * @param a0vec    fine grid of power prior parameters for which we have estimates
  // * @param lognca0  estimate lognc pertaining to fine grid a0vec
  //
  // * @return linearly interpolated log normalizing constant.
  real pp_lognc(real a0, vector a0vec, vector lognca0) {
    // find index of a0vec closest to a0
    int i = findClosestIndex(a0, a0vec);
    // if not exact match, use linear interpolation to get estimated lognc
    if ( a0 != a0vec[i] ) {
      real x1 = a0vec[i];
      real x2 = a0vec[i + 1];
      real y1 = lognca0[i];
      real y2 = lognca0[i + 1];
      return y1 + (y2 - y1) * (a0 - x1) / (x2 - x1);
    }
    return lognca0[i];
  }

  real logit_beta_lpdf(real x, real shape1, real shape2) {
    return
      -lbeta(shape1, shape2) - shape2 * x - (shape1 + shape2) * log1p_exp(-x);
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
  vector[p]                           mean_beta; // mean for normal initial prior on coefficients
  vector<lower=0>[p]                  sd_beta; //sd for normal initial prior on coefficients
  real                                disp_mean; // mean for the half-normal prior on dispersion parameter
  real<lower=0>                       disp_sd; // sd for the half-normal prior on dispersion parameter
  int<lower=0>                        s; // number of a0s for which we have log nc
  vector<lower=0,upper=1>[s]          a0_lognc;
  matrix[s,K-1]                       lognc; // the j-th column is the log nc for a0_lognc using the j-th historical datasets
  real<lower=0>                       a0_shape1;
  real<lower=0>                       a0_shape2;
  vector<lower=0,upper=1>[K-1]        a0_lower; // lower bounds for a0s
  vector<lower=a0_lower,upper=1>[K-1] a0_upper; // upper bounds for a0s
  int<lower=1,upper=5>                dist;
  int<lower=1,upper=9>                link;
  vector[N]                           offs; // offset
}
transformed data{
  real lognc_disp      = normal_lccdf(0 | disp_mean, disp_sd);
  real lognc_logit_a0s = 0;

  for ( i in 1:(K-1) ) {
    if( a0_upper[i] != 1 || a0_lower[i] != 0 ) {
      lognc_logit_a0s = lognc_logit_a0s + log_diff_exp( beta_lcdf(a0_upper[i] | a0_shape1, a0_shape2), beta_lcdf(a0_lower[i] | a0_shape1, a0_shape2) );
    }
  }
}
parameters {
  vector[p] beta;
  vector<lower=0>[(dist > 2) ? 1 :  0] dispersion;
  vector<lower=logit(a0_lower),upper=logit(a0_upper)>[K-1] logit_a0s;
}
transformed parameters {
  vector<lower=a0_lower,upper=a0_upper>[K-1] a0s;
  a0s = inv_logit(logit_a0s);
}
model {
  // prior on beta
  target += normal_lpdf(beta  | mean_beta, sd_beta);
  if ( dist <= 2 ) {
    target += glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, 1.0, X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]); // current data likelihood

    for ( k in 2:K ) {
      // prior on logit(a0)
      target += logit_beta_lpdf(logit_a0s[k-1] | a0_shape1, a0_shape2);

      target += a0s[k-1] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, 1.0, X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]); // power prior
    }
  }
  else {
    target += normal_lpdf(dispersion | disp_mean, disp_sd) - lognc_disp;  // half-normal prior for dispersion
    target += glm_lp(y[ start_idx[1]:end_idx[1] ],
    beta, dispersion[1], X[ start_idx[1]:end_idx[1], ], dist, link,
    offs[ start_idx[1]:end_idx[1] ]); // current data likelihood

    for ( k in 2:K ) {
      // prior on logit(a0)
      target += logit_beta_lpdf(logit_a0s[k-1] | a0_shape1, a0_shape2);

      target += a0s[k-1] * glm_lp(y[ start_idx[k]:end_idx[k] ],
      beta, dispersion[1], X[ start_idx[k]:end_idx[k], ], dist, link,
      offs[ start_idx[k]:end_idx[k] ]);  // power prior
    }
  }
  target += -lognc_logit_a0s;

  // Subtract log nc from power prior
  for ( k in 2:K ) {
      target += -pp_lognc(a0s[k-1], a0_lognc, lognc[, k-1]);
  }
}

