#' Helper functions for computing log likelihood for a mixture of generalized linear models

#' compute log sum of exponentials
#' @param x vector
#' @noRd
log_sum_exp = function(x){
  x_max = max(x)
  x_max + log( sum(exp(x - x_max)) )
}

#' compute log(1 + exp(x))
#'
#' @param x real value
#' @noRd
log1p_exp = function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' compute log density for dirichlet
#' @param x real number
#' @param conc first shape parameter for Beta distribution
#' @noRd
dirichlet_lp = function(x, conc){
  lgamma( sum(conc) ) - sum( lgamma(conc) ) + as.numeric( (conc - 1) %*% log(x) )
}

#' compute the log likelihood contribution of each individual to each component of a mixture of normal GLMs
#'
#' @include expfam_loglik.R
#'
#' @param y     response vector
#' @param beta  matrix of regression coefficients
#' @param X     design matrix
#' @param link  integer giving link function
#' @param offs  offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
normal_glm_mixture_contrib = function(y, X, beta, disp, log_probs, link, offs) {
  log_2pi = 1.837877066409345483560659
  K         = length(log_probs)
  theta     = X %*% beta + offs
  if ( link != 1 ){
    theta = get_lp2mean(theta, link)
  }

  inv_disp  = 1/disp
  y_sq      = y^2
  contrib   = sapply(1:K, function(k){
    log_probs[k] +
      inv_disp[k] * (y * theta[, k] - 0.5 * theta[, k]^2) -
      0.5 * ( y_sq * inv_disp[k] + log_2pi + log(disp[k]) )
  })
  return( contrib )
}

#' compute the log likelihood contribution of each individual to each component of a mixture of Bernoulli GLMs
#'
#' @include expfam_loglik.R
#'
#' @param y     response vector
#' @param beta  matrix of regression coefficients
#' @param X     design matrix
#' @param link  integer giving link function
#' @param offs  offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
bernoulli_glm_mixture_contrib = function(y, X, beta, disp, log_probs, link, offs) {
  K         = length(log_probs)
  theta     = X %*% beta + offs
  if ( link != 3 ){
    theta = binomial('logit')$linkfun( get_lp2mean(theta, link) )
  }

  contrib   = sapply(1:K, function(k){
    log_probs[k] + y * theta[, k] - log1p_exp( theta[, k] )
  })
  return( contrib )
}

#' compute the log likelihood contribution of each individual to each component of a mixture of Poisson GLMs
#'
#' @include expfam_loglik.R
#'
#' @param y     response vector
#' @param beta  matrix of regression coefficients
#' @param X     design matrix
#' @param link  integer giving link function
#' @param offs  offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
poisson_glm_mixture_contrib = function(y, X, beta, disp, log_probs, link, offs) {
  K         = length(log_probs)
  theta     = X %*% beta + offs
  if ( link != 2 ){
    theta = log( get_lp2mean(theta, link) )
  }

  log_y_factorial = lgamma(y+1)
  contrib         = sapply(1:K, function(k){
    log_probs[k] + y * theta[, k] - exp(theta[, k]) - log_y_factorial
  })
  return( contrib )
}

#' compute the log likelihood contribution of each individual to each component of a mixture of Gamma GLMs
#'
#' @include expfam_loglik.R
#'
#' @param y     response vector
#' @param beta  matrix of regression coefficients
#' @param X     design matrix
#' @param link  integer giving link function
#' @param offs  offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
gamma_glm_mixture_contrib = function(y, X, beta, disp, log_probs, link, offs) {
  K         = length(log_probs)
  theta     = X %*% beta + offs
  if ( link != 4 ){
    theta = 1 / get_lp2mean(theta, link)
  }

  log_y     = log(y)
  inv_disp  = 1/disp
  contrib   = sapply(1:K, function(k){
    log_probs[k] +
      inv_disp[k] * ( -y * theta[, k] + log(theta[, k]) ) +
      (inv_disp[k] - 1) * log_y - inv_disp[k] * log(disp[k]) - lgamma(inv_disp[k])
  })
  return( contrib )
}

#' compute the log likelihood contribution of each individual to each component of a mixture of inverse-Gaussian GLMs
#'
#' @include expfam_loglik.R
#'
#' @param y     response vector
#' @param beta  matrix of regression coefficients
#' @param X     design matrix
#' @param link  integer giving link function
#' @param offs  offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
invgauss_glm_mixture_contrib = function(y, X, beta, disp, log_probs, link, offs) {
  log_2pi   = 1.837877066409345483560659
  K         = length(log_probs)
  theta     = X %*% beta + offs
  if ( link != 9 ){
    theta = 1 / ( get_lp2mean(theta, link)^2 )
  }

  log_y_cubed = 3 * log(y)
  inv_y       = 1/y
  inv_disp    = 1/disp
  contrib     = sapply(1:K, function(k){
    log_probs[k] +
      + inv_disp[k] * ( -0.5 * y * theta[, k] + sqrt(theta[, k]) ) -
      0.5 * ( inv_disp[k] * inv_y + log(disp[k]) + log_2pi + log_y_cubed )
  })
  return( contrib )
}

#' wrapper function to compute a n x K matrix of log likelihood contributions for a mixture of GLMs
#'
#' @param y    response vector
#' @param beta matrix of regression coefficients
#' @param X    design matrix
#' @param dist integer giving distribution
#' @param link integer giving link function
#' @param offs offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
glm_mixture_contrib = function(y, beta, disp, log_probs, X, dist, link, offs) {
  if (dist == 1) {     # Bernoulli
    return( bernoulli_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs) )
  }else if (dist == 2) {  # Poisson
    return( poisson_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs) )
  }else if (dist == 3) {  # Normal
    return( normal_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs) )
  }else if (dist == 4) { # Gamma
    return( gamma_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs) )
  }else if (dist == 5) { # Inverse-Gaussian
    return( invgauss_glm_mixture_contrib(y, X, beta, disp, log_probs, link, offs) )
  }else{
    stop("Distribution not supported")
  }
}

#' compute the log likelihood for a mixture of GLMs
#'
#' @param y    response vector
#' @param beta matrix of regression coefficients
#' @param X    design matrix
#' @param dist integer giving distribution
#' @param link integer giving link function
#' @param offs offset
#' @param disp  vector of dispersion parameters
#' @param log_probs vector of log of component probabilities
#' @noRd
glm_mixture_lp = function(y, beta, disp, log_probs, X, dist, link, offs) {
  n        = length(y)
  contribs = glm_mixture_contrib(y, beta, disp, log_probs, X, dist, link, offs)
  contrib  = sapply(1:n, function(i){
    log_sum_exp(contribs[i,])
  })
  return( sum(contrib) )
}
