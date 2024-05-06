#' Compute mean from linear predictor in a GLM
#' @param eta  linear predictor
#' @param link integer giving link function
#' @noRd
get_lp2mean = function(eta, link) {
  if (link == 1){
    return(eta) # identity link
  }else if(link == 2){
    return(exp(eta)) # log link
  }else if(link == 3){
    return(binomial('logit')$linkinv(eta)) # logit link
  }else if(link == 4){
    return(1/eta) # inverse link
  }else if(link == 5){
    return(binomial('probit')$linkinv(eta)) # probit link
  }else if(link == 6){
    return(binomial('cauchit')$linkinv(eta)) # cauchit link
  }else if(link == 7){
    return(binomial('cloglog')$linkinv(eta)) # complementary log-log link
  }else if(link == 8){
    return(eta^2) # sqrt link
  }else if(link == 9){
    return(1/sqrt(eta)) # 1/mu^2 link
  }else{
    stop("Link not supported")
  }
}

#' compute density for normal GLM
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter (variance)
#' @noRd
normal_glm_lp = function(y, beta, X, link, offs, phi) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 1 ){
    theta = get_lp2mean(theta, link)
  }
  return( sum( dnorm(y, mean = theta, sd = sqrt(phi), log = T) ) )
}

#' compute density for bernoulli GLM
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter (phi = 1)
#' @noRd
bernoulli_glm_lp = function(y, beta, X, link, offs, phi = 1) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 3 ){
    theta = binomial('logit')$linkfun( get_lp2mean(theta, link) )
  }
  return( sum( y*theta - log1p(exp(theta)) ) )
}

#' compute density for poisson GLM
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter (phi = 1)
#' @noRd
poisson_glm_lp = function(y, beta, X, link, offs, phi = 1) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 2 ){
    theta = log( get_lp2mean(theta, link) )
  }
  return( sum( y*theta - exp(theta) - lgamma(y + 1) ) )
}

#' compute density for gamma GLM
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter
#' @noRd
gamma_glm_lp = function(y, beta, X, link, offs, phi) {
  n     = length(y)
  tau   = 1 / phi # shape parameter
  theta = X %*% beta + offs
  if ( link != 4 ){
    theta = 1 / get_lp2mean(theta, link)
  }
  return( sum( dgamma(y, shape = tau, rate = tau * theta, log = T) ) )
}

#' compute density for inverse-gaussian GLM
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter
#' @noRd
invgauss_glm_lp = function(y, beta, X, link, offs, phi) {
  n       = length(y)
  tau     = 1 / phi # shape parameter
  log_2pi = log(2*pi)
  theta   = X %*% beta + offs
  if ( link != 9 ){
    theta = 1 / ( get_lp2mean(theta, link)^2 )
  }
  return(
    0.5 * (
      n * (log(tau) - log_2pi) - 3 * sum(log(y)) - tau * sum( ( (y*sqrt(theta) - 1) * 1/sqrt(y) )^2 )
    )
  )
}

#' wrapper function compute density for a given link function and a given distribution
#' @param y    response vector
#' @param beta regression coefficients
#' @param X    design matrix
#' @param dist integer giving distribution
#' @param link integer giving link function
#' @param offs offset
#' @param phi  dispersion parameter
#' @noRd
glm_lp = function(y, beta, X, dist, link, offs, phi) {
  # Compute likelihood
  if (dist == 1) {     # Bernoulli
    return( bernoulli_glm_lp(y, beta, X, link, offs, phi) )
  }else if (dist == 2) {  # Poisson
    return( poisson_glm_lp(y, beta, X, link, offs, phi) )
  }else if (dist == 3) {  # Normal
    return( normal_glm_lp(y, beta, X, link, offs, phi) )
  }else if (dist == 4) { # Gamma
    return( gamma_glm_lp(y, beta, X, link, offs, phi) )
  }else if (dist == 5) { # Inverse-Gaussian
    return( invgauss_glm_lp(y, beta, X, link, offs, phi) )
  }else{
    stop("Distribution not supported")
  }
}
