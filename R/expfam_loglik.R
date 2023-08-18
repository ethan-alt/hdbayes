# Compute mean from linear predictor in a GLM
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

normal_glm_lp = function(y, beta, X, link, offs, phi) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 1 ){
    theta = get_lp2mean(theta, link)
  }
  return( dnorm(y, mean = theta, sd = sqrt(phi), log = T) )
}

bernoulli_glm_lp = function(y, beta, X, link, offs, phi = 1) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 3 ){
    theta = binomial('logit')$linkfun( get_lp2mean(theta, link) )
  }
  return( sum( y*theta - log(1 + exp(theta)) ) ) # dot_product(y, theta) - sum( log1p_exp(theta) )
}

poisson_glm_lp = function(y, beta, X, link, offs, phi = 1) {
  n = length(y)
  theta = X %*% beta + offs
  if ( link != 2 ){
    theta = log( get_lp2mean(theta, link) )
  }
  return( sum( y*theta - exp(theta) - lgamma(y + 1) ) ) # dot_product(y, theta) - sum( exp(theta) + lgamma(y + 1) );
}

gamma_glm_lp = function(y, beta, X, link, offs, phi) {
  n     = length(y)
  tau   = 1 / phi # shape parameter
  theta = X %*% beta + offs
  if ( link != 4 ){
    theta = 1 / get_lp2mean(theta, link)
  }
  return( dgamma(y, shape = tau, rate = tau * theta, log = T) ) # gamma_lpdf(y | tau, tau * theta );
}

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
  # 0.5 * (n * (log(tau) - log_2pi) - 3 * sum(log(y)) - tau * dot_self( (y .* sqrt(theta) - 1) .* inv_sqrt(y) ) );
}

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
