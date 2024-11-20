#' compute log density for Gumbel distribution (same as gumbel_lpdf in Stan)
#' density function: f(x) = 1/sigma * exp{-z - exp(-z)}, where z = (x-mu)/sigma
#'
#' @param x     real number
#' @param mu    location parameter for Gumbel distribution
#' @param sigma scale parameter for Gumbel distribution
#' @noRd
log_dgumbel = function(x, mu = 0, sigma = 1) {
  if ( any(sigma <= 0) ){
    stop("Scale parameter sigma must be positive")
  }
  z = (x - mu) / sigma
  return( -log(sigma) - z - exp(-z) )
}

#' compute log(CDF) for Gumbel distribution (same as gumbel_lcdf in Stan)
#' CDF(x) = exp{-exp(-z)}, where z = (x-mu)/sigma
#'
#' @param x     real number
#' @param mu    location parameter for Gumbel distribution
#' @param sigma scale parameter for Gumbel distribution
#' @noRd
log_pgumbel <- function(x, mu = 0, sigma = 1) {
  if (sigma <= 0){
    stop("Scale parameter sigma must be positive")
  }
  z = (x - mu) / sigma
  return( -exp(-z) )
}

#' compute log density for an AFT model
#'
#' @param y_obs      log of observed event time (uncensored)
#' @param y_cen      log of censored time
#' @param eta_obs    linear predictor for uncensored
#' @param eta_cen    linear predictor for censored
#' @param scale      scale parameter of the AFT model
#' @param dist       integer giving distribution
#' @noRd
aft_model_lp = function(y_obs, y_cen, eta_obs, eta_cen, scale, dist) {
  # Compute likelihood
  if ( dist == 1 ) { # log-normal
    loglik = sum( stats::dnorm(y_obs, mean = eta_obs, sd = scale, log = T) ) + # uncensored data
      sum( stats::pnorm(y_cen, mean = eta_cen, sd = scale, lower.tail = F, log.p = T) ) # censored data
  }else if ( dist == 2 ) { # log-logistic
    loglik = sum( stats::dlogis(y_obs, location = eta_obs, scale = scale, log = T) ) + # uncensored data
      sum( stats::plogis(y_cen, location = eta_cen, scale = scale, lower.tail = F, log.p = T) ) # censored data
  }else if ( dist == 3 ) { # weibull
    loglik = sum( log_dgumbel(-y_obs, mu = -eta_obs, sigma = scale) ) + # uncensored data
      sum( log_pgumbel(-y_cen, mu = -eta_cen, sigma = scale) ) # censored data
  }else{
    stop("Distribution not supported.")
  }
  return(loglik)
}
