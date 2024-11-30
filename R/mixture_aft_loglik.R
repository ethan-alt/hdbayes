#' Helper functions for computing log likelihood for a mixture of accelerated failure time (AFT) models

#' compute log-likelihood for an uncensored observation
#'
#' @include aft_loglik.R
#'
#' @param y_obs      log of observed event time (uncensored)
#' @param eta_obs    linear predictor for uncensored
#' @param scale      scale parameter of the AFT model
#' @param dist       integer giving distribution
#' @noRd
aft_model_obs_lpdf = function(y_obs, eta_obs, scale, dist) {
  # Compute likelihood
  if ( dist == 1 ) { # log-normal
    loglik = stats::dnorm(y_obs, mean = eta_obs, sd = scale, log = T) # uncensored data
  }else if ( dist == 2 ) { # log-logistic
    loglik = stats::dlogis(y_obs, location = eta_obs, scale = scale, log = T) # uncensored data
  }else if ( dist == 3 ) { # weibull
    loglik = log_dgumbel(-y_obs, mu = -eta_obs, sigma = scale) # uncensored data
  }else{
    stop("Distribution not supported.")
  }
  return(loglik)
}

#' compute log-likelihood for a censored observation
#'
#' @include aft_loglik.R
#'
#' @param y_cen      log of censored time
#' @param eta_cen    linear predictor for censored
#' @param scale      scale parameter of the AFT model
#' @param dist       integer giving distribution
#' @noRd
aft_model_cen_lpdf = function(y_cen, eta_cen, scale, dist) {
  # Compute likelihood
  if ( dist == 1 ) { # log-normal
    loglik = pnorm(y_cen, mean = eta_cen, sd = scale, lower.tail = F, log.p = T) # censored data
  }else if ( dist == 2 ) { # log-logistic
    loglik = stats::plogis(y_cen, location = eta_cen, scale = scale, lower.tail = F, log.p = T) # censored data
  }else if ( dist == 3 ) { # weibull
    loglik = log_pgumbel(-y_cen, mu = -eta_cen, sigma = scale) # censored data
  }else{
    stop("Distribution not supported.")
  }
  return(loglik)
}

#' compute log-likelihood of a mixture model for an uncensored observation
#'
#' @include aft_loglik.R
#' @include expfam_loglik.R
#'
#' @param y_obs      log of observed event time (uncensored)
#' @param eta_obs    linear predictor for uncensored
#' @param scale      a vector of scale parameters of the AFT models
#' @param log_probs  vector of log of component probabilities
#' @param dist       integer giving distribution
#' @noRd
aft_model_obs_mixture_lp = function(y_obs, eta_obs, scale, log_probs, dist) {
  K = length(log_probs)
  # log likelihood contribution of each component
  contribs = sapply(1:K, function(k){
    log_probs[k] + aft_model_obs_lpdf(y_obs, eta_obs[k], scale[k], dist)
  })

  return(log_sum_exp(contribs))
}

#' compute log-likelihood of a mixture model for a censored observation
#'
#' @include aft_loglik.R
#' @include expfam_loglik.R
#'
#' @param y_cen      log of censored time
#' @param eta_cen    linear predictor for censored
#' @param scale      a vector of scale parameters of the AFT models
#' @param log_probs  vector of log of component probabilities
#' @param dist       integer giving distribution
#' @noRd
aft_model_cen_mixture_lp = function(y_cen, eta_cen, scale, log_probs, dist) {
  K = length(log_probs)
  # log likelihood contribution of each component
  contribs = sapply(1:K, function(k){
    log_probs[k] + aft_model_cen_lpdf(y_cen, eta_cen[k], scale[k], dist)
  })

  return(log_sum_exp(contribs))
}
