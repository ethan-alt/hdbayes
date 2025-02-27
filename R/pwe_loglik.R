#' compute log density for a PWE model
#'
#' @param y         observed event times
#' @param eta       linear predictor for each obs
#' @param lambda    baseline hazards
#' @param breaks    (J+1)-dim vector giving intervals
#' @param j         index giving interval into which obs i failed / was censored
#' @param J         number of time intervals
#' @param death_ind event indicator (1 = event; 0 = censored)
#' @noRd
pwe_lpdf = function(y, eta, lambda, breaks, j, J, death_ind) {
  # Get hazard corresponding to failure / censoring time
  lambda_j = lambda[j]

  # Compute cumulative baseline hazard at each interval
  cumblhaz = cumsum( lambda[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) )
  cumblhaz = c(0, cumblhaz)

  # Compute cumulative hazard for each observation
  cumhaz = lambda_j * (y - breaks[j]) + cumblhaz[j]
  cumhaz = cumhaz * exp(eta)

  # log likelihood = event_ind * log(hazard) - cumhaz
  # log(hazard) = log( lambda * exp(eta) ) = log(lambda) + eta
  loglik = death_ind * ( log(lambda_j) + eta ) - cumhaz
  return(loglik)
}
