
#'
#' Posterior of normalized power prior
#'
#' Sample from the posterior distribution of a normal linear model
#' using the normalized power prior (NPP)
#'
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula        a two-sided formula giving the relationship between the response variable and covariates
#' @param data           a `data.frame` giving the current data
#' @param histdata       a `data.frame` giving the historical data
#' @param beta.mean      mean parameter for initial prior on regression coefficients (including intercept). Defaults to a vector of zeros.
#' @param beta.cov       covariance parameter for initial prior on regression coefficients (including intercept). Defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param sigmasq.shape  shape parameter for inverse-gamma prior on variance parameter
#' @param sigmasq.scale  scale parameter for inverse-gamma prior on variance parameter
#' @param a0.shape1      first shape parameter for beta prior on a0. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2      first shape parameter for beta prior on a0. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param offset         vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0        vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...            arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' lm.npp(
#'    cd4 ~ treatment + age + race,
#'    data = actg019, histdata = actg036, chains = 1
#' )
#'
#'
lm.npp = function(
  formula,
  data,
  histdata,
  beta.mean      = NULL,
  beta.cov       = NULL,
  sigmasq.shape  = 2.1,
  sigmasq.scale  = 1.1,
  a0.shape1      = 1,
  a0.shape2      = 1,
  offset         = NULL,
  offset0        = NULL,
  ...
) {
  ## get model information
  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]
  n  = length(y)
  n0 = length(y0)
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)
  p  = ncol(X)

  ## Default offset is vector of 0s
  if ( is.null(offset) )
    offset = rep(0, n)
  if ( is.null(offset0) )
    offset0 = rep(0, n0)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) )
    beta.mean = rep(0, ncol(X))
  if ( is.null(beta.cov) )
    beta.cov  = diag(100, ncol(X))

  ## perform data checks
  data.checks(
    formula, gaussian('identity'), data, histdata, offset, offset0, check.hist = TRUE
  )

  standat = list(
    'n'              = n,
    'n0'             = n0,
    'p'              = p,
    'y'              = y,
    'X'              = X,
    'y0'             = y0,
    'X0'             = X0,
    'mean_beta'      = beta.mean,
    'cov_beta'       = beta.cov,
    'sigmasq_shape'  = sigmasq.shape,
    'sigmasq_scale'  = sigmasq.scale,
    'a0_shape1'      = a0.shape1,
    'a0_shape2'      = a0.shape2,
    'offset'         = offset,
    'offset0'        = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$lm_npp, data = standat, ...)

  ## rename parameters
  newnames = colnames(X)
  fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}





