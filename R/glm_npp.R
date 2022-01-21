
#'
#' Posterior of normalized power prior
#'
#' Sample from the posterior distribution of a GLM using the normalized power prior
#' (NPP). Before using this function, users must estimate the logarithm of the
#' normalizing constant across a range of power prior parameters (a0), possibly
#' smoothing techniques over a find grid.
#'
#'
#' @include data_checks.R
#' @include glm_npp_lognc.R
#'
#' @export
#'
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data        a `data.frame` giving the current data
#' @param histdata    a `data.frame` giving the historical data
#' @param beta.mean   mean parameter for initial prior on regression coefficients (including intercept). Defaults to a vector of zeros.
#' @param beta.cov    covariance parameter for initial prior on regression coefficients (including intercept). Defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param disp.shape  shape parameter for inverse-gamma prior on dispersion parameter
#' @param disp.scale  scale parameter for inverse-gamma prior on dispersion parameter
#' @param a0          vector giving values of the power prior parameter for which the logarithm of the normalizing constant has been evaluated
#' @param lognc       vector (of same length as a0) giving the logarithm of the normalizing constant (as estimated by \code{\link[hdbayes]{glm.npp.lognc}})
#' @param a0.shape1   first shape parameter for beta prior on a0. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2   first shape parameter for beta prior on a0. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param offset      vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0     vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...         arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#'
glm.npp = function(
  formula,
  family,
  data,
  histdata,
  a0,
  lognc,
  beta.mean   = NULL,
  beta.cov    = NULL,
  disp.shape  = 2.1,
  disp.scale  = 1.1,
  a0.shape1   = 1,
  a0.shape2   = 1,
  offset      = NULL,
  offset0     = NULL,
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
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

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

  ## check a0 and lognc
  if (length(a0) != length(lognc))
    stop('a0 and lognc must have the same length')
  if ( any(is.na(a0) ) )
    stop('a0 must not have missing values')
  if ( any(is.na(lognc)) )
    stop('lognc must not have missing values')
  if ( any( a0 < 0 || a0 > 1 ) )
    stop('each element of a0 should be between 0 and 1')

  ## perform data checks
  data.checks(
    formula, family, data, histdata, offset, offset0, check.hist = TRUE
  )

  standat = list(
    'n'           = n,
    'n0'          = n0,
    'p'           = p,
    'k'           = length(a0),
    'y'           = y,
    'X'           = X,
    'y0'          = y0,
    'X0'          = X0,
    'beta_mean'   = beta.mean,
    'beta_cov'    = beta.cov,
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'a0_vec'      = a0,
    'lognca0_vec' = lognc,
    'a0_shape1'   = a0.shape1,
    'a0_shape2'   = a0.shape2,
    'dist'        = dist,
    'link'        = link,
    'offset'      = offset,
    'offset0'     = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$glm_npp_posterior, data = standat, ...)

  ## rename parameters
  newnames = colnames(X)
  fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}





