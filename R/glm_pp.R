
#'
#' Posterior of power prior with fixed a0
#'
#' Sample from the posterior distribution of a GLM using the power prior with
#' a fixed power prior parameter (a0)
#'
#'
#' @include data_checks.R
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
#' @param a0          scalar between 0 and 1 giving the (fixed) power prior parameter
#' @param offset      vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0     vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...         arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#'
glm.pp = function(
  formula,
  family,
  data,
  histdata,
  a0,
  beta.mean   = NULL,
  beta.cov    = NULL,
  disp.shape  = 2.1,
  disp.scale  = 1.1,
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

  ## check a0
  if ( length(a0) != 1)
    stop('a0 must be a scalar')
  if ( a0 < 0 | a0 > 1 )
    stop("a0 must be a scalar between 0 and 1")

  ## perform data checks
  data.checks(
    formula, family, data, histdata, offset, offset0, check.hist = TRUE
  )

  standat = list(
    'n'           = n,
    'n0'          = n0,
    'p'           = p,
    'y'           = y,
    'X'           = X,
    'y0'          = y0,
    'X0'          = X0,
    'mean_beta'   = beta.mean,
    'cov_beta'    = beta.cov,
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'a0'          = a0,
    'dist'        = dist,
    'link'        = link,
    'offset'      = offset,
    'offset0'     = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$glm_pp, data = standat, ...)

  ## rename parameters
  newnames = colnames(X)
  fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}





