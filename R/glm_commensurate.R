
#'
#' Commensurate prior
#'
#' Sample from the posterior distribution of a GLM using the Commensurate
#' prior of Hobbs et al. This prior assumes that the regression coefficients
#' for the current data set conditional on those for the historical data set are
#' multivariate normal with mean equal to the regression coefficients of the historical data
#' and covariance equal to the inverse of a user-specified precision parameter (tau).
#' Dispersion parameters (if applicable) between the
#' current and historical data sets are treated as independent.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data        a `data.frame` giving the current data
#' @param histdata    a `data.frame` giving the historical data
#' @param tau         a vector giving the commensurate prior parameter. The
#'                    dimension must equal the number of predictors (including the intercept). Each
#'                    element must be positive, corresponding to a normal precision
#'                    parameter.
#' @param beta0.mean  mean parameter for initial prior on regression coefficients (including intercept). Defaults to a vector of zeros.
#' @param beta0.cov   covariance parameter for initial prior on regression coefficients (including intercept). Defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param disp.shape  shape parameter for inverse-gamma prior on dispersion parameter for current data set
#' @param disp.scale  scale parameter for inverse-gamma prior on dispersion parameter for current data set
#' @param disp0.shape shape parameter for inverse-gamma prior on dispersion parameter for historical data set
#' @param disp0.scale scale parameter for inverse-gamma prior on dispersion parameter for historical data set
#' @param offset      vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0     vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...         arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#'
glm.commensurate = function(
  formula,
  family,
  data,
  histdata,
  tau,
  beta0.mean  = NULL,
  beta0.cov   = NULL,
  disp.shape  = 2.1,
  disp.scale  = 1.1,
  disp0.shape = 2.1,
  disp0.scale = 1.1,
  offset  = NULL,
  offset0 = NULL,
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

  ## Check if each element of tau is positive
  if ( any(is.na(tau) ) )
    stop('tau must be a vector of non-missing values')
  if ( any(tau <= 0) )
    stop("Each element of tau must be positive")

  ## Check if tau is of proper length
  if ( length(tau) != ncol(X) )
    stop(
      paste0("The design matrix is of dimension ", ncol(X), " but tau is of dimension ", length(tau))
    )

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta0.mean) )
    beta0.mean = rep(0, ncol(X))
  if ( is.null(beta0.cov) )
    beta0.cov  = diag(100, ncol(X))

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
    'beta0_mean'  = beta0.mean,
    'beta0_cov'   = beta0.cov,
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'disp0_shape' = disp0.shape,
    'disp0_scale' = disp0.scale,
    'tau'         = tau,
    'dist'        = dist,
    'link'        = link,
    'offset'      = offset,
    'offset0'     = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$glm_commensurate, data = standat, ...)

  ## rename parameters
  if ( family$family %in% c('binomial', 'poisson') ) {
    newnames = c(colnames(X), paste0( colnames(X0), '_hist' ), 'lp__')
  } else {
    newnames = c(colnames(X), 'dispersion')
    newnames = c( newnames, paste0(newnames, '_hist'), 'lp__' )
  }
  fit@sim$fnames_oi = newnames
  return(fit)
}
