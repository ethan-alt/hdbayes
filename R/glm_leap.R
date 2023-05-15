
#'
#' Posterior of LEAP
#'
#' ...
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
#' @param mean_beta   mean parameter for initial prior on regression coefficients (including intercept). Defaults to a vector of zeros.
#' @param cov_beta    covariance parameter for initial prior on regression coefficients (including intercept). Defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param disp_shape  shape parameter for inverse-gamma prior on dispersion parameter
#' @param disp_scale  scale parameter for inverse-gamma prior on dispersion parameter
#' @param offset      vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0     vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...         arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' ## take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' glm.leap(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(), data = actg019, histdata = actg036,
#'   chains = 1, warmup = 500, iter = 1000
#' )
#'
#'
#'
glm.leap = function(
    formula,
    family,
    data,
    histdata,
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

  ## perform data checks
  data.checks(
    formula, family, data, histdata, offset, offset0, check.hist = TRUE
  )

  standat = list(
    'n'           = n,
    'n0'          = n0,
    'p'           = p,
    'K'           = 2,
    'y'           = y,
    'X'           = X,
    'y0'          = y0,
    'X0'          = X0,
    'mean_beta'   = beta.mean,
    'cov_beta'    = beta.cov,
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'conc'        = c(0.95, 0.95),
    'gamma_lower' = 0,
    'gamma_upper' = 1,
    'dist'        = dist,
    'link'        = link,
    'offset'      = offset,
    'offset0'     = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$glm_leap, data = standat, ...)

  ## rename parameters
  newnames = colnames(X)
  fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}





