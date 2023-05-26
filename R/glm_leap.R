
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
#' @param data.list   a list of `data.frame` giving the current data followed by historical data
#' @param K           the desired number of classes to identify
#' @param beta.mean   a `p x K` matrix of mean parameter for initial prior on regression coefficients (including intercept). Defaults to a matrix of zeros.
#' @param beta.cov    a list of `K` `data.frame` each size `p x p` giving covariance parameter for initial prior on regression coefficients (including intercept). Each defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param disp.shape  shape parameter for inverse-gamma prior on dispersion parameter
#' @param disp.scale  scale parameter for inverse-gamma prior on dispersion parameter
#' @param offset.list a list of `data.frame` giving the offset for current data followed by historical data. For each `data.frame`, the number of rows correpond to observations and columns correspond to classes. Defaults to matrices of 0s
#' @param ...         arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return            an object of class `stanfit` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' # take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' glm.leap(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(),
#'   data.list = list(actg019, actg036),
#'   K = 2,
#'   chains = 1, warmup = 500, iter = 1000
#' )
#'
glm.leap = function(
    formula,
    family,
    data.list,
    K           = 2,
    beta.mean   = NULL,
    beta.cov    = NULL,
    disp.shape  = 2.1,
    disp.scale  = 1.1,
    offset.list = NULL,
    ...
) {
  ## get model information
  data=data.list[[1]]
  histdata=data.list[[2]] # not yet support multiple histdata; only read 2nd dataset
  offset=offset.list[[1]]
  offset0=offset.list[[2]]

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

  ## Default offset is matrix of 0s
  if ( is.null(offset) )
    offset = matrix(rep(0, n*K), ncol=K)
  if ( is.null(offset0) )
    offset0 = matrix(rep(0, n0*K), ncol=K)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) )
    beta.mean = matrix(rep(0, ncol(X)*K), ncol=K)
  if ( is.null(beta.cov) )
    beta.cov  = replicate(K, diag(100, ncol(X)), simplify=F)

  ## perform data checks
  data.checks(formula, family, data.list, offset.list=NULL)

  standat = list(
    'n'           = n,
    'n0'          = n0,
    'p'           = p,
    'K'           = K,
    'y'           = y,
    'X'           = X,
    'y0'          = y0,
    'X0'          = X0,
    'mean_beta'   = beta.mean,
    'cov_beta'    = do.call(rbind, beta.cov),
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





