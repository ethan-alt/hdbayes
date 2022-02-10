#' Bayesian hierarchical model
#'
#' Sample from the posterior distribution of a GLM using the Bayesian Hierarchical
#' Model (BHM). This model assumes that the regression coefficients of the
#' historical and current data are different, but are correlated through a
#' common distribution, whose hyperparameters (i.e., mean and covariance) are treated as random.
#' The mean hyperparmater is multivariate normal and the covariance hyperparameter is
#' inverse-Wishart.
#' The dispersion parameters (if applicable) are independent.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula      a two-sided formula giving the relationship between the response variable and covariates
#' @param family       an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data         a `data.frame` giving the current data
#' @param histdata     a `data.frame` giving the historical data
#' @param norm.hp.mean a vector whose dimension is equal to the number of regression coefficients giving the
#'                     mean for the normal hyperprior on the regression coefficients. Defaults to a vector of zeros.
#' @param norm.hp.cov  a covariance matrix whose dimension is equal to the number of regression coefficients
#'                     that gives the covariance for the normal hyperprior on the regression coefficients.
#'                     Defaults to a diagonal matrix of 1s.
#' @param iw.hp.df     degrees of freedom for the inverse-wishart prior on the covariance
#'                     matrix of the regression coefficients. Defaults to `num.predictors + 10`.
#' @param iw.hp.scale  scale matrix for inverse-Wishart prior on the covariance of the regression coefficients
#'                     Defaults to a diagonal matrix of 1s.
#' @param disp.shape   shape parameter for inverse-gamma prior on dispersion parameter for current data set
#' @param disp.scale   scale parameter for inverse-gamma prior on dispersion parameter for current data set
#' @param disp0.shape  shape parameter for inverse-gamma prior on dispersion parameter for historical data set
#' @param disp0.scale  scale parameter for inverse-gamma prior on dispersion parameter for historical data set
#' @param offset       vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0      vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...          arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return             an object of class `stanfit` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' ## take subset for speed purposes
#' actg019 = actg019[1:200, ]
#' actg036 = actg036[1:100, ]
#' glm.bhm(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(), data = actg019, histdata = actg036,
#'   chains = 1, warmup = 500, iter = 1000
#' )
glm.bhm = function(
  formula,
  family,
  data,
  histdata,
  norm.hp.mean    = NULL,
  norm.hp.cov     = NULL,
  iw.hp.df        = NULL,
  iw.hp.scale     = NULL,
  disp.shape      = 2.1,
  disp.scale      = 1.1,
  disp0.shape     = 2.1,
  disp0.scale     = 1.1,
  offset          = NULL,
  offset0         = NULL,
  ...
) {
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

  ## perform data checks
  data.checks(
    formula, family, data, histdata, offset, offset0, check.hist = TRUE
  )

  ## Default offset is vector of 0s
  if ( is.null(offset) )
    offset = rep(0, n)
  if ( is.null(offset0) )
    offset0 = rep(0, n0)

  ## Default normal hyperprior on regression coefficients is N(0, 10)
  if ( is.null(norm.hp.mean) )
    norm.hp.mean = rep(0, ncol(X))
  if ( is.null(norm.hp.cov) )
    norm.hp.cov  = diag(1, ncol(X))

  ## Default IW hyperprior on regression coefficients is IW(p + 10, 1 * I_p)
  if ( is.null(iw.hp.df) )
    iw.hp.df = p + 10
  if ( is.null(iw.hp.scale) )
    iw.hp.scale = diag(1, p)

  standat = list(
    'n'               = n,
    'n0'              = n0,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'y0'              = y0,
    'X0'              = X0,
    'hp_norm_mean'    = norm.hp.mean,
    'hp_norm_cov'     = norm.hp.cov,
    'hp_iw_df'        = iw.hp.df,
    'hp_iw_scale'     = iw.hp.scale,
    'disp_shape'      = disp.shape,
    'disp_scale'      = disp.scale,
    'disp0_shape'     = disp0.shape,
    'disp0_scale'     = disp0.scale,
    'dist'            = dist,
    'link'            = link,
    'offset'          = offset,
    'offset0'         = offset0
  )

  ## fit model in stan
  fit = rstan::sampling(stanmodels$glm_bhm, data = standat, ...)

  ## rename parameters
  oldnames = fit@sim$fnames_oi
  if ( family$family %in% c('binomial', 'poisson') ) {
    newnames = c(colnames(X), paste0( colnames(X0), '_hist' ))
  } else {
    newnames = c(colnames(X), 'dispersion')
    newnames = c( newnames, paste0(newnames, '_hist'))
  }
  fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}
