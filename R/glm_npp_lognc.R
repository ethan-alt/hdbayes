
#'
#' Estimate the logarithm of the normalizing constant for normalized power prior
#'
#' Uses Markov chain Monte Carlo and bridge sampling
#' to estimate the logarithm of the normalizing
#' constant for the normalized power prior for a fixed value of a0. The discounts
#' the historical data likelihood by a value a0 between 0 and 1. The initial prior
#' is a multivariate-normal inverse-gamma prior on the regression coefficients
#' (conditional on the dispersion) and the dispersion parameter (if applicable).
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula      a two-sided formula giving the relationship between the response variable and covariates
#' @param family       an object of class `family`. See [stats::family()]
#' @param histdata     a `data.frame` giving the historical data
#' @param a0           the power prior parameter (a scalar between 0 and 1)
#' @param beta.mean    mean parameter for the normal initial prior on the regression coefficients given the dispersion parameter. Defaults to a vector of 0s
#' @param beta.cov     covariance parameter for the conditional multivariate normal initial prior on the regression coefficients. The covariance used is \code{dispersion * beta.cov}. Defaults to a diagonal matrix of 100s
#' @param disp.shape   shape parameter for inverse-gamma initial prior on dispersion parameter
#' @param disp.scale   scale parameter for inverse-gamma initial prior on dispersion parameter
#' @param offset0      vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param bridge.args  <optional> a `list` giving arguments to pass onto [bridgesampling::bridge_sampler()]
#' @param ...          arguments passed to [rstan::sampling()] (e.g. iter, chains)
#'
#' @return             a vector giving the value of a0, the estimated logarithm of the normalizing constant, the minimum effective sample size of
#'                     the MCMC sampling, and the maximum Rhat.
#'
#' @examples
#' data(actg036)
#' ## take subset for speed purposes
#' actg036 = actg036[1:50, ]
#' glm.npp.lognc(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(), histdata = actg036, a0 = 0.5,
#'   chains = 1, warmup = 500, iter = 5000
#' )
#'
#'

glm.npp.lognc = function(
  formula,
  family,
  histdata,
  a0,
  beta.mean   = NULL,
  beta.cov    = NULL,
  disp.shape  = 2.1,
  disp.scale  = 1.1,
  offset0     = NULL,
  bridge.args = NULL,
  ...
) {
  y0 = histdata[, all.vars(formula)[1]]
  n0 = length(y0)
  X0 = model.matrix(formula, histdata)
  p  = ncol(X0)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  if (length(a0) != 1)
    stop('a0 must be a scalar')
  if ( a0 < 0 | a0 > 1 )
    stop("a0 must be between 0 and 1")

  ## Default offset is vector of 0s
  if ( is.null(offset0) )
    offset0 = rep(0, n0)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) )
    beta.mean = rep(0, ncol(X0))
  if ( is.null(beta.cov) )
    beta.cov  = diag(100, ncol(X0))

  standat = list(
    'n0'          = n0,
    'p'           = p,
    'y0'          = y0,
    'X0'          = X0,
    'beta_mean'   = beta.mean,
    'beta_cov'    = beta.cov,
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'a0'          = a0,
    'dist'        = dist,
    'link'        = link,
    'offset0'     = offset0
  )

  ## perform data checks
  data.checks(
    formula, family, data = histdata, histdata = NULL, offset = offset0, offset0 = NULL, check.hist = FALSE
  )

  ## fit in rstan; obtain max Rhat and min n_eff
  fit  = rstan::sampling(stanmodels$glm_npp_prior, data = standat, ...)
  summ = rstan::summary(fit)$summary

  ## estimate log normalizing constant
  bs = do.call(what = bridgesampling::bridge_sampler, args = c('samples' = fit, as.list(bridge.args)))

  ## Return vector of a0, lognc, min_n_eff, max_Rhat
  res = c(
    'a0'        = a0,
    'lognc'     = bs$logml,
    'min_n_eff' = min(summ[, 'n_eff']),
    'max_Rhat'  = max(summ[, 'Rhat'])
  )
  if ( res['min_n_eff'] < 1000 )
    warning(
      paste0(
        'The minimum effective sample size of the MCMC sampling is ',
        res['min_n_eff'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res['max_Rhat'] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        res['max_Rhat'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}


