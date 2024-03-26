#' Posterior of normalized power prior (NPP) for normal linear models
#'
#' Sample from the posterior distribution of a normal linear model using the NPP by Duan et al. (2006) <doi:10.1002/env.752>.
#' The power prior parameters (\eqn{a_0}'s) are treated as random with independent beta priors. The current and historical
#' data sets are assumed to have a common dispersion parameter (\eqn{\sigma^2}) with an inverse-gamma prior. Conditional on
#' \eqn{\sigma^2}, the initial priors on the regression coefficients are independent normal distributions with variance
#' \eqn{\propto (\sigma^2)^{-1}}. In this case, the normalizing constant for the NPP has a closed form.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients.
#'                          Conditional on the variance parameter sigmasq for the outcome, beta.sd * sqrt(sigmasq) gives
#'                          the sd for the initial prior on regression coefficients. If a scalar is provided, same as for
#'                          beta.mean. Defaults to a vector of 10s.
#' @param sigmasq.shape     shape parameter for inverse-gamma prior on variance parameter. Defaults to 2.1.
#' @param sigmasq.scale     scale parameter for inverse-gamma prior on variance parameter. Defaults to 1.1.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.lower          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          lower bounds for each element of the a0 vector. If a scalar is provided, a0.lower will be a
#'                          vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param a0.upper          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          upper bounds for each element of the a0 vector. If a scalar is provided, same as for a0.lower.
#'                          Defaults to a vector of 1s.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns an object of class `draws_df` giving posterior samples.
#'
#' @references
#'  Duan, Y., Ye, K., and Smith, E. P. (2005). Evaluating water quality using power priors to incorporate historical information. Environmetrics, 17(1), 95–106.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   lm.npp(
#'     formula = cd4 ~ treatment + age + race,
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
lm.npp = function(
    formula,
    data.list,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    sigmasq.shape     = 2.1,
    sigmasq.scale     = 1.1,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks(formula, gaussian('identity'), data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  num.obs      = 1 + end.index - start.index
  p            = ncol(X)
  N            = length(y)
  K            = length(data.list) - 1 # number of historical data sets

  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  ## outcome after adjusting for offset
  y_offset = y - offset

  ## Default prior on regression coefficients given sigmasq is N(0, sigmasq * 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = to.vector(param = beta.sd, default.value = 10, len = p)
  beta.cov = diag(beta.sd^2, nrow = p, ncol = p)

  hist.mle  = matrix(NA, nrow = p, ncol = K)
  hist.prec = array(NA, c(K, p, p)) # each element is X0'X0 for a historical data set
  sumy0sq   = vector(length = K) # each element is sum(y0^2) for a historical data set (y0 has been adjusted for offsets)
  for(k in 1:K) {
    histdata       = data.list[[1+k]]
    offs0          = offset[ start.index[1+k]:end.index[1+k] ]
    histdata$offs0 = offs0
    sumy0sq[k]     = sum( y_offset[ start.index[1+k]:end.index[1+k] ]^2 )
    X0             = model.matrix(formula, histdata)
    hist.prec[k, , ] = crossprod(X0)
    fit.lm =
      suppressWarnings(
        lm(formula = formula, data = histdata, offset = offs0)
      )
    hist.mle[, k] = as.numeric(fit.lm$coefficients)
  }

  ## Default lower bound for each a0 is 0; default upper bound for each a0 is 1
  if ( !is.null(a0.lower) ){
    if ( !( is.vector(a0.lower) & (length(a0.lower) %in% c(1, K)) ) )
      stop("a0.lower must be a scalar or a vector of length ", K, " if a0.lower is not NULL")
  }
  a0.lower = to.vector(param = a0.lower, default.value = 0, len = K)
  if ( !is.null(a0.upper) ){
    if ( !( is.vector(a0.upper) & (length(a0.upper) %in% c(1, K)) ) )
      stop("a0.upper must be a scalar or a vector of length ", K, " if a0.upper is not NULL")
  }
  a0.upper = to.vector(param = a0.upper, default.value = 1, len = K)

  standat = list(
    'K'               = K, # number of historical data sets
    'n'               = num.obs[1], # current data sample size
    'n0s'             = num.obs[-1],
    'p'               = p,
    'y'               = y_offset[start.index[1]:end.index[1]], # current data outcome after adjusting for offsets
    'X'               = X[start.index[1]:end.index[1], ], # current data design matrix
    "hist_prec"       = hist.prec,
    "hist_mle"        = hist.mle,
    "sumy0sq"         = sumy0sq,
    'mean_beta'       = beta.mean,
    'cov_beta'        = beta.cov,
    'sigmasq_shape'   = sigmasq.shape,
    'sigmasq_scale'   = sigmasq.scale,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper
  )

  lm_npp     = instantiate::stan_package_model(
    name = "lm_npp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = lm_npp$sample(data = standat,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                      ...)
  d   = fit$draws(format = 'draws_df')

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  oldnames = c(oldnames, paste0('a0s[', 1:K, ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:K))
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}
