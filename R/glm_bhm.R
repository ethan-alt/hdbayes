#' Posterior of Bayesian hierarchical model (BHM)
#'
#' Sample from the posterior distribution of a GLM using the BHM.
#'
#' The BHM assumes that the regression coefficients for the historical and current data are different, but
#' are correlated through a common distribution, whose hyperparameters (i.e., mean and standard deviation (sd)
#' (the covariance matrix is assumed to have a diagonal structure)) are treated as random. The number of
#' of regression coefficients for the current data is assumed to be the same as that for the historical data.
#'
#' The hyperpriors on the mean and the sd hyperparameters are independent normal and independent half-normal
#' distributions, respectively. The priors on the dispersion parameters (if applicable) for the current and
#' historical data sets are independent half-normal distributions.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param meta.mean.mean    a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the normal hyperpriors on the mean hyperparameters of regression coefficients. If
#'                          a scalar is provided, meta.mean.mean will be a vector of repeated elements of the given scalar.
#'                          Defaults to a vector of 0s.
#' @param meta.mean.sd      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the normal hyperpriors on the mean hyperparameters of regression coefficients. If a
#'                          scalar is provided, same as for meta.mean.mean. Defaults to a vector of 10s.
#' @param meta.sd.mean      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for meta.mean.mean. Defaults to a vector of 0s.
#' @param meta.sd.sd        a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for meta.mean.mean. Defaults to a vector of 1s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the means for the half-normal priors on the dispersion parameters. If a scalar is
#'                          provided, same as for meta.mean.mean. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the sds for the half-normal priors on the dispersion parameters. If a scalar is
#'                          provided, same as for meta.mean.mean. Defaults to a vector of 10s.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          [cmdstanr::sample()].
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in [cmdstanr::sample()].
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in [cmdstanr::sample()].
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns an object of class `draws_df` giving posterior samples.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:200, ]
#'   actg036 = actg036[1:100, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.bhm(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.bhm = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks(formula, family, data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  ## Default normal hyperprior on mean of regression coefficients is N(0, 10^2)
  if ( !is.null(meta.mean.mean) ){
    if ( !( is.vector(meta.mean.mean) & (length(meta.mean.mean) %in% c(1, p)) ) )
      stop("meta.mean.mean must be a scalar or a vector of length ", p, " if meta.mean.mean is not NULL")
  }
  meta.mean.mean = to.vector(param = meta.mean.mean, default.value = 0, len = p)
  if ( !is.null(meta.mean.sd) ){
    if ( !( is.vector(meta.mean.sd) & (length(meta.mean.sd) %in% c(1, p)) ) )
      stop("meta.mean.sd must be a scalar or a vector of length ", p, " if meta.mean.sd is not NULL")
  }
  meta.mean.sd = to.vector(param = meta.mean.sd, default.value = 10, len = p)

  ## Default half-normal hyperprior on sd of regression coefficients is N^{+}(0, 1)
  if ( !is.null(meta.sd.mean) ){
    if ( !( is.vector(meta.sd.mean) & (length(meta.sd.mean) %in% c(1, p)) ) )
      stop("meta.sd.mean must be a scalar or a vector of length ", p, " if meta.sd.mean is not NULL")
  }
  meta.sd.mean = to.vector(param = meta.sd.mean, default.value = 0, len = p)
  if ( !is.null(meta.sd.sd) ){
    if ( !( is.vector(meta.sd.sd) & (length(meta.sd.sd) %in% c(1, p)) ) )
      stop("meta.sd.sd must be a scalar or a vector of length ", p, " if meta.sd.sd is not NULL")
  }
  meta.sd.sd = to.vector(param = meta.sd.sd, default.value = 1, len = p)

  ## Default half-normal hyperprior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) %in% c(1, K)) ) )
      stop("disp.mean must be a scalar or a vector of length ", K, " if disp.mean is not NULL")
  }
  disp.mean = to.vector(param = disp.mean, default.value = 0, len = K)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) %in% c(1, K)) ) )
      stop("disp.sd must be a scalar or a vector of length ", K, " if disp.sd is not NULL")
  }
  disp.sd = to.vector(param = disp.sd, default.value = 10, len = K)


  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'meta_mean_mean'  = meta.mean.mean,
    'meta_mean_sd'    = meta.mean.sd,
    'meta_sd_mean'    = meta.sd.mean,
    'meta_sd_sd'      = meta.sd.sd,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )

  glm_bhm    = instantiate::stan_package_model(
    name = "glm_bhm",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_bhm$sample(data = standat,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)
  pars = fit$metadata()$model_params

  ## rename parameters
  oldnames = paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( K == 1 ) {
    newnames = colnames(X)
  }else {
    newnames = c(colnames(X), paste0( colnames(X), '_hist_', rep(1:(K-1), each = p) ))
  }
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    if (K == 1) {
      newnames = c(newnames, 'dispersion')
    }else {
      newnames = c(newnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
    }
  }
  ## reorder parameters so that regression coefficients appear at the top
  pars = c(pars[1], pars[pars %in% oldnames], (pars[!pars %in% oldnames])[-1])
  d   = fit$draws(format = 'draws_df', variables = pars)
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}
