#' Posterior of power prior (PP) with fixed \eqn{a_0}
#'
#' Sample from the posterior distribution of a GLM using the PP by Ibrahim and Chen (2000) <doi:10.1214/ss/1009212673>.
#'
#' The power prior parameters (\eqn{a_0}'s) are treated as fixed. The initial priors on the regression coefficients
#' are independent normal priors. The current and historical data sets are assumed to have a common dispersion parameter
#' with a half-normal prior (if applicable).
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
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for beta.mean. Defaults to a vector of 10s.
#' @param a0.vals           a scalar between 0 and 1 or a vector whose dimension is equal to the number of historical
#'                          data sets giving the (fixed) power prior parameter for each historical data set. Each element of
#'                          vector should be between 0 and 1. If a scalar is provided, same as for beta.mean.
#' @param disp.mean         mean parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           sd parameter for the half-normal prior on dispersion parameter. Defaults to 10.
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
#'  Chen, M.-H. and Ibrahim, J. G. (2000). Power prior distributions for Regression Models. Statistical Science, 15(1).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.pp(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data_list,
#'     a0.vals = 0.5,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.pp = function(
    formula,
    family,
    data.list,
    a0.vals,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
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

  ## Default prior on regression coefficients is N(0, 10^2)
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

  ## check a0 values
  if ( !( is.vector(a0.vals) & (length(a0.vals) %in% c(1, K-1)) ) )
    stop("a0.vals must be a scalar or a vector of length ", K-1)
  a0.vals = to.vector(param = a0.vals, len = K-1)
  if ( any(a0.vals < 0 | a0.vals > 1 ) )
    stop("Each element of a0.vals must be a scalar between 0 and 1")
  a0.vals = c(1, a0.vals) # first element = 1 for current data

  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) == 1) ) )
      stop("disp.mean must be a scalar if disp.mean is not NULL")
  }
  disp.mean = to.vector(param = disp.mean, default.value = 0, len = 1)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) == 1) ) )
      stop("disp.sd must be a scalar if disp.sd is not NULL")
  }
  disp.sd = to.vector(param = disp.sd, default.value = 10, len = 1)


  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta.mean,
    'sd_beta'         = beta.sd,
    'a0_vals'         = a0.vals,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )

  glm_pp     = instantiate::stan_package_model(
    name = "glm_pp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_pp$sample(data = standat,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                      ...)
  d   = fit$draws(format = 'draws_df')

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}
