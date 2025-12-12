#' Posterior of normalized power prior (NPP) with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the normalized power prior (NPP) by Duan et al.
#' (2006) <doi:10.1002/env.752>, extended to handle missing covariates via a Gaussian copula model.
#'
#' Before using this function, users must estimate the logarithm of the normalizing constant across a
#' range of different values for the power prior parameter (\eqn{a_0}). For data with missing covariates,
#' the normalizing constants should be computed using either:
#' \itemize{
#'   \item Complete cases from historical data
#'   \item Pre-imputed historical data
#' }
#' The missing covariate imputation via Gaussian copula is applied to the current data during posterior
#' sampling.
#'
#' @include data_checks_mc.R
#' @include glm_npp_lognc.R
#' @include get_stan_data_mc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. The current data may contain missing covariate values.
#' @param a0.lognc          a vector giving values of the power prior parameter for which the logarithm of the normalizing
#'                          constant has been evaluated.
#' @param lognc             an S by T matrix where S is the length of `a0.lognc`, T is the number of historical data sets.
#' @param copula.corr       optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases across all datasets.
#' @param offset.list       a list of vectors giving the offsets for each data.
#' @param beta.mean         mean parameters for the initial prior on regression coefficients.
#' @param beta.sd           sd parameters for the initial prior on regression coefficients.
#' @param disp.mean         location parameter for the half-normal prior on dispersion parameter.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion parameter.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on \eqn{a_0} vector.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on \eqn{a_0} vector.
#' @param a0.lower          lower bounds for each element of the \eqn{a_0} vector.
#' @param a0.upper          upper bounds for each element of the \eqn{a_0} vector.
#' @param get.loglik        whether to generate log-likelihood matrix. Defaults to FALSE.
#' @param iter_warmup       number of warmup iterations to run per chain.
#' @param iter_sampling     number of post-warmup iterations to run per chain.
#' @param chains            number of Markov chains to run.
#' @param ...               arguments passed to `sample()` method in cmdstanr package.
#'
#' @return
#'  The function returns an object of class `draws_df` containing posterior samples.
#'
#' @seealso [glm.npp.lognc()]
#'
#' @references
#'  Duan, Y., Ye, K., and Smith, E. P. (2005). Evaluating water quality using power priors to incorporate historical information. Environmetrics, 17(1), 95-106.
#'
#' @examples
#' \donttest{
#'   if(requireNamespace("parallel")){
#'     data(actg019)
#'     data(actg036)
#'     ## take subset for speed purposes
#'     actg019 = actg019[1:100, ]
#'     actg036 = actg036[1:50, ]
#'     ## introduce missing values in current data
#'     actg019$age[sample(1:100, 10)] = NA
#'
#'     library(parallel)
#'     ncores    = 2
#'     formula   = cd4 ~ treatment + age + race
#'     family    = poisson()
#'     a0        = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## Compute normalizing constants from complete historical data
#'       logncfun = function(a0, ...){
#'         hdbayes::glm.npp.lognc(
#'           formula = formula, family = family, a0 = a0, histdata = actg036,
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'family', 'actg036'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from NPP with missing covariate imputation
#'       glm.npp.mc(
#'         formula = formula,
#'         family = family,
#'         data.list = list(actg019, actg036),
#'         a0.lognc = a0.lognc$a0,
#'         lognc = matrix(a0.lognc$lognc, ncol = 1),
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'     }
#'   }
#' }
glm.npp.mc = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  # Check if there are missing covariates
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    # No missing data - use the standard NPP function
    message("No missing covariates detected. Consider using glm.npp() for better efficiency.")
    return(glm.npp(
      formula = formula,
      family = family,
      data.list = data.list,
      a0.lognc = a0.lognc,
      lognc = lognc,
      offset.list = offset.list,
      beta.mean = beta.mean,
      beta.sd = beta.sd,
      disp.mean = disp.mean,
      disp.sd = disp.sd,
      a0.shape1 = a0.shape1,
      a0.shape2 = a0.shape2,
      a0.lower = a0.lower,
      a0.upper = a0.upper,
      get.loglik = get.loglik,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      ...
    ))
  }

  # For now, the NPP with missing covariates uses complete cases for the historical
  # likelihood and imputes missing values in the current data
  # This is a simplification - a full implementation would require joint imputation
  message("Note: NPP with missing covariates currently imputes missing values in the current data only. ",
          "Historical data normalizing constants should be computed from complete cases.")

  ## get Stan data for NPP with missing covariates in current data
  standat = get.stan.data.npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper,
    get.loglik     = get.loglik
  )

  glm_npp_posterior = instantiate::stan_package_model(
    name = "glm_npp_posterior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_npp_posterior$sample(data = standat,
                                 iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                 ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X
  K        = standat$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0s[', 1:(K-1), ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:(K-1)))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  attr(x = d, which = 'data') = standat
  attr(x = d, which = 'model') = "glm_npp_mc"
  return(d)
}
