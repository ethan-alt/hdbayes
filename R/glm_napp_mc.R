#' Posterior of normalized asymptotic power prior (NAPP) with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the normalized asymptotic power prior (NAPP) by
#' Ibrahim et al. (2015) <doi:10.1002/sim.6728>, extended to handle missing covariates via a Gaussian copula model.
#'
#' The normalized asymptotic power prior (NAPP) assumes that the regression coefficients and logarithm of the
#' dispersion parameter are a multivariate normal distribution with mean equal to the maximum likelihood
#' estimate of the historical data and covariance matrix equal to \eqn{a_0^{-1}} multiplied by the inverse Fisher
#' information matrix of the historical data.
#'
#' For data with missing covariates, the MLE and Fisher information from historical data are computed using
#' complete cases. The missing covariate imputation via Gaussian copula is applied to the current data during
#' posterior sampling.
#'
#' @include data_checks_mc.R
#' @include get_stan_data_mc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical datasets. The current data may contain missing covariate values.
#' @param copula.corr       optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases across all datasets.
#' @param offset.list       a list of vectors giving the offsets for each data.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector.
#' @param get.loglik        whether to generate log-likelihood matrix. Defaults to FALSE.
#' @param iter_warmup       number of warmup iterations to run per chain.
#' @param iter_sampling     number of post-warmup iterations to run per chain.
#' @param chains            number of Markov chains to run.
#' @param ...               arguments passed to `sample()` method in cmdstanr package.
#'
#' @return
#'  The function returns an object of class `draws_df` containing posterior samples.
#'
#' @references
#'  Ibrahim, J. G., Chen, M., Gwon, Y., and Chen, F. (2015). The power prior: Theory and applications. Statistics in Medicine, 34(28), 3724-3749.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   ## introduce missing values in current data
#'   actg019$age[sample(1:100, 10)] = NA
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.napp.mc(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.napp.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  # Check if there are missing covariates
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    # No missing data - use the standard NAPP function
    message("No missing covariates detected. Consider using glm.napp() for better efficiency.")
    return(glm.napp(
      formula = formula,
      family = family,
      data.list = data.list,
      offset.list = offset.list,
      a0.shape1 = a0.shape1,
      a0.shape2 = a0.shape2,
      get.loglik = get.loglik,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      ...
    ))
  }

  # For NAPP with missing covariates, we use complete cases for the historical
  # data MLE computation and impute missing values in the current data
  message("Note: NAPP with missing covariates computes MLE from complete historical cases and ",
          "imputes missing values in the current data.")

  ## get Stan data for NAPP (uses complete cases for historical data)
  standat = get.stan.data.napp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    offset.list = offset.list,
    a0.shape1   = a0.shape1,
    a0.shape2   = a0.shape2,
    get.loglik  = get.loglik
  )

  glm_napp = instantiate::stan_package_model(
    name = "glm_napp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_napp$sample(data = standat,
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
  oldnames = c(oldnames, paste0('a0s[', 1:K, ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:K))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  attr(x = d, which = 'data') = standat
  attr(x = d, which = 'model') = "glm_napp_mc"
  return(d)
}
