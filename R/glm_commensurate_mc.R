#' Posterior of Commensurate Prior with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the commensurate prior
#' with Gaussian copula imputation for missing covariate values.
#'
#' This function extends \code{\link{glm.commensurate}} to handle missing values in covariates
#' using a Gaussian copula model.
#'
#' @include data_checks_mc.R
#' @include get_stan_data_mc.R
#' @include missing_data_utils.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element is the current data, the rest are historical.
#'                          Missing values (NA) are allowed in covariates.
#' @param copula.corr       an optional correlation matrix for the Gaussian copula. If NULL (default), estimated from data.
#' @param offset.list       a list of vectors giving the offsets. Defaults to NULL.
#' @param beta0.mean        mean for the normal initial prior on historical data coefficients. Defaults to 0.
#' @param beta0.sd          sd for the normal initial prior on historical data coefficients. Defaults to 10.
#' @param disp.mean         location parameter for the half-normal prior on dispersion. Defaults to 0.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion. Defaults to 10.
#' @param p.spike           prior probability for the spike component. Defaults to 0.1.
#' @param spike.mean        mean for the half-normal spike prior on commensurability. Defaults to 200.
#' @param spike.sd          sd for the half-normal spike prior on commensurability. Defaults to 0.1.
#' @param slab.mean         mean for the half-normal slab prior on commensurability. Defaults to 0.
#' @param slab.sd           sd for the half-normal slab prior on commensurability. Defaults to 5.
#' @param get.loglik        whether to generate log-likelihood matrix. Defaults to FALSE.
#' @param iter_warmup       number of warmup iterations. Defaults to 1000.
#' @param iter_sampling     number of post-warmup iterations. Defaults to 1000.
#' @param chains            number of Markov chains. Defaults to 4.
#' @param ...               arguments passed to `sample()` method in cmdstanr package.
#'
#' @return
#'  The function returns an object of class `draws_df` containing posterior samples.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   actg019$cd4[sample(1:100, 10)] = NA
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.commensurate.mc(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.commensurate.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  standat = get.stan.data.cp.mc(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    copula.corr = copula.corr,
    offset.list = offset.list,
    beta0.mean  = beta0.mean,
    beta0.sd    = beta0.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd,
    p.spike     = p.spike,
    spike.mean  = spike.mean,
    spike.sd    = spike.sd,
    slab.mean   = slab.mean,
    slab.sd     = slab.sd,
    get.loglik  = get.loglik
  )

  glm_commensurate_mc = instantiate::stan_package_model(
    name = "glm_commensurate_mc",
    package = "hdbayes"
  )

  fit = glm_commensurate_mc$sample(data = standat,
                                   iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                   ...)

  p = standat$p
  X = standat$X
  K = standat$K
  oldnames = c(
    paste0("beta[", 1:p, "]"),
    paste0("beta0[", 1:p, "]"),
    paste0("tau[", 1:p, "]")
  )
  newnames = c(
    colnames(X),
    paste0(colnames(X), "_hist"),
    paste0("tau_", colnames(X))
  )
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0('dispersion[', 1:K, ']'))
    if (K == 1) {
      newnames = c(newnames, 'dispersion')
    } else {
      newnames = c(newnames, 'dispersion', paste0('dispersion_hist_', 1:(K-1)))
    }
  }
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  attr(x = d, which = 'data') = standat
  attr(x = d, which = 'model') = "glm_commensurate_mc"
  return(d)
}
