#' Posterior with normal/half-normal prior and missing covariates
#'
#' Sample from the posterior distribution of a GLM using the normal/half-normal prior
#' (reference prior) with Gaussian copula imputation for missing covariate values.
#'
#' This function extends \code{\link{glm.post}} to handle missing values in covariates using a
#' Gaussian copula model.
#'
#' @include data_checks_mc.R
#' @include get_stan_data_mc.R
#' @include missing_data_utils.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. Only the first element is used (current data).
#'                          Missing values (NA) are allowed in covariates.
#' @param copula.corr       an optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases.
#' @param offset.list       a list of vectors giving the offsets. Defaults to NULL.
#' @param beta.mean         mean for the normal prior on regression coefficients. Defaults to 0.
#' @param beta.sd           sd for the normal prior on regression coefficients. Defaults to 10.
#' @param disp.mean         location parameter for the half-normal prior on dispersion. Defaults to 0.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion. Defaults to 10.
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
#'   actg019 = actg019[1:100, ]
#'   actg019$cd4[sample(1:100, 10)] = NA
#'   data_list = list(currdata = actg019)
#'   glm.post.mc(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.post.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  standat = get.stan.data.post.mc(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    copula.corr = copula.corr,
    offset.list = offset.list,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd,
    get.loglik  = get.loglik
  )

  glm_post_mc = instantiate::stan_package_model(
    name = "glm_post_mc",
    package = "hdbayes"
  )

  fit = glm_post_mc$sample(data = standat,
                           iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                           ...)

  p = standat$p
  X = standat$X
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  attr(x = d, which = 'data') = standat
  attr(x = d, which = 'model') = "glm_post_mc"
  return(d)
}
