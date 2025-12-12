#' Posterior of Bayesian hierarchical model (BHM) with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the Bayesian hierarchical model (BHM)
#' with Gaussian copula imputation for missing covariate values.
#'
#' This function extends \code{\link{glm.bhm}} to handle missing values in covariates using a
#' Gaussian copula model. The copula approach:
#' \itemize{
#'   \item Models dependencies between covariates via a correlation matrix
#'   \item Handles mixed variable types (continuous, binary, ordinal, count)
#'   \item Propagates imputation uncertainty through the posterior
#' }
#'
#' The Bayesian hierarchical model (BHM) assumes that the regression coefficients for the historical and
#' current data are different, but are correlated through a common distribution, whose hyperparameters
#' (i.e., mean and standard deviation (sd)) are treated as random.
#'
#' @include data_checks_mc.R
#' @include get_stan_data_mc.R
#' @include missing_data_utils.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. Missing values (NA) are allowed in covariates but not in
#'                          the response variable.
#' @param copula.corr       an optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases in the data.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of `offset.list` is equal to
#'                          the length of `data.list`. The length of each element of `offset.list` is equal to the number
#'                          of rows in the corresponding element of `data.list`. Defaults to a list of vectors of 0s.
#' @param meta.mean.mean    a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the normal hyperpriors on the mean hyperparameters of regression coefficients.
#'                          If a scalar is provided, `meta.mean.mean` will be a vector of repeated elements of the given
#'                          scalar. Defaults to a vector of 0s.
#' @param meta.mean.sd      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the normal hyperpriors on the mean hyperparameters of regression coefficients. If
#'                          a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 10s.
#' @param meta.sd.mean      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 0s.
#' @param meta.sd.sd        a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 1s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the location parameters for the half-normal priors on the dispersion parameters.
#'                          If a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the scale parameters for the half-normal priors on the dispersion parameters. If a
#'                          scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 10s.
#' @param get.loglik        whether to generate log-likelihood matrix. Defaults to FALSE.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
#' @return
#'  The function returns an object of class `draws_df` containing posterior samples. The object has two attributes:
#'
#'  \describe{
#'    \item{data}{a list of variables specified in the data block of the Stan program}
#'
#'    \item{model}{a character string indicating the model name ("glm_bhm_mc")}
#'  }
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   ## introduce some missing values
#'   actg019$cd4[sample(1:100, 10)] = NA
#'   actg036$age[sample(1:50, 5)] = NA
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.bhm.mc(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.bhm.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for BHM with missing covariates
  standat = get.stan.data.bhm.mc(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    copula.corr    = copula.corr,
    offset.list    = offset.list,
    meta.mean.mean = meta.mean.mean,
    meta.mean.sd   = meta.mean.sd,
    meta.sd.mean   = meta.sd.mean,
    meta.sd.sd     = meta.sd.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    get.loglik     = get.loglik
  )

  glm_bhm_mc = instantiate::stan_package_model(
    name = "glm_bhm_mc",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_bhm_mc$sample(data = standat,
                          iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                          ...)

  ## rename parameters
  p        = standat$p
  K        = standat$K
  X        = standat$X
  oldnames = paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( K == 1 ) {
    newnames = colnames(X)
  } else {
    newnames = c(colnames(X), paste0( colnames(X), '_hist_', rep(1:(K-1), each = p) ))
  }
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    if (K == 1) {
      newnames = c(newnames, 'dispersion')
    } else {
      newnames = c(newnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
    }
  }
  ## reorder parameters so that regression coefficients appear at the top
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  ## add model name as an attribute
  attr(x = d, which = 'model') = "glm_bhm_mc"
  return(d)
}
