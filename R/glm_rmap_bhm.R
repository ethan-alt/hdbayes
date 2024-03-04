#' Posterior of robust meta-analytic predictive prior (RMAP)
#'
#' First step for sampling from the posterior distribution of a GLM using the RMAP by Schmidli et al. (2014) <doi:10.1111/biom.12242>.
#'
#' The RMAP is a mixture prior of two components where one component is a prior induced by the Bayesian hierarchical model (BHM),
#' and the other is a vague (noninformative) prior. This function samples from the prior induced by the BHM.
#'
#' @include data_checks.R
#' @include glm_bhm.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param hist.data.list    a list of `data.frame`s. Each element in the list is a historical data set.
#' @param hist.offset.list  a list of vectors giving the offsets for each historical data. The length of hist.offset.list is
#'                          equal to the length of hist.data.list. The length of each element of hist.offset.list is equal
#'                          to the number of rows in the corresponding element of hist.data.list. Defaults to a list of
#'                          vectors of 0s.
#' @param meta.mean.mean    a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the normal hyperpriors on the mean hyperparameters of regression coefficients. If
#'                          a scalar is provided, meta.mean.mean will be a vector of repeated elements of the given scalar.
#'                          Defaults to a vector of 0s.
#' @param meta.mean.sd      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the normal hyperpriors on the mean hyperparameters of regression coefficients. If a
#'                          scalar is provided, same as for meta.mean.mean. Defaults to a vector of 1s.
#' @param meta.sd.mean      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for meta.mean.mean. Defaults to a vector of 0s.
#' @param meta.sd.sd        a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for meta.mean.mean. Defaults to a vector of 10s.
#' @param hist.disp.mean    a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          means for the half-normal hyperpriors on the dispersion parameters. If a scalar is provided,
#'                          same as for meta.mean.mean. Defaults to a vector of 0s.
#' @param hist.disp.sd      a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          sds for the half-normal hyperpriors on the dispersion parameters. If a scalar is provided, same
#'                          as for meta.mean.mean. Defaults to a vector of 10s.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          [cmdstanr::sample()].
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in [cmdstanr::sample()].
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in [cmdstanr::sample()].
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns a matrix of the samples of regression coefficients from the prior induced by the BHM.
#'  The number of columns is equal to the number of regression coefficients, and the number of rows is equal to
#'  the number of MCMC samples.
#'
#' @seealso [glm.rmap.bhm.approx()] for the second step and [glm.rmap()] for the final step of implementing RMAP.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036) ## historical data
#'   ## take subset for speed purposes
#'   actg036 = actg036[1:50, ]
#'   hist_data_list = list(actg036)
#'   glm.rmap.bhm(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     hist.data.list = hist_data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.rmap.bhm = function(
    formula,
    family,
    hist.data.list,
    hist.offset.list  = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    hist.disp.mean    = NULL,
    hist.disp.sd      = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## perform data checks
  suppressMessages(data.checks(formula, family, hist.data.list, hist.offset.list))

  ## fit BHM using all historical data
  hist.bhm = glm.bhm(
    formula           = formula,
    family            = family,
    data.list         = hist.data.list,
    offset.list       = hist.offset.list,
    meta.mean.mean    = meta.mean.mean,
    meta.mean.sd      = meta.mean.sd,
    meta.sd.mean      = meta.sd.mean,
    meta.sd.sd        = meta.sd.sd,
    disp.mean         = hist.disp.mean,
    disp.sd           = hist.disp.sd,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    chains            = chains,
    ...
  )

  varnames = posterior::variables(hist.bhm)
  suppressWarnings({
    beta_mean_samples = as.matrix( hist.bhm[, varnames[grepl("beta_mean", varnames)]] )
    beta_sd_samples   = as.matrix( hist.bhm[, varnames[grepl("beta_sd", varnames)]] )
  })
  p = ncol(beta_mean_samples)

  ## sampling from the prior induced by the BHM
  ## equivalent to sampling from the meta-analytic predictive (MAP) prior
  beta_pred = sapply(1:p, function(j){
    stats::rnorm(
      n = nrow(beta_mean_samples),
      mean = beta_mean_samples[, j],
      sd = beta_sd_samples[, j]
    )
  })
  beta_pred           = as.matrix(beta_pred)
  colnames(beta_pred) = varnames[2:(1+p)]
  return(
    list(
      beta_pred = beta_pred,
      hist_bhm  = hist.bhm ## can be used for assessing MCMC convergence
      )
  )
}
