#' Log marginal likelihood of a GLM under meta-analytic predictive (MAP) prior
#'
#' @description Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under the meta-analytic predictive (MAP) prior. The MAP prior is equivalent to the prior
#' induced by the Bayesian hierarchical model (BHM).
#'
#' @description The arguments related to MCMC sampling are utilized to draw samples from the MAP prior. These
#' samples are then used to compute the logarithm of the normalizing constant of the BHM using only historical
#' data sets.
#'
#' @include data_checks.R
#' @include glm_bhm_lognc.R
#'
#' @export
#'
#' @param post.samples      output from [glm.bhm()] giving posterior samples of a GLM under the Bayesian hierarchical
#'                          model (BHM), with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup`
#'                          in `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method
#'                          in cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"glm_bhm"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood of the meta-analytic predictive (MAP) prior}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the normalizing constant of the Bayesian hierarchical model (BHM) using all data sets}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` containing the output from using
#'    [bridgesampling::bridge_sampler()] to compute the logarithm of the normalizing constant of the BHM using historical
#'    data sets}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'  }
#'
#' @references
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   formula = outcome ~ scale(age) + race + treatment + scale(cd4)
#'   family = binomial('logit')
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.bhm = glm.bhm(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.map(
#'     post.samples = d.bhm,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.map = function(
    post.samples,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')
  K         = stan.data$K
  if ( K == 1 ){
    stop("data.list should include at least one historical data set")
  }
  ## computing log normalizing constant for BHM using all data sets
  res.all = glm.bhm.lognc(
    post.samples   = post.samples,
    bridge.args    = bridge.args
  )

  ## get Stan data for BHM using historical data sets
  hist.stan.data           = stan.data
  hist.stan.data$K         = K - 1
  n                        = stan.data$end_idx[1] ## current data sample size
  hist.stan.data$N         = stan.data$N - n
  hist.stan.data$start_idx = stan.data$start_idx[-1] - n
  hist.stan.data$end_idx   = stan.data$end_idx[-1] - n
  hist.stan.data$y         = stan.data$y[-(1:n)]
  hist.stan.data$X         = stan.data$X[-(1:n), ]
  hist.stan.data$disp_mean = stan.data$disp_mean[-1]
  hist.stan.data$disp_sd   = stan.data$disp_sd[-1]
  hist.stan.data$offs      = stan.data$offs[-(1:n)]
  hist.stan.data$log_lik   = 0

  ## fit BHM using historical data sets
  glm_bhm = instantiate::stan_package_model(
    name = "glm_bhm",
    package = "hdbayes"
  )
  fit = glm_bhm$sample(data = hist.stan.data,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)
  summ = posterior::summarise_draws(fit)

  if ( hist.stan.data$dist > 2 ) {
    ## rename parameters
    K        = hist.stan.data$K
    oldnames = paste0( 'dispersion[', 1:K, ']' )
    if (K == 1) {
      newnames = 'dispersion'
    }else {
      newnames = c('dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
    }
    hist.post.samples = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  }else {
    hist.post.samples = fit$draws(format = 'draws_df')
  }
  attr(x = hist.post.samples, which = 'data') = hist.stan.data

  ## compute log normalizing constant for BHM using historical data sets
  res.hist = glm.bhm.lognc(
    post.samples   = hist.post.samples,
    bridge.args    = bridge.args
  )

  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## the minimum estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat
  res = list(
    'model'        = "glm_bhm",
    'logml'        = res.all$lognc - res.hist$lognc,
    'bs'           = res.all$bs,
    'bs.hist'      = res.hist$bs,
    'min_ess_bulk' = min(summ[, 'ess_bulk']),
    'max_Rhat'     = max(summ[, 'rhat'])
  )

  if ( res[['min_ess_bulk']] < 1000 )
    warning(
      paste0(
        'The minimum bulk effective sample size of the MCMC sampling is ',
        round(res[['min_ess_bulk']], 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res[['max_Rhat']] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        round(res[['max_Rhat']], 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}
