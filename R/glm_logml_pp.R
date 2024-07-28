#' Log marginal likelihood of a GLM under power prior (PP)
#'
#' @description Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under the power prior (PP).
#'
#' @description The arguments related to MCMC sampling are utilized to draw samples from the power prior (PP).
#' These samples are then used to compute the logarithm of the normalizing constant of the PP using only historical
#' data sets.
#'
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include glm_logml_post.R
#' @include glm_pp_lognc.R
#'
#' @export
#'
#' @param post.samples      output from [glm.pp()] giving posterior samples of a GLM under the power prior (PP), with
#'                          an attribute called 'data' which includes the list of variables specified in the data block
#'                          of the Stan program.
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
#'  If all of the power prior parameters (\eqn{a_0}'s) are equal to zero, or if the posterior samples are obtained from using only
#'  one data set (the current data), then the function will return the same result as the output from [glm.logml.post()].
#'
#'  If at least one of the power prior parameters (\eqn{a_0}'s) is non-zero, the function will return a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"PP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the normalizing constant of the power prior (PP) using all data sets}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` containing the output from using
#'    [bridgesampling::bridge_sampler()] to compute the logarithm of the normalizing constant of the PP using historical
#'    data sets}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'  }
#'
#' @references
#'  Chen, M.-H. and Ibrahim, J. G. (2000). Power prior distributions for Regression Models. Statistical Science, 15(1).
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   formula = cd4 ~ treatment + age + race
#'   family = poisson('log')
#'   a0 = 0.5
#'   d.pp = glm.pp(
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     a0.vals = a0,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.pp(
#'     post.samples = d.pp,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.pp = function(
    post.samples,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')
  K         = stan.data$K
  a0_vals   = stan.data$a0_vals

  if ( K == 1 || sum(a0_vals) == 1 ){
    ## PP is equivalent to normal/half-normal prior in this case
    stan.data.post   = stan.data[c('p', 'mean_beta', 'sd_beta', 'disp_mean', 'disp_sd',
                                'dist', 'link')]
    n                = stan.data$end_idx[1] ## current data sample size
    stan.data.post$n = n
    stan.data.post$y = stan.data$y[1:n]
    stan.data.post$X = stan.data$X[1:n, ]
    stan.data.post$offs = stan.data$offs[1:n]

    attr(post.samples, 'data') = stan.data.post
    res = glm.logml.post(
      post.samples = post.samples,
      bridge.args = bridge.args
    )
    return(res)
  }

  ## computing log normalizing constant for power prior using all data sets
  res.all = glm.pp.lognc(
    post.samples   = post.samples,
    bridge.args    = bridge.args
  )

  ## get Stan data for PP using historical data sets
  hist.stan.data           = stan.data
  hist.stan.data$K         = K - 1
  n                        = stan.data$end_idx[1] ## current data sample size
  hist.stan.data$N         = stan.data$N - n
  hist.stan.data$start_idx = stan.data$start_idx[-1] - n
  hist.stan.data$end_idx   = stan.data$end_idx[-1] - n
  hist.stan.data$y         = stan.data$y[-(1:n)]
  hist.stan.data$X         = stan.data$X[-(1:n), ]
  hist.stan.data$offs      = stan.data$offs[-(1:n)]
  hist.stan.data$a0_vals   = a0_vals[-1]

  ## fit PP using historical data sets
  glm_pp     = instantiate::stan_package_model(
    name = "glm_pp",
    package = "hdbayes"
  )
  fit = glm_pp$sample(data = hist.stan.data,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                      ...)
  summ = posterior::summarise_draws(fit)

  if ( hist.stan.data$dist > 2 ) {
    ## rename parameters
    oldnames = 'dispersion[1]'
    newnames = 'dispersion'
    hist.post.samples = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  }else {
    hist.post.samples = fit$draws(format = 'draws_df')
  }
  attr(x = hist.post.samples, which = 'data') = hist.stan.data

  ## compute log normalizing constant for PP using historical data sets
  res.hist = glm.pp.lognc(
    post.samples   = hist.post.samples,
    bridge.args    = bridge.args
  )

  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## the minimum estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat
  res = list(
    'model'        = "PP",
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
