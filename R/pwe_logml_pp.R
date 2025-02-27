#' Log marginal likelihood of a piecewise exponential (PWE) model under the power prior (PP)
#'
#' @description Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a PWE model under the power prior (PP).
#'
#' @description The arguments related to MCMC sampling are utilized to draw samples from the power prior (PP).
#' These samples are then used to compute the logarithm of the normalizing constant of the PP using only historical
#' data set.
#'
#' @include pwe_logml_post.R
#' @include pwe_pp_lognc.R
#'
#' @export
#'
#' @param post.samples      output from [pwe.pp()] giving posterior samples of a PWE model under the power prior (PP), with
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
#'  If the power prior parameter (\eqn{a_0}) is equal to zero, then the function will return the same result as the output
#'  from [pwe.logml.post()].
#'
#'  If the power prior parameters (\eqn{a_0}) is non-zero, the function will return a `list` with the following objects
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
#'    data set}
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
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1684)
#'     data(E1690)
#'     ## take subset for speed purposes
#'     E1684 = E1684[1:100, ]
#'     E1690 = E1690[1:50, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'     E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'     E1684$cage = as.numeric(scale(E1684$age))
#'     E1690$cage = as.numeric(scale(E1690$age))
#'     data_list = list(currdata = E1690, histdata = E1684)
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     d.pp = pwe.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       a0 = 0.5,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     pwe.logml.pp(
#'       post.samples = d.pp,
#'       bridge.args = list(silent = TRUE),
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'   }
#' }
pwe.logml.pp = function(
    post.samples,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')
  a0        = stan.data$a0

  if ( a0 == 0 ){
    ## PP is equivalent to normal/half-normal prior in this case
    res = pwe.logml.post(
      post.samples = post.samples,
      bridge.args = bridge.args
    )
    return(res)
  }

  ## computing log normalizing constant for PP using all data sets
  res.all = pwe.pp.lognc(
    post.samples   = post.samples,
    is.prior       = FALSE,
    bridge.args    = bridge.args
  )

  ## sample from PP
  pwe_pp_prior = instantiate::stan_package_model(
    name = "pwe_pp_prior",
    package = "hdbayes"
  )

  fit = pwe_pp_prior$sample(data = stan.data,
                            iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                            ...)
  summ = posterior::summarise_draws(fit)

  hist.post.samples = fit$draws(format = 'draws_df')
  attr(x = hist.post.samples, which = 'data') = stan.data

  ## compute log normalizing constant for PP using historical data set
  res.hist = pwe.pp.lognc(
    post.samples   = hist.post.samples,
    is.prior       = TRUE,
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
