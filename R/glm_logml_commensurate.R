#' Log marginal likelihood of a GLM under commensurate prior (CP)
#'
#' @description Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under the commensurate prior (CP).
#'
#' @description The arguments related to MCMC sampling are utilized to draw samples from the commensurate prior.
#' These samples are then used to compute the logarithm of the normalizing constant of the commensurate prior using
#' historical data sets.
#'
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#' @include glm_commensurate_lognc.R
#'
#' @export
#'
#' @param post.samples      output from [glm.commensurate()] giving posterior samples of a GLM under the commensurate
#'                          prior (CP), with an attribute called 'data' which includes the list of variables specified
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
#'    \item{model}{"glm_commensurate"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the normalizing constant of the commensurate prior (CP) using all data sets}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` containing the output from using
#'    [bridgesampling::bridge_sampler()] to compute the logarithm of the normalizing constant of the CP using historical
#'    data sets}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'  }
#'
#' @references
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019   = actg019[1:100, ]
#'   actg036   = actg036[1:50, ]
#'   formula   = cd4 ~ treatment + age + race
#'   family    = poisson()
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.cp      = glm.commensurate(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     p.spike = 0.1,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.commensurate(
#'     post.samples = d.cp,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.logml.commensurate = function(
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

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = stan.data$p
  X        = stan.data$X
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"))
  newnames = c(colnames(X), paste0( colnames(X), '_hist') )
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
  }
  oldnames = c(oldnames, paste0("comm_prec[", 1:p,"]"))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_spike = pnorm(0, mean = stan.data$mu_spike, sd = stan.data$sigma_spike, lower.tail = F, log.p = T)
  stan.data$lognc_slab  = pnorm(0, mean = stan.data$mu_slab, sd = stan.data$sigma_slab, lower.tail = F, log.p = T)
  stan.data$lognc_disp  = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    N          = data$N
    beta       = pars[paste0("beta[", 1:p,"]")]
    beta0      = pars[paste0("beta0[", 1:p,"]")]
    comm_prec  = pars[paste0("comm_prec[", 1:p,"]")]
    comm_sd    = 1/sqrt(comm_prec)

    ## prior on beta0 and beta
    prior_lp   = sum( dnorm(beta0, mean = data$beta0_mean, sd = data$beta0_sd, log = T) ) +
      sum( dnorm(beta, mean = beta0, sd = comm_sd, log = T) )

    ## spike and slab prior on commensurability
    prior_lp   = prior_lp + sum( sapply(1:p, function(i){
      p_spike    = data$p_spike
      spike_lp   = dnorm(comm_prec[i], mean = data$mu_spike, sd = data$sigma_spike, log = T) - data$lognc_spike
      slab_lp    = dnorm(comm_prec[i], mean = data$mu_slab, sd = data$sigma_slab, log = T) - data$lognc_slab
      log_sum_exp( c( log(p_spike) + spike_lp, log(1 - p_spike) + slab_lp ) )
    }) )

    dist       = data$dist
    link       = data$link
    start.idx  = data$start_idx
    end.idx    = data$end_idx
    dispersion = rep(1.0, K)
    if ( dist > 2 ){
      ## prior on dispersion
      dispersion = pars[c("dispersion", paste0( "dispersion", "_hist_", 1:(K-1) ))]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
      ## historical data likelihood
      prior_lp    = prior_lp + sum( sapply(2:K, function(k){
        y         = data$y[ start.idx[k]:end.idx[k] ]
        X         = data$X[ start.idx[k]:end.idx[k], ]
        offs      = data$offs[ start.idx[k]:end.idx[k] ]
        glm_lp(y, beta0, X, dist, link, offs, dispersion[k])
      }) )
    }else {
      ## historical data likelihood
      prior_lp    = prior_lp +
        glm_lp(data$y[ start.idx[2]:N ], beta0,
               data$X[ start.idx[2]:N, ], dist, link,
               data$offs[ start.idx[2]:N ], 1.0)
    }
    ## current data likelihood
    y         = data$y[ start.idx[1]:end.idx[1] ]
    X         = data$X[ start.idx[1]:end.idx[1], ]
    offs      = data$offs[ start.idx[1]:end.idx[1] ]
    data_lp   = glm_lp(y, beta, X, dist, link, offs, dispersion[1])
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p*2)
  ub = rep(Inf, p*2)
  if( stan.data$dist > 2 ) {
    lb = c(lb, rep(0, K))
    ub = c(ub, rep(Inf, K))
  }
  lb = c(lb, rep(0, p))
  ub = c(ub, rep(Inf, p))
  names(ub) = colnames(d)
  names(lb) = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = d,
        'log_posterior' = log_density,
        'data' = stan.data,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## get Stan data for CP using historical data sets
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

  ## sample from CP using historical data sets
  glm_commensurate_prior = instantiate::stan_package_model(
    name = "glm_commensurate_prior",
    package = "hdbayes"
  )
  fit = glm_commensurate_prior$sample(data = hist.stan.data,
                                iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                ...)
  summ = posterior::summarise_draws(fit)

  hist.post.samples = fit$draws(format = 'draws_df')
  attr(x = hist.post.samples, which = 'data') = hist.stan.data

  ## compute log normalizing constant for CP using historical data sets
  res.hist = glm.commensurate.lognc(
    post.samples   = hist.post.samples,
    bridge.args    = bridge.args
  )

  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## the minimum estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat
  res = list(
    'model'        = "glm_commensurate",
    'logml'        = bs$logml - res.hist$lognc,
    'bs'           = bs,
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
