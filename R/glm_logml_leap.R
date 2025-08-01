#' Log marginal likelihood of a GLM under latent exchangeability prior (LEAP)
#'
#' @description Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under the latent exchangeability prior (LEAP).
#'
#' @description The arguments related to MCMC sampling are utilized to draw samples from the LEAP. These
#' samples are then used to compute the logarithm of the normalizing constant of the LEAP using historical
#' data sets.
#'
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#' @include glm_leap_lognc.R
#'
#' @export
#'
#' @param post.samples      output from [glm.leap()] giving posterior samples of a GLM under the latent exchangeability
#'                          prior (LEAP), with an attribute called 'data' which includes the list of variables specified
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
#'    \item{model}{"glm_leap"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using
#'    [bridgesampling::bridge_sampler()] to compute the logarithm of the normalizing constant of the
#'    latent exchangeability prior (LEAP) using all data sets}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` containing the output from using
#'    [bridgesampling::bridge_sampler()] to compute the logarithm of the normalizing constant of the
#'    LEAP using historical data sets}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'  }
#'
#' @references
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2024). LEAP: The latent exchangeability prior for borrowing information from historical data. Biometrics, 80(3).
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
#'   formula   = outcome ~ scale(age) + race + treatment + scale(cd4)
#'   family    = binomial('logit')
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.leap    = glm.leap(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.leap(
#'     post.samples = d.leap,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.logml.leap = function(
    post.samples,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  stan.data = attr(post.samples, 'data')
  K         = stan.data$K

  d        = as.matrix(post.samples)
  p        = stan.data$p
  oldnames = paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  oldnames = c(oldnames, "logit_gamma")
  if ( K > 2 ){
    oldnames = c(oldnames, paste0("delta_raw[", 1:(K-2), "]"))
  }
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  stan.data$lognc_disp  = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

  ## compute log normalizing constants for logit(gamma)
  gamma_shape1    = stan.data$conc[1]
  gamma_shape2    = sum(stan.data$conc[2:K])
  stan.data$lognc_logit_gamma = 0

  if( stan.data$gamma_lower != 0 || stan.data$gamma_upper != 1 ) {
    stan.data$lognc_logit_gamma = log( pbeta(stan.data$gamma_upper, shape1 = gamma_shape1, shape2 = gamma_shape2) -
                                   pbeta(stan.data$gamma_lower, shape1 = gamma_shape1, shape2 = gamma_shape2) )
  }

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    betaMat    = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    prior_lp   = sum( dnorm(betaMat, mean = as.numeric(data$mean_beta),
                            sd = as.numeric(data$sd_beta), log = T) )
    betaMat    = matrix(betaMat, nrow = p, ncol = K)

    ## prior on logit(gamma)
    conc         = data$conc
    gamma_shape1 = conc[1]
    gamma_shape2 = sum(conc[2:K])
    logit_gamma  = pars[["logit_gamma"]]
    log1m_gamma  = -log1p_exp(logit_gamma) # log(1 - gamma)
    log_probs    = c(logit_gamma, 0) + log1m_gamma

    prior_lp     = prior_lp + logit_beta_lp(logit_gamma, gamma_shape1, gamma_shape2) -
      data$lognc_logit_gamma

    if( K > 2 ){
      delta_raw = as.numeric(pars[paste0("delta_raw[", 1:(K-2), "]")])
      delta_raw = c(delta_raw, 1 - sum(delta_raw))
      prior_lp  = prior_lp + dirichlet_lp(delta_raw, conc[2:K])
      log_probs = c(logit_gamma, log(delta_raw)) + log1m_gamma
    }

    dist         = data$dist
    link         = data$link
    if ( dist > 2 ) {
      dispersion = pars[paste0("dispersion[", 1:K,"]")]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
    }else {
      dispersion = rep(1.0, K)
    }
    ## historical data likelihood
    prior_lp     = prior_lp  + glm_mixture_lp(data$y0, betaMat, dispersion, log_probs, data$X0, dist, link, data$offs0)
    ## current data likelihood
    data_lp      = glm_lp(data$y, betaMat[, 1], data$X, dist, link, data$offs[,1], dispersion[1])
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p*K)
  ub = rep(Inf, p*K)
  if( stan.data$dist > 2 ) {
    lb = c(lb, rep(0, K) )
    ub = c(ub, rep(Inf, K) )
  }
  lb = c(lb, binomial('logit')$linkfun(stan.data$gamma_lower))
  ub = c(ub, binomial('logit')$linkfun(stan.data$gamma_upper))
  if( K > 2 ){
    lb = c(lb, rep(0, K-2))
    ub = c(ub, rep(1, K-2))
  }
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

  ## get Stan data for LEAP using historical data sets
  hist.stan.data = stan.data[c('n0', 'p', 'K', 'y0', 'X0', 'mean_beta',
                                'sd_beta', 'disp_mean', 'disp_sd', 'conc',
                                'gamma_lower', 'gamma_upper', 'dist',
                                'link', 'offs0')]

  ## sample from LEAP using historical data sets
  glm_leap_prior = instantiate::stan_package_model(
    name = "glm_leap_prior",
    package = "hdbayes"
  )
  fit = glm_leap_prior$sample(data = hist.stan.data,
                              iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                              ...)
  summ = posterior::summarise_draws(fit)

  hist.post.samples = fit$draws(format = 'draws_df')
  attr(x = hist.post.samples, which = 'data') = hist.stan.data

  ## compute log normalizing constant for LEAP using historical data sets
  res.hist = glm.leap.lognc(
    post.samples   = hist.post.samples,
    bridge.args    = bridge.args
  )

  ## Return a list of model name, estimated log marginal likelihood, outputs from bridgesampling::bridge_sampler,
  ## the minimum estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat
  res = list(
    'model'        = "glm_leap",
    'logml'        = bs$logml - res.hist$lognc,
    'bs'           = bs,
    'bs.hist'      = res.hist$bs,
    'min_ess_bulk' = min(summ[, 'ess_bulk'], na.rm = T),
    'max_Rhat'     = max(summ[, 'rhat'], na.rm = T)
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
