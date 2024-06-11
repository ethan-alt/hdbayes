#' Estimate the logarithm of the normalizing constant for commensurate prior (CP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the normalizing
#' constant for the CP using historical data sets.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param hist.data.list    a list of `data.frame`s giving historical data sets.
#' @param hist.offset.list  a list of vectors giving the offsets for each historical data. The length of hist.offset.list
#'                          is equal to the length of hist.data.list. The length of each element of offset.list is equal to
#'                          the number of rows in the corresponding element of hist.data.list. Defaults to a list of vectors
#'                          of 0s.
#' @param beta0.mean        a scalar or a vector whose dimension is equal to the number of regression coefficients
#'                          giving the mean parameters for the prior on the historical data regression coefficients. If a
#'                          scalar is provided, same as for tau. Defaults to a vector of 0s.
#' @param beta0.sd          a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the prior on the historical data regression coefficients. If a scalar is
#'                          provided, same as for tau. Defaults to a vector of 10s.
#' @param hist.disp.mean    a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          means for the half-normal priors on the dispersion parameters. If a scalar is provided, same as
#'                          for tau. Defaults to a vector of 0s.
#' @param hist.disp.sd      a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          sds for the half-normal priors on the dispersion parameters. If a scalar is provided, same as
#'                          for tau. Defaults to a vector of 10s.
#' @param p.spike           a scalar between 0 and 1 giving the probability of the spike component in spike-and-slab prior
#'                          on commensurability parameter \eqn{\tau}. Defaults to 0.1.
#' @param spike.mean        a scalar giving the location parameter for the half-normal prior (spike component) on \eqn{\tau}.
#'                          Defaults to 200.
#' @param spike.sd          a scalar giving the scale parameter for the half-normal prior (spike component) on \eqn{\tau}.
#'                          Defaults to 0.1.
#' @param slab.mean         a scalar giving the location parameter for the half-normal prior (slab component) on \eqn{\tau}.
#'                          Defaults to 0.
#' @param slab.sd           a scalar giving the scale parameter for the half-normal prior (slab component) on \eqn{\tau}.
#'                          Defaults to 5.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{lognc}{the estimated logarithm of the normalizing constant}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from [bridgesampling::bridge_sampler()]}
#'  }
#'
#' @references
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg036 = actg036[1:50, ]
#'   glm.commensurate.lognc(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson(),
#'     hist.data.list = list(histdata = actg036),
#'     p.spike = 0.1,
#'     chains = 1, iter_warmup = 500, iter_sampling = 2000
#'   )
#' }
glm.commensurate.lognc = function(
    formula,
    family,
    hist.data.list,
    hist.offset.list  = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    hist.disp.mean    = NULL,
    hist.disp.sd      = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for CP
  standat = get.stan.data.cp(
    formula        = formula,
    family         = family,
    data.list      = hist.data.list,
    offset.list    = hist.offset.list,
    beta0.mean     = beta0.mean,
    beta0.sd       = beta0.sd,
    disp.mean      = hist.disp.mean,
    disp.sd        = hist.disp.sd,
    p.spike        = p.spike,
    spike.mean     = spike.mean,
    spike.sd       = spike.sd,
    slab.mean      = slab.mean,
    slab.sd        = slab.sd
  )

  glm_commensurate_prior = instantiate::stan_package_model(
    name = "glm_commensurate_prior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit      = glm_commensurate_prior$sample(data = standat,
                                 iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                 ...)
  d        = fit$draws(format = 'draws_df')
  summ     = posterior::summarise_draws(d)

  p        = standat$p
  K        = standat$K
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"), paste0("comm_prec[", 1:p,"]"))
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  d = suppressWarnings(
    as.matrix(d[, oldnames, drop=F])
  )

  standat$lognc_spike = pnorm(0, mean = standat$mu_spike, sd = standat$sigma_spike, lower.tail = F, log.p = T)
  standat$lognc_slab  = pnorm(0, mean = standat$mu_slab, sd = standat$sigma_slab, lower.tail = F, log.p = T)

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
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
    if ( dist > 2 ) {
      dispersion = pars[paste0("dispersion[", 1:K,"]")]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) -
        sum( pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T) )
      hist_lp    = sum( sapply(1:data$K, function(k){
        start.idx = data$start_idx[k]
        end.idx   = data$end_idx[k]
        y         = data$y[ start.idx:end.idx ]
        X         = data$X[ start.idx:end.idx, ]
        offs      = data$offs[ start.idx:end.idx ]
        glm_lp(y, beta0, X, dist, link, offs, dispersion[k])
      }) )
    }else {
      dispersion = 1.0
      hist_lp    = glm_lp(data$y, beta0, data$X, dist, link, data$offs, dispersion)
    }
    return(hist_lp + prior_lp)
  }

  lb = c(rep(-Inf, p*2), rep(0, p))
  ub = rep(Inf, p*3)
  if( standat$dist > 2 ) {
    lb = c(lb, rep(0, K) )
    ub = c(ub, rep(Inf, K) )
  }
  names(ub) = colnames(d)
  names(lb) = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = d,
        'log_posterior' = log_density,
        'data' = standat,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## Return a list of lognc, min_n_eff, max_Rhat, and output from bridgesampling::bridge_sampler
  res = list(
    'lognc'        = bs$logml,
    'min_ess_bulk' = min(summ[, 'ess_bulk']),
    'max_Rhat'     = max(summ[, 'rhat']),
    'bs'           = bs
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
