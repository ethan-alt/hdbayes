#' Estimate the logarithm of the normalizing constant for commensurate prior (CP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the commensurate prior (CP) using
#' historical data sets.
#'
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [glm.commensurate()] giving posterior samples of a GLM under the commensurate
#'                          prior (CP), with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`)
#'                          to pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{lognc}{the estimated logarithm of the normalizing constant}
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
#'   hist.data.list = list(histdata = actg036)
#'   standat = hdbayes:::get.stan.data.cp(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson(),
#'     data.list = hist.data.list
#'   )
#'   glm_commensurate_prior = instantiate::stan_package_model(
#'     name = "glm_commensurate_prior",
#'     package = "hdbayes"
#'   )
#'   fit = glm_commensurate_prior$sample(
#'     data = standat,
#'     iter_warmup = 500, iter_sampling = 1000, chains = 1
#'   )
#'   d.cp.hist  = fit$draws(format = 'draws_df')
#'   attr(x = d.cp.hist, which = 'data') = standat
#'   glm.commensurate.lognc(
#'     post.samples = d.cp.hist
#'   )
#' }
glm.commensurate.lognc = function(
    post.samples,
    bridge.args       = NULL
) {
  ## get Stan data for CP
  stan.data = attr(post.samples, 'data')

  p        = stan.data$p
  K        = stan.data$K
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"), paste0("comm_prec[", 1:p,"]"))
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  d = suppressWarnings(
    as.matrix(post.samples[, oldnames, drop=F])
  )

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_spike = pnorm(0, mean = stan.data$mu_spike, sd = stan.data$sigma_spike, lower.tail = F, log.p = T)
  stan.data$lognc_slab  = pnorm(0, mean = stan.data$mu_slab, sd = stan.data$sigma_slab, lower.tail = F, log.p = T)
  stan.data$lognc_disp  = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

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
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
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
  if( stan.data$dist > 2 ) {
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
        'data' = stan.data,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## Return a list of lognc and output from bridgesampling::bridge_sampler
  res = list(
    'lognc'        = bs$logml,
    'bs'           = bs
  )
  return(res)
}
