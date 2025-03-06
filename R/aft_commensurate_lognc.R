#' Estimate the logarithm of the normalizing constant for commensurate prior (CP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the commensurate prior (CP) using
#' all data sets or using historical data set only.
#'
#' @include aft_loglik.R
#' @include mixture_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of an AFT model under the commensurate prior (CP) or samples from the CP,
#'                          with an attribute called 'data' which includes the list of variables specified in the data
#'                          block of the Stan program.
#' @param is.prior          whether the samples are from the CP (using historical data set only). Defaults to FALSE.
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
#'     d.cp = aft.commensurate(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       p.spike = 0.1,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     aft.commensurate.lognc(
#'       post.samples = d.cp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
aft.commensurate.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  ## get Stan data for CP
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)

  ## rename parameters
  p        = stan.data$p
  X        = stan.data$X0_obs
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  if( !is.prior ){
    oldnames = c(oldnames, "scale")
    newnames = c(newnames, "scale")
  }
  oldnames = c(oldnames, paste0("beta0[", 1:p, "]"), "scale0")
  newnames = c(newnames, paste0( colnames(X), '_hist'), "scale_hist")
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames = c(oldnames, paste0("comm_prec[", 1:p,"]"))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants for half-normal priors
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)
  stan.data$lognc_spike       = pnorm(0, mean = stan.data$mu_spike, sd = stan.data$sigma_spike, lower.tail = F, log.p = T)
  stan.data$lognc_slab        = pnorm(0, mean = stan.data$mu_slab, sd = stan.data$sigma_slab, lower.tail = F, log.p = T)
  stan.data$is_prior          = is.prior

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p          = data$p
    beta       = as.numeric( pars[paste0("beta[", 1:p,"]")] )
    beta0      = as.numeric( pars[paste0("beta0[", 1:p,"]")] )
    scale0     = as.numeric( pars['scale0'] )
    comm_prec  = as.numeric( pars[paste0("comm_prec[", 1:p,"]")] )
    comm_sd    = 1/sqrt(comm_prec)

    ## prior on beta0 and beta
    prior_lp   = sum( dnorm(beta0, mean = data$beta0_mean, sd = data$beta0_sd, log = T) ) +
      sum( dnorm(beta, mean = beta0, sd = comm_sd, log = T) )

    ## prior on scale0
    prior_lp   = prior_lp + dnorm(scale0, mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc

    ## spike and slab prior on commensurability
    prior_lp   = prior_lp + sum( sapply(1:p, function(i){
      p_spike    = data$p_spike
      spike_lp   = dnorm(comm_prec[i], mean = data$mu_spike, sd = data$sigma_spike, log = T) - data$lognc_spike
      slab_lp    = dnorm(comm_prec[i], mean = data$mu_slab, sd = data$sigma_slab, log = T) - data$lognc_slab
      log_sum_exp( c( log(p_spike) + spike_lp, log1p(-p_spike) + slab_lp ) )
    }) )

    eta0_obs   = data$X0_obs %*% beta0
    eta0_cen   = data$X0_cen %*% beta0
    data_lp    = aft_model_lp(data$y0_obs, data$y0_cen, eta0_obs, eta0_cen, scale0, data$dist)

    if ( !data$is_prior ) {
      scale      = as.numeric( pars['scale'] )
      ## prior on scale
      prior_lp   = prior_lp + dnorm(scale, mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc

      eta_obs    = data$X_obs %*% beta
      eta_cen    = data$X_cen %*% beta
      data_lp    = data_lp + aft_model_lp(data$y_obs, data$y_cen, eta_obs, eta_cen, scale, data$dist)
    }
    return(data_lp + prior_lp)
  }

  if( !is.prior ){
    lb = c(rep(-Inf, p), 0, rep(-Inf, p), rep(0, p+1))
  }else{
    lb = c(rep(-Inf, 2*p), rep(0, p+1))
  }
  ub        = rep(Inf, length(lb))
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
