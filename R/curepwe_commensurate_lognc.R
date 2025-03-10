#' Estimate the logarithm of the normalizing constant for commensurate prior (CP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the commensurate prior (CP) using
#' all data sets or using historical data set only.
#'
#' @include pwe_loglik.R
#' @include mixture_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a CurePWE model under the commensurate prior (CP) or samples from
#'                          the CP, with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     d.cp = curepwe.commensurate(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       p.spike = 0.1,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     curepwe.commensurate.lognc(
#'       post.samples = d.cp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
curepwe.commensurate.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  ## get Stan data for CP
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)

  ## rename parameters
  p        = stan.data$p
  X0       = stan.data$X0
  J        = stan.data$J

  if( !is.prior ){
    oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"),
                 paste0("lambda[", 1:J, "]"), paste0("lambda0[", 1:J, "]"))
    newnames = c(colnames(X0), paste0(colnames(X0), '_hist'),
                 paste0("basehaz[", 1:J, "]"), paste0("basehaz_hist[", 1:J, "]"))
    colnames(d)[colnames(d) %in% newnames] = oldnames
    oldnames = c(oldnames, "logit_p_cured", "logit_p_cured0")
  }else{
    oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"),
                 paste0("lambda0[", 1:J, "]"), "logit_p_cured0")
  }
  oldnames = c(oldnames, paste0("comm_prec[", 1:p,"]"))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_hazard      = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )
  stan.data$lognc_spike       = pnorm(0, mean = stan.data$mu_spike, sd = stan.data$sigma_spike, lower.tail = F, log.p = T)
  stan.data$lognc_slab        = pnorm(0, mean = stan.data$mu_slab, sd = stan.data$sigma_slab, lower.tail = F, log.p = T)
  stan.data$is_prior          = is.prior

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p          = data$p
    beta       = as.numeric( pars[paste0("beta[", 1:p,"]")] )
    beta0      = as.numeric( pars[paste0("beta0[", 1:p,"]")] )
    lambda0    = as.numeric( pars[paste0("lambda0[", 1:data$J,"]")] )
    comm_prec  = as.numeric( pars[paste0("comm_prec[", 1:p,"]")] )
    comm_sd    = 1/sqrt(comm_prec)

    logit_p_cured0 = as.numeric( pars[["logit_p_cured0"]] )
    log1m_p_cured0 = -log1p_exp(logit_p_cured0) # log(1 - p_cured0)
    log_p_cured0   = logit_p_cured0 + log1m_p_cured0 # log(p_cured0)

    ## prior on beta0 and beta
    prior_lp   = sum( dnorm(beta0, mean = data$beta0_mean, sd = data$beta0_sd, log = T) ) +
      sum( dnorm(beta, mean = beta0, sd = comm_sd, log = T) )

    ## priors on lambda0 and logit_p_cured0
    prior_lp   = prior_lp + sum( dnorm(lambda0, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard +
      dnorm(logit_p_cured0, mean = data$logit_p_cured_mean, sd = data$logit_p_cured_sd, log = T)

    ## spike and slab prior on commensurability
    prior_lp   = prior_lp + sum( sapply(1:p, function(i){
      p_spike    = data$p_spike
      spike_lp   = dnorm(comm_prec[i], mean = data$mu_spike, sd = data$sigma_spike, log = T) - data$lognc_spike
      slab_lp    = dnorm(comm_prec[i], mean = data$mu_slab, sd = data$sigma_slab, log = T) - data$lognc_slab
      log_sum_exp( c( log(p_spike) + spike_lp, log1p(-p_spike) + slab_lp ) )
    }) )

    eta0       = data$X0 %*% beta0
    contribs0  = cbind(log_p_cured0 + log(1 - data$death_ind0),
                       log1m_p_cured0 + pwe_lpdf(data$y0, eta0, lambda0, data$breaks, data$intindx0, data$J, data$death_ind0))
    data_lp    = sum( apply(contribs0, 1, log_sum_exp) )

    if ( !data$is_prior ) {
      lambda     = as.numeric( pars[paste0("lambda[", 1:data$J,"]")] )

      logit_p_cured = as.numeric( pars[["logit_p_cured"]] )
      log1m_p_cured = -log1p_exp(logit_p_cured) # log(1 - p_cured)
      log_p_cured   = logit_p_cured + log1m_p_cured # log(p_cured)

      ## priors on lambda and logit_p_cured
      prior_lp   = prior_lp + sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard +
        dnorm(logit_p_cured, mean = data$logit_p_cured_mean, sd = data$logit_p_cured_sd, log = T)

      eta        = data$X1 %*% beta
      contribs   = cbind(log_p_cured + log(1 - data$death_ind),
                         log1m_p_cured + pwe_lpdf(data$y1, eta, lambda, data$breaks, data$intindx, data$J, data$death_ind))
      data_lp    = data_lp + sum( apply(contribs, 1, log_sum_exp) )
    }
    return(data_lp + prior_lp)
  }

  if( !is.prior ){
    lb = c(rep(-Inf, p*2), rep(0, J*2), rep(-Inf, 2), rep(0, p))
  }else{
    lb = c(rep(-Inf, p*2), rep(0, J), -Inf, rep(0, p))
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
