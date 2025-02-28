#' Estimate the logarithm of the normalizing constant for Bayesian hierarchical model (BHM)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for Bayesian hierarchical
#' model (BHM) using all data sets or using historical data set only.
#'
#' @include pwe_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a PWE model under the Bayesian hierarchical model (BHM) or samples from the
#'                          prior induced by the BHM, with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param is.prior          whether the samples are from the prior induced by the BHM (using historical data set only).
#'                          Defaults to FALSE.
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
#'     d.bhm = pwe.bhm(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     pwe.bhm.lognc(
#'       post.samples = d.bhm,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
pwe.bhm.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  ## get Stan data for BHM
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)

  ## rename parameters
  p        = stan.data$p
  J        = stan.data$J

  if( is.prior ){
    oldnames = c(paste0("lambda0[", 1:J, "]"), paste0('beta0_raw[', 1:p, ']'))
  }else{
    oldnames = c(paste0("lambda[", 1:J, "]"), paste0("lambda0[", 1:J, "]"))
    newnames = c(paste0("basehaz[", 1:J, "]"), paste0("basehaz_hist[", 1:J, "]"))
    colnames(d)[colnames(d) %in% newnames] = oldnames

    oldnames = c(paste0("lambda[", 1:J, "]"), paste0("lambda0[", 1:J, "]"),
                 paste0('beta_raw[', 1:p, ']'), paste0('beta0_raw[', 1:p, ']'))
  }
  oldnames = c(oldnames, paste0('beta_mean[', 1:p, ']'), paste0('beta_sd[', 1:p, ']'))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_hazard      = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )
  stan.data$lognc_beta_sd     = sum( pnorm(0, mean = stan.data$meta_sd_mean, sd = stan.data$meta_sd_sd, lower.tail = F, log.p = T) )
  stan.data$is_prior          = is.prior

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    beta0_raw  = as.numeric( pars[paste0('beta0_raw[', 1:p, ']')] )
    beta_mean  = as.numeric( pars[paste0('beta_mean[', 1:p, ']')] )
    beta_sd    = as.numeric( pars[paste0('beta_sd[', 1:p, ']')] )
    lambda0    = as.numeric( pars[paste0("lambda0[", 1:data$J,"]")] )

    ## prior on beta_mean and beta_sd
    prior_lp   = sum( dnorm(beta_mean, mean = data$meta_mean_mean, sd = data$meta_mean_sd, log = T) ) +
      sum( dnorm(beta_sd, mean = data$meta_sd_mean, sd = data$meta_sd_sd, log = T) ) - data$lognc_beta_sd

    ## prior on beta0_raw (equivalent to prior on beta0)
    prior_lp   = prior_lp + sum( dnorm(beta0_raw, mean = 0, sd = 1, log = T) )

    ## prior on lambda0
    prior_lp   = prior_lp + sum( dnorm(lambda0, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard

    beta0      = beta_mean + beta0_raw * beta_sd
    eta0       = data$X0 %*% beta0
    data_lp    = sum( pwe_lpdf(data$y0, eta0, lambda0, data$breaks, data$intindx0, data$J, data$death_ind0) )

    if( !data$is_prior ){
      beta_raw   = as.numeric( pars[paste0('beta_raw[', 1:p, ']')] )
      lambda     = as.numeric( pars[paste0("lambda[", 1:data$J,"]")] )

      ## prior on beta_raw (equivalent to prior on beta)
      prior_lp   = prior_lp + sum( dnorm(beta_raw, mean = 0, sd = 1, log = T) )

      ## prior on lambda
      prior_lp   = prior_lp + sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard

      beta       = beta_mean + beta_raw * beta_sd
      eta        = data$X1 %*% beta
      data_lp    = data_lp + sum( pwe_lpdf(data$y1, eta, lambda, data$breaks, data$intindx, data$J, data$death_ind) )
    }
    return(data_lp + prior_lp)
  }

  if( is.prior ){
    lb           = c(rep(0, J), rep(-Inf, 2*p), rep(0, p))
  }else{
    lb           = c(rep(0, 2*J), rep(-Inf, 3*p), rep(0, p))
  }
  ub           = rep(Inf, length(lb))
  names(ub)    = colnames(d)
  names(lb)    = names(ub)

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
