#' Estimate the logarithm of the normalizing constant for power prior (PP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the power prior (PP)
#' using all data sets or using historical data set only. Note that the power prior parameter (\eqn{a_0})
#' is treated as fixed.
#'
#' @include pwe_loglik.R
#' @include mixture_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a CurePWE model under the power prior (PP) or samples from the PP, with an
#'                          attribute called 'data' which includes the list of variables specified in the data block
#'                          of the Stan program.
#' @param is.prior          whether the samples are from the PP (using historical data set only). Defaults to FALSE.
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
#'     d.pp = curepwe.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       a0 = 0.5,
#'       breaks = breaks,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     curepwe.pp.lognc(
#'       post.samples = d.pp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
curepwe.pp.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X0        = stan.data$X0
  J         = stan.data$J
  oldnames  = c(paste0("beta[", 1:p, "]"), paste0("lambda[", 1:J, "]"))
  newnames  = c(colnames(X0), paste0("basehaz[", 1:J, "]"))
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames  = c(oldnames, "logit_p_cured")
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on baseline hazards
  stan.data$lognc_hazard = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )
  stan.data$is_prior     = is.prior

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    beta          = as.numeric( pars[paste0("beta[", 1:data$p,"]")] )
    lambda        = as.numeric( pars[paste0("lambda[", 1:data$J,"]")] )
    logit_p_cured = as.numeric( pars[["logit_p_cured"]] )
    log1m_p_cured = -log1p_exp(logit_p_cured) # log(1 - p_cured)
    log_p_cured   = logit_p_cured + log1m_p_cured # log(p_cured)

    prior_lp   = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
      sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard +
      dnorm(logit_p_cured, mean = data$logit_p_cured_mean, sd = data$logit_p_cured_sd, log = T)

    eta0      = data$X0 %*% beta
    contribs0 = cbind(log_p_cured + log(1 - data$death_ind0),
                      log1m_p_cured + pwe_lpdf(data$y0, eta0, lambda, data$breaks, data$intindx0, data$J, data$death_ind0))
    data_lp   = apply(contribs0, 1, log_sum_exp)
    data_lp   = data$a0 * sum( data_lp )

    if( !data$is_prior ){
      eta      = data$X1 %*% beta
      contribs = cbind(log_p_cured + log(1 - data$death_ind),
                       log1m_p_cured + pwe_lpdf(data$y1, eta, lambda, data$breaks, data$intindx, data$J, data$death_ind))
      data_lp  = data_lp + sum( apply(contribs, 1, log_sum_exp) )
    }
    return(data_lp + prior_lp)
  }

  lb           = c(rep(-Inf, p), rep(0, J), -Inf)
  ub           = rep(Inf, p+J+1)
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
