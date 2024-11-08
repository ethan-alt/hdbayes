#' Estimate the logarithm of the normalizing constant for power prior (PP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the power prior (PP)
#' using all data sets or using historical data set only. Note that the power prior parameter (\eqn{a_0})
#' is treated as fixed.
#'
#' @include aft_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of an AFT model under the power prior (PP) or samples from the PP, with an
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
#'     d.pp = aft.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       a0 = 0.5,
#'       dist = "weibull",
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     aft.pp.lognc(
#'       post.samples = d.pp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
aft.pp.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X0_obs
  oldnames  = paste0("beta[", 1:p, "]")
  newnames  = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames  = c('scale', oldnames)
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on scale
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)
  stan.data$is_prior          = is.prior

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    beta       = as.numeric( pars[paste0("beta[", 1:data$p,"]")] )
    scale      = as.numeric( pars['scale'] )
    prior_lp   = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
      dnorm(scale, mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc

    eta0_obs = data$X0_obs %*% beta
    eta0_cen = data$X0_cen %*% beta
    data_lp  = data$a0 * sum( aft_model_lp(data$y0_obs, data$y0_cen, eta0_obs, eta0_cen, scale, data$dist) )

    if( !data$is_prior ){
      eta_obs = data$X_obs %*% beta
      eta_cen = data$X_cen %*% beta
      data_lp = data_lp + sum( aft_model_lp(data$y_obs, data$y_cen, eta_obs, eta_cen, scale, data$dist) )
    }
    return(data_lp + prior_lp)
  }

  lb           = c(0, rep(-Inf, p))
  ub           = rep(Inf, 1+p)
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
