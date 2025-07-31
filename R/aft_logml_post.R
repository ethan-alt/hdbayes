#' Log marginal likelihood of an accelerated failure time (AFT) model under a normal/half-normal prior
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of an AFT model under the normal/half-normal prior.
#'
#' @include aft_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [aft.post()] giving posterior samples of an AFT model under the normal/half-normal
#'                          prior, with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"aft_post"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the marginal likelihood of the AFT model under the normal/half-normal prior}
#'  }
#'
#' @references
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1690)
#'     ## take subset for speed purposes
#'     E1690 = E1690[1:100, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'     E1690$cage = as.numeric(scale(E1690$age))
#'     data_list = list(currdata = E1690)
#'     d.post = aft.post(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       beta.sd = 10,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     aft.logml.post(
#'       post.samples = d.post,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
aft.logml.post = function(
    post.samples,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X_obs
  oldnames  = paste0("beta[", 1:p, "]")
  newnames  = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames = c('scale', oldnames)
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on scale
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    beta       = as.numeric( pars[paste0("beta[", 1:data$p,"]")] )
    scale      = as.numeric( pars['scale'] )
    prior_lp   = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
      dnorm(scale, mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc

    eta_obs = data$X_obs %*% beta
    eta_cen = data$X_cen %*% beta
    data_lp = sum( aft_model_lp(data$y_obs, data$y_cen, eta_obs, eta_cen, scale, data$dist) )

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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "aft_post",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
