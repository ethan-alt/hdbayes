#' Log marginal likelihood of a piecewise exponential (PWE) model under a normal/half-normal prior
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of a PWE model under the normal/half-normal prior.
#'
#' @include pwe_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [pwe.post()] giving posterior samples of a PWE model under the normal/half-normal
#'                          prior, with an attribute called 'data' which includes the list of variables specified
#'                          in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"Normal/Half-Normal"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the marginal likelihood of the PWE model under the normal/half-normal prior}
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     d.post = pwe.post(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     pwe.logml.post(
#'       post.samples = d.post,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
pwe.logml.post = function(
    post.samples,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X1        = stan.data$X1
  J         = stan.data$J
  if( p > 0 ){
    oldnames  = c(paste0("beta[", 1:p, "]"), paste0("lambda[", 1:J, "]"))
    newnames  = c(colnames(X1), paste0("basehaz[", 1:J, "]"))
    lb        = c(rep(-Inf, p), rep(0, J))

  }else{
    oldnames  = paste0("lambda[", 1:J, "]")
    newnames  = paste0("basehaz[", 1:J, "]")
    lb        = rep(0, J)
  }
  colnames(d)[colnames(d) %in% newnames] = oldnames
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on baseline hazards
  stan.data$lognc_hazard = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    lambda     = as.numeric( pars[paste0("lambda[", 1:data$J,"]")] )

    if( p > 0 ){
      beta     = as.numeric( pars[paste0("beta[", 1:p,"]")] )
      prior_lp = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
        sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard
      eta      = data$X1 %*% beta
      data_lp  = sum( pwe_lpdf(data$y1, eta, lambda, data$breaks, data$intindx, data$J, data$death_ind) )

    }else{
      prior_lp = sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard
      data_lp  = sum( pwe_lpdf(data$y1, 0, lambda, data$breaks, data$intindx, data$J, data$death_ind) )
    }

    return(data_lp + prior_lp)
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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "Normal/Half-Normal",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
