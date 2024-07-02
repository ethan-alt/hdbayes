#' Estimate the logarithm of the normalizing constant for power prior (PP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the power prior (PP)
#' using all data sets or using historical data sets only. Note that the power prior parameters (\eqn{a_0}'s)
#' are treated as fixed.
#'
#' @include expfam_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a GLM under the power prior (PP) or samples from the PP, with an
#'                          attribute called 'data' which includes the list of variables specified in the data block
#'                          of the Stan program.
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
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.pp = glm.pp(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data_list,
#'     a0.vals = 0.5,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.pp.lognc(
#'     post.samples = d.pp,
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.pp.lognc = function(
    post.samples,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = stan.data$p
  X        = stan.data$X
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, 'dispersion')
  }
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  stan.data$lognc_disp  = pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T)

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    beta       = pars[paste0("beta[", 1:data$p,"]")]
    prior_lp   = sum( dnorm(beta, mean = data$mean_beta, sd = data$sd_beta, log = T) )
    dist       = data$dist
    link       = data$link
    dispersion = 1.0
    if ( dist > 2 ){
      dispersion = pars[["dispersion"]]
      prior_lp   = prior_lp +
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) - data$lognc_disp
    }
    data_lp = as.numeric( data$a0_vals %*% sapply(1:data$K, function(k){
      start.idx = data$start_idx[k]
      end.idx   = data$end_idx[k]
      y         = data$y[ start.idx:end.idx ]
      X         = data$X[ start.idx:end.idx, ]
      offs      = data$offs[ start.idx:end.idx ]
      glm_lp(y, beta, X, dist, link, offs, dispersion)
    }) )
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p)
  ub = rep(Inf, p)
  if( stan.data$dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
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
