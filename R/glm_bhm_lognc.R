#' Estimate the logarithm of the normalizing constant for Bayesian hierarchical model (BHM)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for Bayesian hierarchical
#' model (BHM) using all data sets or using historical data sets only.
#'
#' @include expfam_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a GLM under the Bayesian hierarchical model (BHM) or samples from the
#'                          prior induced by the BHM, with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
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
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   formula = outcome ~ scale(age) + race + treatment + scale(cd4)
#'   family = binomial('logit')
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.bhm = glm.bhm(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.bhm.lognc(
#'     post.samples = d.bhm,
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.bhm.lognc = function(
    post.samples,
    bridge.args       = NULL
) {
  ## get Stan data for BHM
  stan.data = attr(post.samples, 'data')

  ## rename parameters
  p        = stan.data$p
  K        = stan.data$K
  oldnames = paste0("beta_raw[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  oldnames = c(oldnames, paste0("beta_mean[", 1:p, "]"), paste0("beta_sd[", 1:p, "]"))
  if ( stan.data$dist > 2 ) {
    if (K == 1) {
      oldnames = c(oldnames, 'dispersion')
    }else {
      oldnames = c(oldnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
    }
  }
  d    = suppressWarnings(
    as.matrix( post.samples[, oldnames, drop=F] )
  )

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_beta_sd = sum( pnorm(0, mean = stan.data$meta_sd_mean, sd = stan.data$meta_sd_sd, lower.tail = F, log.p = T) )
  stan.data$lognc_disp    = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta_raw   = pars[paste0("beta_raw[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta_mean  = pars[paste0("beta_mean[", 1:p, "]")]
    beta_sd    = pars[paste0("beta_sd[", 1:p, "]")]
    ## prior on beta_mean and beta_sd
    prior_lp   = sum( dnorm(beta_mean, mean = data$meta_mean_mean, sd = data$meta_mean_sd, log = T) ) +
      sum( dnorm(beta_sd, mean = data$meta_sd_mean, sd = data$meta_sd_sd, log = T) ) - data$lognc_beta_sd
    ## prior on beta_raw (equivalent to prior on beta)
    prior_lp   = prior_lp + sum( dnorm(beta_raw, mean = 0, sd = 1, log = T) )
    beta_raw   = matrix(beta_raw, nrow = p, ncol = K)
    beta       = apply(beta_raw, 2, function(c){ as.numeric(beta_mean + beta_sd * c) })
    dist       = data$dist
    link       = data$link
    dispersion = rep(1.0, K)
    if ( dist > 2 ){
      ## prior on dispersion
      if (K == 1) {
        dispersion = pars[["dispersion"]]
      }else {
        dispersion = pars[c("dispersion", paste0( "dispersion", "_hist_", 1:(K-1) ))]
      }
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
    }
    data_lp = sum(
      sapply(1:K, function(k){
        start.idx = data$start_idx[k]
        end.idx   = data$end_idx[k]
        y         = data$y[ start.idx:end.idx ]
        X         = data$X[ start.idx:end.idx, ]
        offs      = data$offs[ start.idx:end.idx ]
        glm_lp(y, beta[, k], X, dist, link, offs, dispersion[k])
      })
    )
    return(data_lp + prior_lp)
  }

  lb           = c(rep(-Inf, p*(K+1)), rep(0, p))
  ub           = rep(Inf, p*(K+2))
  if( stan.data$dist > 2 ) {
    lb = c(lb, rep(0, K))
    ub = c(ub, rep(Inf, K))
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
