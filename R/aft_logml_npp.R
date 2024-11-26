#' Log marginal likelihood of an accelerated failure time (AFT) model under normalized power prior (NPP)
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of an AFT model under the
#' normalized power prior (NPP).
#'
#' @include aft_loglik.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [aft.npp()] giving posterior samples of an AFT model under the normalized
#'                          power prior (NPP), with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"NPP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the marginal likelihood of the normalized power prior (NPP)}
#'  }
#'
#' @references
#'  Duan, Y., Ye, K., and Smith, E. P. (2005). Evaluating water quality using power priors to incorporate historical information. Environmetrics, 17(1), 95–106.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' \donttest{
#'   if(requireNamespace("parallel")){
#'     library(parallel)
#'     ncores    = 2
#'
#'     if(requireNamespace("survival")){
#'       library(survival)
#'       data(E1684)
#'       data(E1690)
#'       ## take subset for speed purposes
#'       E1684 = E1684[1:100, ]
#'       E1690 = E1690[1:50, ]
#'       ## replace 0 failure times with 0.50 days
#'       E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'       E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'       E1684$cage = as.numeric(scale(E1684$age))
#'       E1690$cage = as.numeric(scale(E1690$age))
#'       data_list = list(currdata = E1690, histdata = E1684)
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
#'     }
#'
#'     a0 = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## call created function
#'       ## wrapper to obtain log normalizing constant in parallel package
#'       logncfun = function(a0, ...){
#'         hdbayes::aft.npp.lognc(
#'           formula = formula, histdata = data_list[[2]], a0 = a0, dist = "weibull",
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'data_list'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from normalized power prior
#'       d.npp = aft.npp(
#'         formula = formula,
#'         data.list = data_list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = a0.lognc$lognc,
#'         dist = "weibull",
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'       aft.logml.npp(
#'         post.samples = d.npp,
#'         bridge.args = list(silent = TRUE)
#'       )
#'     }
#'   }
#' }
aft.logml.npp = function(
    post.samples,
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
  oldnames  = c(oldnames, "scale", "logit_a0")
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on scale
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)

  ## compute log normalizing constant for a0
  a0_shape1       = stan.data$a0_shape1
  a0_shape2       = stan.data$a0_shape2

  stan.data$lognc_logit_a0 = 0
  if( stan.data$a0_lower != 0 || stan.data$a0_upper != 1 ) {
    stan.data$lognc_logit_a0 = log( pbeta(stan.data$a0_upper, shape1 = a0_shape1, shape2 = a0_shape2) -
                                      pbeta(stan.data$a0_lower, shape1 = a0_shape1, shape2 = a0_shape2) )
  }

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    a0_shape1  = data$a0_shape1
    a0_shape2  = data$a0_shape2
    a0_lower   = data$a0_lower
    a0_upper   = data$a0_upper

    beta       = as.numeric( pars[paste0("beta[", 1:data$p,"]")] )
    scale      = as.numeric( pars['scale'] )
    prior_lp   = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
      dnorm(scale, mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc

    logit_a0   = as.numeric(pars["logit_a0"])
    a0         = binomial('logit')$linkinv(logit_a0)
    ## prior on logit(a0)
    prior_lp   = prior_lp + logit_beta_lp(logit_a0, shape1 = a0_shape1, shape2 = a0_shape2) - data$lognc_logit_a0

    eta_obs  = data$X_obs %*% beta
    eta_cen  = data$X_cen %*% beta
    eta0_obs = data$X0_obs %*% beta
    eta0_cen = data$X0_cen %*% beta
    data_lp  = a0 * aft_model_lp(data$y0_obs, data$y0_cen, eta0_obs, eta0_cen, scale, data$dist) +
      aft_model_lp(data$y_obs, data$y_cen, eta_obs, eta_cen, scale, data$dist)

    ## subtract log nc from power prior
    a0_lognc = data$a0_lognc
    lognc    = data$lognc
    prior_lp = prior_lp - pp_lognc(a0, a0_lognc, lognc)
    return(data_lp + prior_lp)
  }

  lb        = c(rep(-Inf, p), 0)
  ub        = rep(Inf, length(lb))
  lb        = c(lb, binomial('logit')$linkfun(stan.data$a0_lower))
  ub        = c(ub, binomial('logit')$linkfun(stan.data$a0_upper))
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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "NPP",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
