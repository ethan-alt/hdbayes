#' Log marginal likelihood of a GLM under normalized power prior (NPP)
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of a GLM under the
#' normalized power prior (NPP).
#'
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [glm.npp()] giving posterior samples of a GLM under the normalized power
#'                          prior (NPP), with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass onto
#'                          [bridgesampling::bridge_sampler()].
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
#'     data(actg019)
#'     data(actg036)
#'     ## take subset for speed purposes
#'     actg019 = actg019[1:100, ]
#'     actg036 = actg036[1:50, ]
#'
#'     library(parallel)
#'     ncores    = 2
#'     data.list = list(data = actg019, histdata = actg036)
#'     formula   = cd4 ~ treatment + age + race
#'     family    = poisson()
#'     a0        = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## call created function
#'       ## wrapper to obtain log normalizing constant in parallel package
#'       logncfun = function(a0, ...){
#'         hdbayes::glm.npp.lognc(
#'           formula = formula, family = family, a0 = a0, histdata = data.list[[2]],
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'family', 'data.list'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from normalized power prior
#'       d.npp = glm.npp(
#'         formula = formula,
#'         family = family,
#'         data.list = data.list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = matrix(a0.lognc$lognc, ncol = 1),
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'       glm.logml.npp(
#'         post.samples = d.npp,
#'         bridge.args = list(silent = TRUE)
#'       )
#'     }
#'   }
#' }
glm.logml.npp = function(
    post.samples,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X
  K         = stan.data$K
  oldnames  = paste0("beta[", 1:p, "]")
  newnames  = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('logit_a0s[', 1:(K-1), "]"))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  stan.data$lognc_disp  = pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T)

  ## compute log normalizing constants (lognc) for logit(a0s)
  a0_shape1       = stan.data$a0_shape1
  a0_shape2       = stan.data$a0_shape2
  lognc_logit_a0s = sum( sapply(1:(K-1), function(k){
    lower = stan.data$a0_lower[k]
    upper = stan.data$a0_upper[k]
    if( lower != 0 || upper != 1 ) {
      log( pbeta(upper, shape1 = a0_shape1, shape2 = a0_shape2) - pbeta(lower, shape1 = a0_shape1, shape2 = a0_shape2) )
    }else{
      0
    }
  }) )
  stan.data$lognc_logit_a0s = lognc_logit_a0s

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    a0_shape1  = data$a0_shape1
    a0_shape2  = data$a0_shape2
    a0_lower   = data$a0_lower
    a0_upper   = data$a0_upper
    y          = data$y
    X          = data$X
    dist       = data$dist
    link       = data$link
    offs       = data$offs

    beta       = pars[paste0("beta[", 1:p,"]")]
    logit_a0s  = pars[paste0('logit_a0s[', 1:(K-1), "]")]
    a0s        = c(1, binomial('logit')$linkinv(logit_a0s))

    ## initial prior on beta
    prior_lp   = sum( dnorm(beta, mean = data$mean_beta, sd = data$sd_beta, log = T) )

    ## prior on logit(a0)
    prior_lp   = prior_lp +
      sum( sapply(logit_a0s, logit_beta_lp, shape1 = a0_shape1, shape2 = a0_shape2) ) - data$lognc_logit_a0s

    dispersion = 1
    if ( dist > 2 ){
      ## prior on dispersion
      dispersion = pars[["dispersion"]]
      prior_lp   = prior_lp +
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) - data$lognc_disp
    }

    data_lp = as.numeric( a0s %*% sapply(1:K, function(k){
      start.idx = data$start_idx[k]
      end.idx   = data$end_idx[k]
      y         = data$y[ start.idx:end.idx ]
      X         = data$X[ start.idx:end.idx, ]
      offs      = data$offs[ start.idx:end.idx ]
      glm_lp(y, beta, X, dist, link, offs, dispersion)
    }) )

    ## subtract log nc from power prior
    a0_lognc = data$a0_lognc
    lognc    = data$lognc
    prior_lp = prior_lp - sum( sapply(1:(K-1), function(k){
      pp_lognc(a0s[k+1], a0_lognc, lognc[, k])
    }) )
    return(data_lp + prior_lp)
  }

  lb           = rep(-Inf, p)
  ub           = rep(Inf, p)
  if( stan.data$dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
  }
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
