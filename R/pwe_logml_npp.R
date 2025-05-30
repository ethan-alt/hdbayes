#' Log marginal likelihood of a piecewise exponential (PWE) model under normalized power prior (NPP)
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of a PWE model under the
#' normalized power prior (NPP).
#'
#' @include pwe_loglik.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [pwe.npp()] giving posterior samples of a PWE model under the normalized
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
#'       nbreaks = 3
#'       probs   = 1:nbreaks / nbreaks
#'       breaks  = as.numeric(
#'         quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'       )
#'       breaks  = c(0, breaks)
#'       breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
#'     }
#'
#'     a0 = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## call created function
#'       ## wrapper to obtain log normalizing constant in parallel package
#'       logncfun = function(a0, ...){
#'         hdbayes::pwe.npp.lognc(
#'           formula = formula, histdata = data_list[[2]], breaks = breaks, a0 = a0,
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'data_list', 'breaks'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from normalized power prior
#'       d.npp = pwe.npp(
#'         formula = formula,
#'         data.list = data_list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = a0.lognc$lognc,
#'         breaks = breaks,
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'       pwe.logml.npp(
#'         post.samples = d.npp,
#'         bridge.args = list(silent = TRUE)
#'       )
#'     }
#'   }
#' }
pwe.logml.npp = function(
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
    oldnames = c(paste0("beta[", 1:p, "]"), paste0("lambda[", 1:J, "]"))
    newnames = c(colnames(X1), paste0("basehaz[", 1:J, "]"))
    lb        = c(rep(-Inf, p), rep(0, J), binomial('logit')$linkfun(stan.data$a0_lower))

  }else{
    oldnames = paste0("lambda[", 1:J, "]")
    newnames = paste0("basehaz[", 1:J, "]")
    lb        = c(rep(0, J), binomial('logit')$linkfun(stan.data$a0_lower))
  }
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames  = c(oldnames, "logit_a0")
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on baseline hazards
  stan.data$lognc_hazard = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )

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
    p          = data$p
    lambda     = as.numeric( pars[paste0("lambda[", 1:data$J,"]")] )
    logit_a0   = as.numeric(pars["logit_a0"])
    a0         = binomial('logit')$linkinv(logit_a0)
    ## prior on logit(a0)
    prior_lp   = logit_beta_lp(logit_a0, shape1 = a0_shape1, shape2 = a0_shape2) - data$lognc_logit_a0

    if( p > 0 ){
      beta       = as.numeric( pars[paste0("beta[", 1:p,"]")] )
      prior_lp   = prior_lp + sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) ) +
        sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard
      eta        = data$X1 %*% beta
      eta0       = data$X0 %*% beta
      data_lp    = a0 * sum( pwe_lpdf(data$y0, eta0, lambda, data$breaks, data$intindx0, data$J, data$death_ind0) ) +
        sum( pwe_lpdf(data$y1, eta, lambda, data$breaks, data$intindx, data$J, data$death_ind) )

    }else{
      prior_lp   = prior_lp + sum( dnorm(lambda, mean = data$hazard_mean, sd = data$hazard_sd, log = T) ) - data$lognc_hazard
      data_lp    = a0 * sum( pwe_lpdf(data$y0, 0, lambda, data$breaks, data$intindx0, data$J, data$death_ind0) ) +
        sum( pwe_lpdf(data$y1, 0, lambda, data$breaks, data$intindx, data$J, data$death_ind) )
    }

    ## subtract log nc from power prior
    prior_lp = prior_lp - pp_lognc(a0, data$a0_lognc, data$lognc)
    return(data_lp + prior_lp)
  }

  ub        = c(rep(Inf, length(lb) - 1), binomial('logit')$linkfun(stan.data$a0_upper))
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
