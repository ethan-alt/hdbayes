#' Log marginal likelihood of a GLM under power prior (PP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under PP. The function's arguments mirror those of the function [glm.pp()],
#' including two additional arguments: `post.samples` and `bridge.args`. These provide posterior samples
#' from the GLM under PP (for example, the output from [glm.pp()]) and specify arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, `ub`), respectively.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include glm_npp_lognc.R
#'
#' @export
#'
#' @inheritParams glm.pp
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM using the PP, such as the output from [glm.pp()]. Each row corresponds to the
#'                          posterior samples obtained from one iteration of MCMC. The column names of `post.samples` should
#'                          include the names of covariates for regression coefficients, such as "(Intercept)", and
#'                          "dispersion" for the dispersion parameter, if applicable.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"PP" if at least one of the power prior parameters (\eqn{a_0}'s) is non-zero, and "PP (No Borrow)" otherwise.}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
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
#'   formula = cd4 ~ treatment + age + race
#'   family = poisson('log')
#'   a0 = 0.5
#'   d.pp = glm.pp(
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     a0.vals = a0,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.pp(
#'     post.samples = d.pp,
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     a0.vals = a0,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.logml.pp = function(
    post.samples,
    formula,
    family,
    data.list,
    a0.vals,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for PP
  standat = get.stan.data.pp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    a0.vals     = a0.vals,
    offset.list = offset.list,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd
  )

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family)

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = standat$p
  X        = standat$X
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion')
  }
  d = d[, oldnames, drop=F]

  ## estimate log normalizing constant for PP
  lognc   = 0
  K       = standat$K
  a0_vals = standat$a0_vals
  if ( K >= 2 & sum(a0_vals) > 1 ){
    lognc = sapply(2:K, function(k){
      start.index = standat$start_idx[k]
      end.index   = standat$end_idx[k]
      offs0       = standat$offs[start.index:end.index]

      glm.npp.lognc(
        formula = formula,
        family = family,
        histdata = data.list[[k]],
        a0 = a0_vals[k],
        offset0 = offs0,
        beta.mean = standat$mean_beta,
        beta.sd = standat$sd_beta,
        disp.mean = standat$disp_mean,
        disp.sd = standat$disp_sd,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        chains = chains,
        ...
      )["lognc"]
    })
    lognc = sum(lognc)
  }
  standat$lognc = lognc

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
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) -
        pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T)
    }
    data_lp = 0
    for ( k in 1:data$K ){
      start.idx = data$start_idx[k]
      end.idx   = data$end_idx[k]
      y         = data$y[ start.idx:end.idx ]
      X         = data$X[ start.idx:end.idx, ]
      offs      = data$offs[ start.idx:end.idx ]
      data_lp   = data_lp + data$a0_vals[k] * glm_lp(y, beta, X, dist, link, offs, dispersion)
    }
    return(data_lp + prior_lp - data$lognc)
  }

  lb           = rep(-Inf, p)
  ub           = rep(Inf, p)
  if( standat$dist > 2 ) {
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
        'data' = standat,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## Return a list of model name and log marginal likelihood
  res = list(
    'model' = "PP",
    'logml' = bs$logml
  )
  if( sum(a0_vals) == 1 )
    res[["model"]] = "PP (No borrow)"
  return(res)
}



#' Log marginal likelihood of a GLM under Bayesian hierarchical model (BHM)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under BHM. The function's arguments mirror those of the function [glm.bhm()],
#' including two additional arguments: `post.samples` and `bridge.args`. These provide posterior samples
#' from the GLM under BHM (for example, the output from [glm.bhm()]) and specify arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, `ub`), respectively.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @inheritParams glm.bhm
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM using the PP, such as the output from [glm.pp()]. Each row corresponds to the
#'                          posterior samples obtained from one iteration of MCMC. The column names of `post.samples` should
#'                          include the names of covariates for regression coefficients, such as "(Intercept)", and
#'                          "dispersion" for the dispersion parameter, if applicable.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"BHM"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
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
#'   glm.logml.bhm(
#'     post.samples = d.bhm,
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.logml.bhm = function(
    post.samples,
    formula,
    family,
    data.list,
    offset.list       = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for BHM
  standat = get.stan.data.bhm(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    offset.list    = offset.list,
    meta.mean.mean = meta.mean.mean,
    meta.mean.sd   = meta.mean.sd,
    meta.sd.mean   = meta.sd.mean,
    meta.sd.sd     = meta.sd.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd
  )

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family)

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = standat$p
  K        = standat$K
  X        = standat$X
  oldnames = paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( K == 1 ) {
    newnames = colnames(X)
  }else {
    newnames = c(colnames(X), paste0( colnames(X), '_hist_', rep(1:(K-1), each = p) ))
  }
  colnames(d)[colnames(d) %in% newnames] = oldnames
  oldnames = c(oldnames, paste0("beta_mean[", 1:p, "]"), paste0("beta_sd[", 1:p, "]"))
  if ( !family$family %in% c('binomial', 'poisson') ) {
    if (K == 1) {
      oldnames = c(oldnames, 'dispersion')
    }else {
      oldnames = c(oldnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
    }
  }
  d = d[, oldnames, drop=F]

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta       = pars[paste0("beta[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta       = matrix(beta, nrow = p, ncol = K)
    beta_mean  = pars[paste0("beta_mean[", 1:p, "]")]
    beta_sd    = pars[paste0("beta_sd[", 1:p, "]")]
    ## prior on beta_mean
    prior_lp   = sum( dnorm(beta_mean, mean = data$meta_mean_mean, sd = data$meta_mean_sd, log = T) )
    ## prior on beta_sd
    prior_lp   = prior_lp + sum(
      dnorm(beta_sd, mean = data$meta_sd_mean, sd = data$meta_sd_sd, log = T) -
        pnorm(0, mean = data$meta_sd_mean, sd = data$meta_sd_sd, lower.tail = F, log.p = T)
    )
    ## prior on beta
    for ( j in 1:p ){
      prior_lp = prior_lp + sum( dnorm(beta[j, ], mean = beta_mean[j], sd = beta_sd[j], log = T) )
    }
    dist       = data$dist
    link       = data$link
    dispersion = rep(1.0, K)
    if ( dist > 2 ){
      ## prior on dispersion
      dispersion = pars[c("dispersion", paste0( "dispersion", "_hist_", 1:(K-1) ))]
      prior_lp   = prior_lp + sum(
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) -
          pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T)
      )
    }
    data_lp = 0
    for ( k in 1:K ){
      start.idx = data$start_idx[k]
      end.idx   = data$end_idx[k]
      y         = data$y[ start.idx:end.idx ]
      X         = data$X[ start.idx:end.idx, ]
      offs      = data$offs[ start.idx:end.idx ]
      data_lp   = data_lp + glm_lp(y, beta[, k], X, dist, link, offs, dispersion[k])
    }
    return(data_lp + prior_lp)
  }

  lb           = rep(-Inf, p*(K+2))
  ub           = rep(Inf, p*(K+2))
  if( standat$dist > 2 ) {
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
        'data' = standat,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## Return a list of model name and log marginal likelihood
  res = list(
    'model' = "BHM",
    'logml' = bs$logml
  )
  return(res)
}
