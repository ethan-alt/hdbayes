#' Log marginal likelihood of a GLM under Bayesian hierarchical model (BHM))
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under BHM.
#'
#' This function shares the same arguments as [glm.bhm()], while introducing two additional parameters:
#' `post.samples` and `bridge.args`. `post.samples`, if not NULL, provides posterior samples under BHM
#' (e.g., the output from [glm.bhm()]). In case `post.samples` is set to NULL, this function will draw
#' posterior samples under BHM using [glm.bhm()]. `bridge.args` specifies arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`).
#'
#' If `post.samples` is not NULL, it is important to ensure that the values assigned to the shared
#' arguments (excluding those relevant for MCMC sampling) in this function and [glm.bhm()] match with
#' those used in generating `post.samples`. When `post.samples` is set to NULL, the arguments pertinent
#' to MCMC sampling are utilized to draw posterior samples of the GLM under BHM.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include glm_bhm.R
#'
#' @export
#'
#' @inheritParams glm.bhm
#' @param post.samples      If not NULL, an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame`
#'                          giving posterior samples of a GLM under BHM, such as the output from [glm.bhm()].
#'                          Each row of `post.samples` corresponds to the posterior samples obtained from one
#'                          iteration of MCMC. The column names of `post.samples` should include the names of
#'                          covariates for regression coefficients, such as "(Intercept)", and "dispersion" for
#'                          the dispersion parameter, if applicable. Defaults to NULL.
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
#'   glm.logml.bhm(
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     post.samples = d.bhm,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.bhm = function(
    formula,
    family,
    data.list,
    post.samples      = NULL,
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

  if ( is.null(post.samples) ){
    ## fit BHM if post.samples is NULL
    post.samples = glm.bhm(
      formula           = formula,
      family            = family,
      data.list         = data.list,
      offset.list       = offset.list,
      meta.mean.mean    = meta.mean.mean,
      meta.mean.sd      = meta.mean.sd,
      meta.sd.mean      = meta.sd.mean,
      meta.sd.sd        = meta.sd.sd,
      disp.mean         = disp.mean,
      disp.sd           = disp.sd,
      iter_warmup       = iter_warmup,
      iter_sampling     = iter_sampling,
      chains            = chains,
      ...
    )
  }

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family)

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = standat$p
  K        = standat$K ## number of historical data sets
  oldnames = paste0("beta_raw[", rep(1:p, K), ',', rep(1:K, each = p), "]")
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
    beta_raw   = pars[paste0("beta_raw[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta_mean  = pars[paste0("beta_mean[", 1:p, "]")]
    beta_sd    = pars[paste0("beta_sd[", 1:p, "]")]
    ## prior on beta_mean and beta_sd
    prior_lp   = sum( dnorm(beta_mean, mean = data$meta_mean_mean, sd = data$meta_mean_sd, log = T) +
      dnorm(beta_sd, mean = data$meta_sd_mean, sd = data$meta_sd_sd, log = T) -
        pnorm(0, mean = data$meta_sd_mean, sd = data$meta_sd_sd, lower.tail = F, log.p = T)
    )
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
      prior_lp   = prior_lp + sum(
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) -
          pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T)
      )
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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "BHM",
    'logml' = bs$logml,
    'bs'    = bs
  )
}
