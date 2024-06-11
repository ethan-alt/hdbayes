#' Log marginal likelihood of a GLM under normalized power prior (NPP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under NPP.
#'
#' This function mirrors the arguments of [glm.npp()] (except for those relevant for MCMC sampling),
#' and introduces two additional arguments: `post.samples` and `bridge.args`. `post.samples` provides
#' posterior samples from the GLM under NPP (e.g., the output from [glm.npp()]), while `bridge.args`
#' specifies arguments to pass onto [bridgesampling::bridge_sampler()] (other than `samples`,
#' `log_posterior`, `data`, `lb`, `ub`).
#'
#' It is crucial to ensure that the values assigned to the shared arguments in this function and
#' [glm.npp()] correspond to those used in generating `post.samples`.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under NPP, such as the output from [glm.npp()]. Each row corresponds to
#'                          the posterior samples obtained from one iteration of MCMC. The column names of `post.samples`
#'                          should include the names of covariates for regression coefficients, such as "(Intercept)", and
#'                          "dispersion" for the dispersion parameter, if applicable.
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for beta.mean. Defaults to a vector of 10s.
#' @param disp.mean         location parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion parameter. Defaults to 10.
#' @param a0.lognc          a vector giving values of the power prior parameter for which the logarithm of the normalizing
#'                          constant has been evaluated.
#' @param lognc             an S by T matrix where S is the length of a0.lognc, T is the number of historical data sets, and
#'                          the j-th column, j = 1, ..., T, is a vector giving the logarithm of the normalizing constant (as
#'                          estimated by [glm.npp.lognc()] for a0.lognc using the j-th historical data set.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.lower          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          lower bounds for each element of the a0 vector. If a scalar is provided, a0.lower will be a
#'                          vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param a0.upper          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          upper bounds for each element of the a0 vector. If a scalar is provided, same as for a0.lower.
#'                          Defaults to a vector of 1s.
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
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from [bridgesampling::bridge_sampler()]}
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
#'         formula = formula, family = family,
#'         data.list = data.list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = matrix(a0.lognc$lognc, ncol = 1),
#'         bridge.args = list(silent = TRUE)
#'       )
#'     }
#'   }
#' }
glm.logml.npp = function(
    post.samples,
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    bridge.args       = NULL
) {
  ## get Stan data for NPP
  standat = get.stan.data.npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper
  )

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family)

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = standat$p
  X        = standat$X
  K        = standat$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('logit_a0s[', 1:(K-1), "]"))
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  standat$lognc_disp  = pnorm(0, mean = standat$disp_mean, sd = standat$disp_sd, lower.tail = F, log.p = T)

  ## compute log normalizing constants (lognc) for logit(a0s)
  a0_shape1       = standat$a0_shape1
  a0_shape2       = standat$a0_shape2
  lognc_logit_a0s = sum( sapply(1:(K-1), function(k){
    lower = standat$a0_lower[k]
    upper = standat$a0_upper[k]
    if( lower != 0 || upper != 1 ) {
      log( pbeta(upper, shape1 = a0_shape1, shape2 = a0_shape2) - pbeta(lower, shape1 = a0_shape1, shape2 = a0_shape2) )
    }else{
      0
    }
  }) )
  standat$lognc_logit_a0s = lognc_logit_a0s

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
  if( standat$dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
  }
  lb        = c(lb, binomial('logit')$linkfun(standat$a0_lower))
  ub        = c(ub, binomial('logit')$linkfun(standat$a0_upper))
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
    'model' = "NPP",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
