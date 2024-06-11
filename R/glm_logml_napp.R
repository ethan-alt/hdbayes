#' Log marginal likelihood of a GLM under normalized asymptotic power prior (NAPP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under NAPP.
#'
#' This function mirrors the arguments of [glm.napp()] (except for those relevant for MCMC sampling),
#' and introduces two additional arguments: `post.samples` and `bridge.args`. `post.samples` provides
#' posterior samples from the GLM under NAPP (e.g., the output from [glm.napp()]), while `bridge.args`
#' specifies arguments to pass onto [bridgesampling::bridge_sampler()] (other than `samples`,
#' `log_posterior`, `data`, `lb`, `ub`).
#'
#' It is crucial to ensure that the values assigned to the shared arguments in this function and
#' [glm.napp()] correspond to those used in generating `post.samples`.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under NAPP, such as the output from [glm.napp()]. Each row corresponds to
#'                          the posterior samples obtained from one iteration of MCMC. The column names of `post.samples`
#'                          should include the names of covariates for regression coefficients, such as "(Intercept)", and
#'                          "dispersion" for the dispersion parameter, if applicable.
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical datasets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass onto
#'                          [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"NAPP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from [bridgesampling::bridge_sampler()]}
#'  }
#'
#' @references
#'  Ibrahim, J. G., Chen, M., Gwon, Y., and Chen, F. (2015). The power prior: Theory and applications. Statistics in Medicine, 34(28), 3724–3749.
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
#'   d.napp = glm.napp(
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.napp(
#'     post.samples = d.napp,
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.logml.napp = function(
    post.samples,
    formula,
    family,
    data.list,
    offset.list       = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    bridge.args       = NULL
) {
  ## get Stan data for NAPP
  standat = get.stan.data.napp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    offset.list = offset.list,
    a0.shape1   = a0.shape1,
    a0.shape2   = a0.shape2
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
  oldnames = c(oldnames, paste0('logit_a0s[', 1:K, "]"))
  d = d[, oldnames, drop=F]

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    a0_shape1  = data$a0_shape1
    a0_shape2  = data$a0_shape2
    y          = data$y
    X          = data$X
    dist       = data$dist
    link       = data$link
    offs       = data$offs

    beta       = pars[paste0("beta[", 1:p,"]")]
    logit_a0s  = pars[paste0('logit_a0s[', 1:K, "]")]
    a0s        = binomial('logit')$linkinv(logit_a0s)

    ## prior on logit(a0)
    prior_lp   = sum( sapply(logit_a0s, logit_beta_lp, shape1 = a0_shape1, shape2 = a0_shape2) )
    dispersion = 1
    theta      = beta
    if ( dist > 2 ){
      dispersion = pars[["dispersion"]]
      log_disp   = log( dispersion )
      theta      = c(beta, log_disp)
      prior_lp   = prior_lp - log_disp ## jacobian
    }
    data_lp    = glm_lp(y, beta, X, dist, link, offs, dispersion)
    prior_lp   = prior_lp + sum( sapply(1:K, function(k){
      theta_mean  = as.numeric(data$theta_hats[, k])
      theta_covar = as.matrix( 1/a0s[k] * data$theta_covars[k, , ] )
      mvtnorm::dmvnorm(theta, mean = theta_mean, sigma = theta_covar, log = T)
    }) )
    return(data_lp + prior_lp)
  }

  lb           = rep(-Inf, p)
  ub           = rep(Inf, p)
  if( standat$dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
  }
  lb        = c(lb, rep(-Inf, standat$K))
  ub        = c(ub, rep(Inf, standat$K))
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
    'model' = "NAPP",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
