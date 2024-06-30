#' Log marginal likelihood of a GLM under normalized asymptotic power prior (NAPP)
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of a GLM under the
#' normalized asymptotic power prior (NAPP).
#'
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [glm.napp()] giving posterior samples of a GLM under the normalized asymptotic
#'                          power prior (NAPP), with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"NAPP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the marginal likelihood of the normalized asymptotic power prior (NAPP)}
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
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.logml.napp = function(
    post.samples,
    bridge.args       = NULL
) {
  ## get Stan data for NAPP
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)

  ## rename parameters
  p        = stan.data$p
  X        = stan.data$X
  K        = stan.data$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( stan.data$dist > 2 ) {
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
  if( stan.data$dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
  }
  lb        = c(lb, rep(-Inf, stan.data$K))
  ub        = c(ub, rep(Inf, stan.data$K))
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
    'model' = "NAPP",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
