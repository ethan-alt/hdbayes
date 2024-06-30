#' Log marginal likelihood of a GLM under non-informative reference prior
#'
#' Uses bridge sampling to estimate the logarithm of the marginal likelihood of a GLM under the non-informative reference
#' prior (also referred to as the vague prior).
#'
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      output from [glm.reference()] giving posterior samples of a GLM under the non-informative
#'                          reference prior, with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"Reference"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` containing the output from using [bridgesampling::bridge_sampler()]
#'    to compute the logarithm of the marginal likelihood of the non-informative reference prior}
#'  }
#'
#' @references
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   actg019 = actg019[1:100, ]
#'   data.list = list(currdata = actg019)
#'   formula = cd4 ~ treatment + age + race
#'   family = poisson('log')
#'   d.ref = glm.reference(
#'     formula = formula, family = family,
#'     data.list = data.list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.reference(
#'     post.samples = d.ref,
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.logml.reference = function(
    post.samples,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X
  oldnames  = paste0("beta[", 1:p, "]")
  newnames  = colnames(X)
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
    data_lp = glm_lp(data$y, beta, data$X, dist, link, data$offs, dispersion)
    return(data_lp + prior_lp)
  }

  lb           = rep(-Inf, p)
  ub           = rep(Inf, p)
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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "Reference",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
