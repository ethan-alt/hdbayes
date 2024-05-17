#' Log marginal likelihood of a GLM under meta-analytic predictive (MAP) prior (the prior induced by
#' Bayesian hierarchical model (BHM))
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under MAP.
#'
#' This function shares the same arguments as [glm.bhm()], while introducing two additional parameters:
#' `post.samples` and `bridge.args`. `post.samples` provides posterior samples of the GLM under BHM
#' (e.g., the output from [glm.bhm()]), whereas `bridge.args` specifies arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`).
#'
#' It is important to ensure that the values assigned to the shared arguments (excluding those relevant
#' for MCMC sampling) in this function and [glm.bhm()] align with those used in generating `post.samples`.
#' The arguments pertinent to MCMC sampling are utilized to compute the log marginal likelihood under
#' BHM using only historical data sets.
#'
#' @include glm_logml_bhm.R
#'
#' @export
#'
#' @inheritParams glm.bhm
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under BHM, such as the output from [glm.bhm()]. Each row corresponds to the
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
#'    \item{model}{"MAP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using all data sets (including current and historical data)}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using historical data sets}
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
#'   glm.logml.map(
#'     post.samples = d.bhm,
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.map = function(
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
  K = length(data.list)
  if ( K == 1 ){
    stop("data.list should include at least one historical data set")
  }
  ## computing normalizing constant using all data sets
  res.all = glm.logml.bhm(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    post.samples   = post.samples,
    offset.list    = offset.list,
    meta.mean.mean = meta.mean.mean,
    meta.mean.sd   = meta.mean.sd,
    meta.sd.mean   = meta.sd.mean,
    meta.sd.sd     = meta.sd.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    bridge.args    = bridge.args
  )

  ## computing normalizing constant using historical data sets only
  if ( !is.null(disp.mean) && (length(disp.mean) != 1) ){
    disp.mean = disp.mean[-1]
  }
  if ( !is.null(disp.sd) && (length(disp.sd) != 1) ){
    disp.sd = disp.sd[-1]
  }
  res.hist = glm.logml.bhm(
    formula        = formula,
    family         = family,
    data.list      = data.list[-1],
    post.samples   = NULL,
    offset.list    = offset.list[-1],
    meta.mean.mean = meta.mean.mean,
    meta.mean.sd   = meta.mean.sd,
    meta.sd.mean   = meta.sd.mean,
    meta.sd.sd     = meta.sd.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    bridge.args    = bridge.args,
    iter_warmup    = iter_warmup,
    iter_sampling  = iter_sampling,
    chains         = chains,
    ...
  )

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model'   = "MAP",
    'logml'   = res.all$logml - res.hist$logml,
    'bs'      = res.all$bs,
    'bs.hist' = res.hist$bs
  )
  return(res)
}
