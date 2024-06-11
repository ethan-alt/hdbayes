#' Log marginal likelihood of a GLM under a non-informative reference prior
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under a non-informative reference prior (referred to as Ref prior).
#'
#' This function mirrors the same arguments as [glm.reference()] (except for those relevant for MCMC sampling),
#' and introduces two additional parameters: `post.samples` and `bridge.args`. `post.samples` provides posterior
#' samples generated from the GLM under a Ref prior (e.g., the output from [glm.reference()]), whereas
#' `bridge.args` specifies arguments to pass onto [bridgesampling::bridge_sampler()] (other than `samples`,
#' `log_posterior`, `data`, `lb`, and `ub`).
#'
#' It is important to ensure that the values assigned to the shared arguments in this function and
#' [glm.reference()] align with those used in generating `post.samples`.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under a Ref prior, such as the output from [glm.reference()]. Each row
#'                          corresponds to the posterior samples obtained from one iteration of MCMC. The column names
#'                          of `post.samples` should include the names of covariates for regression coefficients, such
#'                          as "(Intercept)", and "dispersion" for the dispersion parameter, if applicable.
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list consisting of only one `data.frame` containing the current data.
#' @param offset.list       a list consisting of only one vector giving the offset for the current data. The length of
#'                          the vector is equal to the number of rows in the current data. The vector has all values
#'                          set to 0 by default.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the normal prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the normal prior on regression coefficients. If a scalar is provided,
#'                          same as for beta.mean. Defaults to a vector of 10s.
#' @param disp.mean         location parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion parameter. Defaults to 10.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"Reference"}
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
#'     formula = formula, family = family,
#'     data.list = data.list,
#'     bridge.args = list(silent = TRUE)
#'   )
#' }
glm.logml.reference = function(
    post.samples,
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    bridge.args       = NULL
) {
  ## get Stan data for reference prior
  standat = get.stan.data.ref(
    formula     = formula,
    family      = family,
    data.list   = data.list,
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

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  standat$lognc_disp  = pnorm(0, mean = standat$disp_mean, sd = standat$disp_sd, lower.tail = F, log.p = T)

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

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model' = "Reference",
    'logml' = bs$logml,
    'bs'    = bs
  )
  return(res)
}
