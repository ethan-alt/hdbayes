#' Log marginal likelihood of a GLM under Latent Exchangeability Prior (LEAP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under LEAP.
#'
#' This function shares the same arguments as [glm.leap()], while introducing two additional
#' parameters: `post.samples` and `bridge.args`. `post.samples` provides posterior samples under LEAP
#' (e.g., the output from [glm.leap()]), whereas `bridge.args` specifies arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`).
#'
#' It is important to ensure that the values assigned to the shared arguments (excluding those relevant
#' for MCMC sampling) in this function and [glm.leap()] match with those used in generating
#' `post.samples`. The arguments pertinent to MCMC sampling are utilized to compute the normalizing
#' constants for LEAP.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#' @include glm_leap_lognc.R
#'
#' @export
#'
#' @inheritParams glm.leap
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under LEAP, such as the output from [glm.leap()]. Each row corresponds
#'                          to the posterior samples obtained from one iteration of MCMC. The column names of `post.samples`
#'                          should include the names of covariates for regression coefficients, such as "(Intercept)",
#'                          and "dispersion" for the dispersion parameter, if applicable.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"LEAP"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using all data sets (including current and historical data)}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using historical data sets (for computing the log normalizing
#'    constant for LEAP)}
#'
#'  }
#'
#' @references
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2023). LEAP: The latent exchangeability prior for borrowing information from historical data. arXiv preprint.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019   = actg019[1:100, ]
#'   actg036   = actg036[1:50, ]
#'   formula   = outcome ~ scale(age) + race + treatment + scale(cd4)
#'   family    = binomial('logit')
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.leap    = glm.leap(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.leap(
#'     post.samples = d.leap,
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     K = 2,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.leap = function(
    post.samples,
    formula,
    family,
    data.list,
    K                 = 2,
    prob.conc         = NULL,
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
  if ( length(data.list) == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for LEAP
  standat = get.stan.data.leap(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    K              = K,
    prob.conc      = prob.conc,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd
  )

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family, is.LEAP = TRUE)

  d        = as.matrix(post.samples)

  p        = standat$p
  K        = standat$K
  oldnames = paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  oldnames = c(oldnames, "gamma")
  if ( K > 2 ){
    oldnames = c(oldnames, paste0("delta_raw[", 1:(K-2), "]"))
  }
  d = d[, oldnames, drop=F]

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    betaMat    = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    prior_lp   = sum( dnorm(betaMat, mean = as.numeric(data$mean_beta),
                            sd = as.numeric(data$sd_beta), log = T) )
    betaMat    = matrix(betaMat, nrow = p, ncol = K)

    ## prior on gamma
    conc         = data$conc
    gamma_shape1 = conc[1]
    gamma_shape2 = sum(conc[2:K])
    gamma        = pars[["gamma"]]
    probs        = c(gamma, 1 - gamma)
    if ( gamma_shape1 != 1 || gamma_shape2 != 1 ){
      prior_lp = prior_lp + dbeta(gamma, gamma_shape1, gamma_shape2, log = T)
    }

    if( K > 2 ){
      delta_raw = pars[paste0("delta_raw[", 1:(K-2), "]")]
      delta_raw = c(delta_raw, 1 - sum(delta_raw))
      prior_lp  = prior_lp + dirichlet_lp(delta_raw, conc[2:K])
      probs     = c(gamma, (1 - gamma) * delta_raw)
    }

    dist         = data$dist
    link         = data$link
    if ( dist > 2 ) {
      dispersion = pars[paste0("dispersion[", 1:K,"]")]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) -
        sum( pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T) )
    }else {
      dispersion = rep(1.0, K)
    }
    ## historical data likelihood
    prior_lp     = prior_lp  + glm_mixture_lp(data$y0, betaMat, dispersion, probs, data$X0, dist, link, data$offs0)
    ## current data likelihood
    data_lp      = glm_lp(data$y, betaMat[, 1], data$X, dist, link, data$offs[,1], dispersion[1])
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p*K)
  ub = rep(Inf, p*K)
  if( standat$dist > 2 ) {
    lb = c(lb, rep(0, K) )
    ub = c(ub, rep(Inf, K) )
  }
  lb = c(lb, standat$gamma_lower)
  ub = c(ub, standat$gamma_upper)
  if( K > 2 ){
    lb = c(lb, rep(0, K-2))
    ub = c(ub, rep(1, K-2))
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

  ## estimate log normalizing constant for LEAP using historical data sets
  res.hist = glm.leap.lognc(
    formula           = formula,
    family            = family,
    hist.data.list    = data.list[-1],
    K                 = K,
    prob.conc         = prob.conc,
    hist.offset.list  = offset.list[-1],
    beta.mean         = beta.mean,
    beta.sd           = beta.sd,
    disp.mean         = disp.mean,
    disp.sd           = disp.sd,
    bridge.args       = bridge.args,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    chains            = chains,
    ...
  )

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model'   = "LEAP",
    'logml'   = bs$logml - res.hist$lognc,
    'bs'      = bs,
    'bs.hist' = res.hist$bs
  )
  return(res)
}
