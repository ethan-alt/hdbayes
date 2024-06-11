#' Estimate the logarithm of the normalizing constant for Latent Exchangeability Prior (LEAP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the normalizing
#' constant for LEAP using historical data sets.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param hist.data.list    a list of `data.frame`s giving historical data sets.
#' @param K                 the desired number of classes to identify. Defaults to 2.
#' @param prob.conc         a scalar or a vector of length `K` giving the concentration parameters for Dirichlet prior.
#'                          If length == 2, a `beta(prob.conc[1], prob.conc[2])` prior is used. If a scalar is provided,
#'                          prob.conc will be a vector of repeated elements of the given scalar. Defaults to a vector of 1s.
#' @param hist.offset.list  a list of matrices giving the offset for historical data. For each matrix, the number of
#'                          rows corresponds to observations and columns correspond to classes. Defaults to a list of
#'                          matrices of 0s.
#' @param beta.mean         a scalar or a `p x K` matrix of mean parameters for initial prior on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept). If a scalar is provided,
#'                          beta.mean will be a matrix of repeated elements of the given scalar. Defaults to a matrix of 0s.
#' @param beta.sd           a scalar or a `p x K` matrix of sd parameters for the initial prior on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept). If a scalar is provided,
#'                          same as for beta.mean. Defaults to a matrix of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the means for
#'                          the half-normal hyperpriors on the dispersion parameters. If a scalar is provided, disp.mean will
#'                          be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the sds for
#'                          the half-normal hyperpriors on the dispersion parameters. If a scalar is provided, same as for
#'                          disp.mean. Defaults to a vector of 10s.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{lognc}{the estimated logarithm of the normalizing constant}
#'
#'    \item{min_ess_bulk}{the minimum estimated bulk effective sample size of the MCMC sampling}
#'
#'    \item{max_Rhat}{the maximum Rhat}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from [bridgesampling::bridge_sampler()]}
#'  }
#'
#' @references
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg036 = actg036[1:50, ]
#'   glm.leap.lognc(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     hist.data.list = list(histdata = actg036),
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 2000
#'   )
#' }
glm.leap.lognc = function(
    formula,
    family,
    hist.data.list,
    K                 = 2,
    prob.conc         = NULL,
    hist.offset.list  = NULL,
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
  ## get Stan data for LEAP
  standat = get.stan.data.leap(
    formula        = formula,
    family         = family,
    data.list      = hist.data.list,
    K              = K,
    prob.conc      = prob.conc,
    offset.list    = hist.offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    all.hist       = TRUE
  )

  glm_leap_prior = instantiate::stan_package_model(
    name = "glm_leap_prior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit      = glm_leap_prior$sample(data = standat,
                                   iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                   ...)
  d        = fit$draws(format = 'draws_df')
  summ     = posterior::summarise_draws(d)

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
  d = suppressWarnings(
    as.matrix(d[, oldnames, drop=F])
  )

  ## compute log normalizing constants (lognc) for half-normal prior on dispersion
  standat$lognc_disp  = sum( pnorm(0, mean = standat$disp_mean, sd = standat$disp_sd, lower.tail = F, log.p = T) )

  ## estimate log normalizing constant
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
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
    }else {
      dispersion = rep(1.0, K)
    }
    hist_lp      = glm_mixture_lp(data$y0, betaMat, dispersion, probs, data$X0, dist, link, data$offs0)
    return(hist_lp + prior_lp)
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
        "samples"       = d,
        'log_posterior' = log_density,
        'data'          = standat,
        'lb'            = lb,
        'ub'            = ub
        ),
      bridge.args
    )
  )

  ## Return a list of lognc, min_n_eff, max_Rhat, and output from bridgesampling::bridge_sampler
  res = list(
    'lognc'        = bs$logml,
    'min_ess_bulk' = min(summ[, 'ess_bulk'], na.rm = T),
    'max_Rhat'     = max(summ[, 'rhat'], na.rm = T),
    'bs'           = bs
  )

  if ( res[['min_ess_bulk']] < 1000 )
    warning(
      paste0(
        'The minimum bulk effective sample size of the MCMC sampling is ',
        round(res[['min_ess_bulk']], 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res[['max_Rhat']] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        round(res[['max_Rhat']], 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}

