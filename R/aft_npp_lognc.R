#' Estimate the logarithm of the normalizing constant for normalized power prior (NPP) for one data set
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the normalizing
#' constant for the NPP for a fixed value of the power prior parameter \eqn{a_0 \in (0, 1)} for one data
#' set. The initial priors are independent normal priors on the regression coefficients and a half-normal
#' prior on the scale parameter.
#'
#' @include data_checks_aft.R
#' @include aft_pp_lognc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#'                          The response is a survival object as returned by the `survival::Surv(time, event)` function,
#'                          where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type of
#'                          censoring is assumed to be right-censoring.
#' @param histdata          a `data.frame` giving the historical data.
#' @param a0                a scalar between 0 and 1 giving the (fixed) power prior parameter for the historical data.
#' @param dist              a character indicating the distribution of survival times. Currently, `dist` can be one of the
#'                          following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          `beta.mean` will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for `beta.mean`. Defaults to a vector of 10s.
#' @param scale.mean        location parameter for the half-normal prior on the scale parameter of the AFT model. Defaults to 0.
#' @param scale.sd          scale parameter for the half-normal prior on the scale parameter of the AFT model. Defaults to 10.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`) to
#'                          pass onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup`
#'                          in `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method
#'                          in cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
#' @return
#'  The function returns a vector giving the value of a0, the estimated logarithm of the normalizing constant, the minimum
#'  estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat.
#'
#' @references
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1684)
#'     ## take subset for speed purposes
#'     E1684 = E1684[1:100, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'     E1684$cage = as.numeric(scale(E1684$age))
#'     aft.npp.lognc(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       histdata = E1684,
#'       a0 = 0.5,
#'       dist = "weibull",
#'       bridge.args = list(silent = TRUE),
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
aft.npp.lognc = function(
    formula,
    histdata,
    a0,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks.aft(formula, list(histdata), dist)

  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  t0            = histdata[, time.name]
  y0            = log(t0)
  eventind0     = as.integer( histdata[, eventind.name] )
  X0            = stats::model.matrix(formula, histdata)
  p             = ncol(X0)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = to.vector(param = beta.sd, default.value = 10, len = p)

  ## Default half-normal prior on scale parameter is N^{+}(0, 10^2)
  if ( !is.null(scale.mean) ){
    if ( !( is.vector(scale.mean) & (length(scale.mean) == 1) ) )
      stop("scale.mean must be a scalar if scale.mean is not NULL")
  }
  scale.mean = to.vector(param = scale.mean, default.value = 0, len = 1)
  if ( !is.null(scale.sd) ){
    if ( !( is.vector(scale.sd) & (length(scale.sd) == 1) ) )
      stop("scale.sd must be a scalar if scale.sd is not NULL")
  }
  scale.sd = to.vector(param = scale.sd, default.value = 10, len = 1)

  ## check a0 value
  a0 = as.numeric(a0)
  if ( any(a0 < 0 | a0 > 1 ) )
    stop("a0 must be a scalar between 0 and 1")

  standat = list(
    'dist'            = dist.to.integer(dist),
    'n0_obs'          = sum(eventind0),
    'n0_cen'          = sum(1 - eventind0),
    'p'               = p,
    'y0_obs'          = y0[which(eventind0 == 1)],
    'y0_cen'          = y0[which(eventind0 == 0)],
    'X0_obs'          = X0[which(eventind0 == 1), ],
    'X0_cen'          = X0[which(eventind0 == 0), ],
    'a0'              = a0,
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd
  )

  aft_pp_prior = instantiate::stan_package_model(
    name = "aft_pp_prior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = aft_pp_prior$sample(data = standat,
                            iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                            ...)
  d    = fit$draws(format = 'draws_df')
  attr(x = d, which = 'data') = standat
  summ = posterior::summarise_draws(d)

  ## compute log normalizing constant
  res.hist = aft.pp.lognc(
    post.samples   = d,
    is.prior       = TRUE,
    bridge.args    = bridge.args
  )

  ## Return vector of a0, lognc, min_n_eff, max_Rhat
  res = c(
    'a0'           = a0,
    'lognc'        = res.hist$lognc,
    'min_ess_bulk' = min(summ[, 'ess_bulk']),
    'max_Rhat'     = max(summ[, 'rhat'])
  )

  if ( res['min_ess_bulk'] < 1000 )
    warning(
      paste0(
        'The minimum bulk effective sample size of the MCMC sampling is ',
        res['min_ess_bulk'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res['max_Rhat'] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        res['max_Rhat'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}
