#' Estimate the logarithm of the normalizing constant for normalized power prior (NPP) for one data set
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the normalizing
#' constant of a mixture cure rate (CurePWE) model under the NPP for a fixed value of the power prior
#' parameter \eqn{a_0 \in (0, 1)} for one data set. The initial priors are independent normal priors on the
#' regression coefficients and half-normal priors on the baseline hazard parameters. Additionally, a normal
#' prior is specified for the logit of the cure fraction \eqn{\pi}.
#'
#' @include data_checks_pwe.R
#' @include curepwe_pp_lognc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates in
#'                          the PWE model. The response is a survival object as returned by the `survival::Surv(time, event)`
#'                          function, where event is a binary indicator for event (0 = no event, 1 = event has occurred).
#'                          The type of censoring is assumed to be right-censoring.
#' @param histdata          a `data.frame` giving the historical data.
#' @param breaks            a numeric vector specifying the time points that define the boundaries of the piecewise
#'                          intervals. The values should be in ascending order, with the final value being greater than
#'                          or equal to the maximum observed time.
#' @param a0                a scalar between 0 and 1 giving the (fixed) power prior parameter for the historical data.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          `beta.mean` will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for `beta.mean`. Defaults to a vector of 10s.
#' @param base.hazard.mean  a scalar or a vector whose dimension is equal to the number of intervals giving the location
#'                          parameters for the half-normal priors on the baseline hazards of the PWE model. If a scalar is
#'                          provided, same as for `beta.mean`. Defaults to 0.
#' @param base.hazard.sd    a scalar or a vector whose dimension is equal to the number of intervals giving the scale
#'                          parameters for the half-normal priors on the baseline hazards of the PWE model. If a scalar is
#'                          provided, same as for `beta.mean`. Defaults to 10.
#' @param logit.pcured.mean mean parameter for the normal prior on the logit of the cure fraction \eqn{\pi}. Defaults to 0.
#' @param logit.pcured.sd   sd parameter for the normal prior on the logit of the cure fraction \eqn{\pi}. Defaults to 3.
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1684[E1684$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     curepwe.npp.lognc(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       histdata = E1684,
#'       breaks = breaks,
#'       a0 = 0.5,
#'       logit.pcured.mean = 0, logit.pcured.sd = 3,
#'       bridge.args = list(silent = TRUE),
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
curepwe.npp.lognc = function(
    formula,
    histdata,
    breaks,
    a0,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    logit.pcured.mean = NULL,
    logit.pcured.sd   = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks.pwe(formula, list(histdata), breaks)

  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  y0            = histdata[, time.name]
  eventind0     = as.integer( histdata[, eventind.name] )
  X0            = stats::model.matrix(formula, histdata)
  if ( '(Intercept)' %in% colnames(X0) )
    X0 = X0[, -1, drop = F]

  p = ncol(X0)
  J = length(breaks) - 1  ## number of intervals

  ## create index giving interval into which obs failed / was censored
  intindx0 = sapply(y0, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })

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

  ## Default half-normal priors on baseline hazards are N^{+}(0, 10^2)
  if ( !is.null(base.hazard.mean) ){
    if ( !( is.vector(base.hazard.mean) & (length(base.hazard.mean) %in% c(1, J)) ) )
      stop("base.hazard.mean must be a scalar or a vector of length ", J, " if base.hazard.mean is not NULL")
  }
  base.hazard.mean = to.vector(param = base.hazard.mean, default.value = 0, len = J)
  if ( !is.null(base.hazard.sd) ){
    if ( !( is.vector(base.hazard.sd) & (length(base.hazard.sd) %in% c(1, J)) ) )
      stop("base.hazard.sd must be a scalar or a vector of length ", J, " if base.hazard.sd is not NULL")
  }
  base.hazard.sd = to.vector(param = base.hazard.sd, default.value = 10, len = J)

  ## Default prior on logit(cure fraction) is N(0, 3^2)
  if ( !is.null(logit.pcured.mean) ){
    if ( !( is.vector(logit.pcured.mean) & (length(logit.pcured.mean) == 1) ) )
      stop("logit.pcured.mean must be a scalar if logit.pcured.mean is not NULL")
  }
  logit.pcured.mean = to.vector(param = logit.pcured.mean, default.value = 0, len = 1)
  if ( !is.null(logit.pcured.sd) ){
    if ( !( is.vector(logit.pcured.sd) & (length(logit.pcured.sd) == 1) ) )
      stop("logit.pcured.sd must be a scalar if logit.pcured.sd is not NULL")
  }
  logit.pcured.sd = to.vector(param = logit.pcured.sd, default.value = 3, len = 1)

  ## check a0 value
  a0 = as.numeric(a0)
  if ( any(a0 < 0 | a0 > 1 ) )
    stop("a0 must be a scalar between 0 and 1")

  standat = list(
    'n0'                 = nrow(X0),
    'J'                  = J,
    'p'                  = p,
    'y0'                 = y0,
    'X0'                 = X0,
    'intindx0'           = intindx0,
    'death_ind0'         = eventind0,
    'breaks'             = breaks,
    'a0'                 = a0,
    'beta_mean'          = beta.mean,
    'beta_sd'            = beta.sd,
    'hazard_mean'        = base.hazard.mean,
    'hazard_sd'          = base.hazard.sd,
    'logit_p_cured_mean' = logit.pcured.mean,
    'logit_p_cured_sd'   = logit.pcured.sd
  )

  curepwe_pp_prior = instantiate::stan_package_model(
    name = "curepwe_pp_prior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = curepwe_pp_prior$sample(data = standat,
                                iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                ...)
  d    = fit$draws(format = 'draws_df')
  attr(x = d, which = 'data') = standat
  summ = posterior::summarise_draws(d)

  ## compute log normalizing constant
  res.hist = curepwe.pp.lognc(
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
