#' Posterior of Bayesian hierarchical model (BHM)
#'
#' Sample from the posterior distribution of a standard cure rate model (referred to as the **CurePWE model**)
#' using the Bayesian hierarchical model (BHM). The CurePWE model assumes that a fraction \eqn{\pi} of the
#' population is "cured", while the remaining \eqn{1 - \pi} are susceptible to the event of interest. The
#' survival function for the entire population is given by:
#' \deqn{S_{\text{pop}}(t) = \pi + (1 - \pi) S(t),}
#' where \eqn{S(t)} represents the survival function of the non-cured individuals. We model \eqn{S(t)} using a
#' piecewise exponential (PWE) model (i.e., a proportional hazards model with a piecewise constant baseline hazard).
#' Covariates are incorporated through the PWE model.
#'
#' The Bayesian hierarchical model (BHM) assumes that the regression coefficients in the PWE models for the
#' historical and current data are different, but are correlated through a common distribution, whose hyperparameters
#' (i.e., mean and standard deviation (sd) (the covariance matrix is assumed to have a diagonal structure))
#' are treated as random. The number of regression coefficients for the current data is assumed to be the
#' same as that for the historical data.
#'
#' The hyperpriors on the mean and the sd hyperparameters are independent normal and independent half-normal
#' distributions, respectively. The baseline hazard parameters for both current and historical data models are
#' assumed to be independent and identically distributed (i.i.d.), each assigned a half-normal prior. Similarly,
#' the cure fractions for both models are treated as i.i.d., with a normal prior specified on the logit of the
#' cure fraction.
#'
#' @include data_checks_pwe.R
#' @include get_stan_data_pwe.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates in
#'                          the PWE model. The response is a survival object as returned by the `survival::Surv(time, event)`
#'                          function, where event is a binary indicator for event (0 = no event, 1 = event has occurred).
#'                          The type of censoring is assumed to be right-censoring.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For fitting CurePWE models, all historical data sets will be
#'                          stacked into one historical data set.
#' @param breaks            a numeric vector specifying the time points that define the boundaries of the piecewise
#'                          intervals. The values should be in ascending order, with the final value being greater than
#'                          or equal to the maximum observed time.
#' @param meta.mean.mean    a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the normal hyperpriors on the mean hyperparameters of regression coefficients.
#'                          If a scalar is provided, `meta.mean.mean` will be a vector of repeated elements of the given
#'                          scalar. Defaults to a vector of 0s.
#' @param meta.mean.sd      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the normal hyperpriors on the mean hyperparameters of regression coefficients. If
#'                          a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 10s.
#' @param meta.sd.mean      a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 0s.
#' @param meta.sd.sd        a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#'                          If a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 1s.
#' @param base.hazard.mean  a scalar or a vector whose dimension is equal to the number of intervals giving the location
#'                          parameters for the half-normal priors on the baseline hazards of the PWE model. If a scalar is
#'                          provided, same as for `meta.mean.mean`. Defaults to 0.
#' @param base.hazard.sd    a scalar or a vector whose dimension is equal to the number of intervals giving the scale
#'                          parameters for the half-normal priors on the baseline hazards of the PWE model. If a scalar is
#'                          provided, same as for `meta.mean.mean`. Defaults to 10.
#' @param logit.pcured.mean mean parameter for the normal prior on the logit of the cure fraction \eqn{\pi}. Defaults to 0.
#' @param logit.pcured.sd   sd parameter for the normal prior on the logit of the cure fraction \eqn{\pi}. Defaults to 3.
#' @param get.loglik        whether to generate log-likelihood matrix. Defaults to FALSE.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`).
#'
#' @return
#'  The function returns an object of class `draws_df` giving posterior samples, with an attribute called 'data' which includes
#'  the list of variables specified in the data block of the Stan program.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1684)
#'     data(E1690)
#'     ## take subset for speed purposes
#'     E1684 = E1684[1:100, ]
#'     E1690 = E1690[1:50, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'     E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'     E1684$cage = as.numeric(scale(E1684$age))
#'     E1690$cage = as.numeric(scale(E1690$age))
#'     data_list = list(currdata = E1690, histdata = E1684)
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     curepwe.bhm(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       logit.pcured.mean = 0, logit.pcured.sd = 3,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
curepwe.bhm = function(
    formula,
    data.list,
    breaks,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    logit.pcured.mean = NULL,
    logit.pcured.sd   = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for BHM
  standat = get.pwe.stan.data.bhm(
    formula          = formula,
    data.list        = data.list,
    breaks           = breaks,
    meta.mean.mean   = meta.mean.mean,
    meta.mean.sd     = meta.mean.sd,
    meta.sd.mean     = meta.sd.mean,
    meta.sd.sd       = meta.sd.sd,
    base.hazard.mean = base.hazard.mean,
    base.hazard.sd   = base.hazard.sd,
    get.loglik       = get.loglik
  )

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
  standat[["logit_p_cured_mean"]] = logit.pcured.mean
  standat[["logit_p_cured_sd"]]   = logit.pcured.sd

  curepwe_bhm = instantiate::stan_package_model(
    name = "curepwe_bhm",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = curepwe_bhm$sample(data = standat,
                           iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                           ...)

  ## rename parameters
  p        = standat$p
  X1       = standat$X1
  J        = standat$J
  oldnames = c("p_cured", "p_cured0",
               paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"),
               paste0("lambda[", 1:J, "]"), paste0("lambda0[", 1:J, "]"))
  newnames = c("p_cured", "p_cured_hist",
               colnames(X1), paste0(colnames(X1), '_hist'),
               paste0("basehaz[", 1:J, "]"), paste0("basehaz_hist[", 1:J, "]"))

  ## reorder parameters so that regression coefficients appear at the top
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
