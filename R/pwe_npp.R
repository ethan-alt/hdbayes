#' Posterior of normalized power prior (NPP)
#'
#' Sample from the posterior distribution of a piecewise exponential (PWE) model using the normalized
#' power prior (NPP) by Duan et al. (2006) <doi:10.1002/env.752>.
#'
#' Before using this function, users must estimate the logarithm of the normalizing constant across a
#' range of different values for the power prior parameter (\eqn{a_0}), possibly smoothing techniques
#' over a find grid. The power prior parameters (\eqn{a_0}'s) are treated as random with independent
#' beta priors. The initial priors on the regression coefficients are independent normal priors. The
#' current and historical data models are assumed to share the baseline hazard parameters with
#' half-normal priors.
#'
#' @include data_checks_pwe.R
#' @include get_stan_data_pwe.R
#' @include pwe_npp_lognc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#'                          The response is a survival object as returned by the `survival::Surv(time, event)` function,
#'                          where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type of
#'                          censoring is assumed to be right-censoring.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For fitting piecewise exponential (PWE) models, all historical
#'                          data sets will be stacked into one historical data set.
#' @param a0.lognc          a vector giving values of the power prior parameter for which the logarithm of the normalizing
#'                          constant has been evaluated.
#' @param lognc             a vector giving the logarithm of the normalizing constant (as estimated by [pwe.npp.lognc()] for
#'                          each value of `a0.lognc` using the historical data set.
#' @param breaks            a numeric vector specifying the time points that define the boundaries of the piecewise
#'                          intervals. The values should be in ascending order, with the final value being greater than
#'                          or equal to the maximum observed time.
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
#' @param a0.shape1         first shape parameter for the beta prior on the power prior parameter (\eqn{a_0}). When
#'                          \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the beta prior on the power prior parameter (\eqn{a_0}). When
#'                          \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.lower          a scalar giving the lower bound for \eqn{a_0}. Defaults to 0.
#' @param a0.upper          a scalar giving the upper bound for \eqn{a_0}. Defaults to 1.
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
#' @seealso [pwe.npp.lognc()]
#'
#' @references
#'  Duan, Y., Ye, K., and Smith, E. P. (2005). Evaluating water quality using power priors to incorporate historical information. Environmetrics, 17(1), 95–106.
#'
#' @examples
#' \donttest{
#'   if(requireNamespace("parallel")){
#'     library(parallel)
#'     ncores    = 2
#'
#'     if(requireNamespace("survival")){
#'       library(survival)
#'       data(E1684)
#'       data(E1690)
#'       ## take subset for speed purposes
#'       E1684 = E1684[1:100, ]
#'       E1690 = E1690[1:50, ]
#'       ## replace 0 failure times with 0.50 days
#'       E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'       E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'       E1684$cage = as.numeric(scale(E1684$age))
#'       E1690$cage = as.numeric(scale(E1690$age))
#'       data_list = list(currdata = E1690, histdata = E1684)
#'       nbreaks = 3
#'       probs   = 1:nbreaks / nbreaks
#'       breaks  = as.numeric(
#'         quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'       )
#'       breaks  = c(0, breaks)
#'       breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin
#'     }
#'
#'     a0 = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## call created function
#'       ## wrapper to obtain log normalizing constant in parallel package
#'       logncfun = function(a0, ...){
#'         hdbayes::pwe.npp.lognc(
#'           formula = formula, histdata = data_list[[2]], breaks = breaks, a0 = a0,
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'data_list', 'breaks'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from normalized power prior
#'       pwe.npp(
#'         formula = formula,
#'         data.list = data_list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = a0.lognc$lognc,
#'         breaks = breaks,
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'     }
#'   }
#' }
pwe.npp = function(
    formula,
    data.list,
    a0.lognc,
    lognc,
    breaks,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = 0,
    a0.upper          = 1,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NPP
  standat = get.pwe.stan.data.npp(
    formula          = formula,
    data.list        = data.list,
    a0.lognc         = a0.lognc,
    lognc            = lognc,
    breaks           = breaks,
    beta.mean        = beta.mean,
    beta.sd          = beta.sd,
    base.hazard.mean = base.hazard.mean,
    base.hazard.sd   = base.hazard.sd,
    a0.shape1        = a0.shape1,
    a0.shape2        = a0.shape2,
    a0.lower         = a0.lower,
    a0.upper         = a0.upper,
    get.loglik       = get.loglik
  )

  pwe_npp = instantiate::stan_package_model(
    name = "pwe_npp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = pwe_npp$sample(data = standat,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)

  ## rename parameters
  p        = standat$p
  X1       = standat$X1
  J        = standat$J
  if( p > 0 ){
    oldnames = c(paste0("beta[", 1:p, "]"), paste0("lambda[", 1:J, "]"), "a0")
    newnames = c(colnames(X1), paste0("basehaz[", 1:J, "]"), "a0")
  }else{
    oldnames = c(paste0("lambda[", 1:J, "]"), "a0")
    newnames = c(paste0("basehaz[", 1:J, "]"), "a0")
  }

  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
