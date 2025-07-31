#' Posterior of commensurate prior (CP)
#'
#' Sample from the posterior distribution of an accelerated failure time (AFT) model using the commensurate prior (CP)
#' by Hobbs et al. (2011) <doi:10.1111/j.1541-0420.2011.01564.x>.
#'
#' The commensurate prior (CP) assumes that the regression coefficients for the current data model conditional on those
#' for the historical data model are independent normal distributions with mean equal to the corresponding regression
#' coefficients for the historical data and variance equal to the inverse of the corresponding elements of a vector of
#' precision parameters (referred to as the commensurability parameter \eqn{\tau}). We regard \eqn{\tau} as random and elicit
#' a spike-and-slab prior, which is specified as a mixture of two half-normal priors, on \eqn{\tau}.
#'
#' The number of current data regression coefficients is assumed to be the same as that of historical data regression
#' coefficients. The scale parameters for both current and historical data models are assumed to be independent and
#' identically distributed, each assigned a half-normal prior.
#'
#' @include data_checks_aft.R
#' @include get_stan_data_aft.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#'                          The response is a survival object as returned by the `survival::Surv(time, event)` function,
#'                          where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type of
#'                          censoring is assumed to be right-censoring.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For fitting accelerated failure time (AFT) models, all historical
#'                          data sets will be stacked into one historical data set.
#' @param dist              a character indicating the distribution of survival times. Currently, `dist` can be one of the
#'                          following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".
#' @param beta0.mean        a scalar or a vector whose dimension is equal to the number of regression coefficients
#'                          giving the mean parameters for the prior on the historical data regression coefficients. If a
#'                          scalar is provided, `beta0.mean` will be a vector of repeated elements of the given scalar.
#'                          Defaults to a vector of 0s.
#' @param beta0.sd          a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the prior on the historical data regression coefficients. If a scalar is
#'                          provided, same as for `beta0.mean`. Defaults to a vector of 10s.
#' @param p.spike           a scalar between 0 and 1 giving the probability of the spike component in spike-and-slab prior
#'                          on commensurability parameter \eqn{\tau}. Defaults to 0.1.
#' @param spike.mean        a scalar giving the location parameter for the half-normal prior (spike component) on \eqn{\tau}.
#'                          Defaults to 200.
#' @param spike.sd          a scalar giving the scale parameter for the half-normal prior (spike component) on \eqn{\tau}.
#'                          Defaults to 0.1.
#' @param slab.mean         a scalar giving the location parameter for the half-normal prior (slab component) on \eqn{\tau}.
#'                          Defaults to 0.
#' @param slab.sd           a scalar giving the scale parameter for the half-normal prior (slab component) on \eqn{\tau}.
#'                          Defaults to 5.
#' @param scale.mean        location parameter for the half-normal prior on the scale parameters of current and historical
#'                          data models. Defaults to 0.
#' @param scale.sd          scale parameter for the half-normal prior on the scale parameters of current and historical data
#'                          models. Defaults to 10.
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
#'  The function returns an object of class `draws_df` containing posterior samples. The object has two attributes:
#'
#'  \describe{
#'    \item{data}{a list of variables specified in the data block of the Stan program}
#'
#'    \item{model}{a character string indicating the model name}
#'  }
#'
#' @references
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
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
#'     aft.commensurate(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       p.spike = 0.1,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
aft.commensurate = function(
    formula,
    data.list,
    dist              = "weibull",
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  K = length(data.list)
  if ( K == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for CP
  standat = get.aft.stan.data.cp(
    formula        = formula,
    data.list      = data.list,
    dist           = dist,
    beta0.mean     = beta0.mean,
    beta0.sd       = beta0.sd,
    p.spike        = p.spike,
    spike.mean     = spike.mean,
    spike.sd       = spike.sd,
    slab.mean      = slab.mean,
    slab.sd        = slab.sd,
    scale.mean     = scale.mean,
    scale.sd       = scale.sd,
    get.loglik     = get.loglik
  )

  aft_commensurate = instantiate::stan_package_model(
    name = "aft_commensurate",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = aft_commensurate$sample(data = standat,
                                iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X_obs
  oldnames = c(paste0("beta[", 1:p, "]"), "scale", paste0("beta0[", 1:p, "]"), "scale0")
  newnames = c(colnames(X), "scale", paste0( colnames(X), '_hist'), "scale_hist")

  ## reorder parameters so that regression coefficients appear at the top
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  ## add model name as an attribute
  attr(x = d, which = 'model') = "aft_commensurate"
  return(d)
}
