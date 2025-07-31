#' Posterior of a normal/half-normal prior
#'
#' Sample from the posterior distribution of an accelerated failure time (AFT) model using a normal/half-normal prior.
#'
#' The priors on the regression coefficients are independent normal distributions. When the normal priors are elicited
#' with large variances, the prior is also referred to as the reference or vague prior. The scale parameter is assumed
#' to be independent of the regression coefficients with a half-normal prior.
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
#' @param data.list         a list consisting of one `data.frame` giving the current data. If `data.list` has more
#'                          than one `data.frame`, only the first element will be used as the current data.
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
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1690)
#'     ## take subset for speed purposes
#'     E1690 = E1690[1:100, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'     E1690$cage = as.numeric(scale(E1690$age))
#'     data_list = list(currdata = E1690)
#'     aft.post(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       beta.sd = 10,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
aft.post = function(
    formula,
    data.list,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for normal/half-normal prior
  standat = get.aft.stan.data.post(
    formula     = formula,
    data.list   = data.list,
    dist        = dist,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    scale.mean  = scale.mean,
    scale.sd    = scale.sd,
    get.loglik  = get.loglik
  )

  aft_post      = instantiate::stan_package_model(
    name = "aft_post",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = aft_post$sample(data = standat,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                      ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X_obs
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  ## add model name as an attribute
  attr(x = d, which = 'model') = "aft_post"
  return(d)
}
