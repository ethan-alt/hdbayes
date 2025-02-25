#' Posterior of latent exchangeability prior (LEAP)
#'
#' Sample from the posterior distribution of an accelerated failure time (AFT) model using the latent exchangeability
#' prior (LEAP) by Alt et al. (2024) <doi:10.1093/biomtc/ujae083>.
#'
#' The latent exchangeability prior (LEAP) discounts the historical data by identifying the most relevant individuals
#' from the historical data. It is equivalent to a prior induced by the posterior of a finite mixture model for the
#' historical data set.
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
#' @param K                 the desired number of classes to identify. Defaults to 2.
#' @param prob.conc         a scalar or a vector of length `K` giving the concentration parameters for Dirichlet prior.
#'                          If length == 2, a `Beta(prob.conc[1], prob.conc[2])` prior is used. If a scalar is provided,
#'                          `prob.conc` will be a vector of repeated elements of the given scalar. Defaults to a vector of 1s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          `beta.mean` will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for `beta.mean`. Defaults to a vector of 10s.
#' @param scale.mean        location parameter for the half-normal prior on the scale parameters for each class. Defaults to 0.
#' @param scale.sd          scale parameter for the half-normal prior on the scale parameters for each class. Defaults to 10.
#' @param gamma.lower       a scalar giving the lower bound for probability of subjects in historical data being exchangeable
#'                          with subjects in current data. Defaults to 0.
#' @param gamma.upper       a scalar giving the upper bound for probability of subjects in historical data being exchangeable
#'                          with subjects in current data. Defaults to 1.
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
#' @references
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2024). LEAP: The latent exchangeability prior for borrowing information from historical data. Biometrics, 80(3).
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
#'     aft.leap(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       K = 2,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
aft.leap = function(
    formula,
    data.list,
    dist              = "weibull",
    K                 = 2,
    prob.conc         = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    gamma.lower       = 0,
    gamma.upper       = 1,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  if ( length(data.list) == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for LEAP
  standat = get.aft.stan.data.leap(
    formula        = formula,
    data.list      = data.list,
    dist           = dist,
    K              = K,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    scale.mean     = scale.mean,
    scale.sd       = scale.sd,
    prob.conc      = prob.conc,
    gamma.lower    = gamma.lower,
    gamma.upper    = gamma.upper,
    get.loglik     = get.loglik
  )

  aft_leap = instantiate::stan_package_model(
    name = "aft_leap",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = aft_leap$sample(data = standat,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                        ...)

  ## rename and reorder parameters so that regression coefficients appear at the top
  p        = standat$p
  X        = standat$X_obs
  K        = standat$K
  oldnames = c(paste0("beta[", 1:p, "]"), "scale")
  newnames = c(colnames(X), "scale")
  oldnames = c(oldnames, paste0( 'probs[', 1:K, ']' ))
  newnames = c(newnames, paste0( 'probs[', 1:K, ']' ))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
