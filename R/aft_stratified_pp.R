#' Posterior of stratified power prior (PP) with fixed \eqn{a_0}
#'
#' Sample from the posterior distribution of an accelerated failure time (AFT) model using the power prior (PP) within
#' predefined strata. If the strata and power prior parameters (\eqn{a_0}'s) are determined based on propensity scores,
#' this function can be used to sample from the posterior of an AFT model under the propensity score-integrated power
#' prior (PSIPP) by Wang et al. (2019) <doi:10.1080/10543406.2019.1657133>.
#'
#' The power prior parameters (\eqn{a_0}'s) are treated as fixed and must be provided for each stratum. Users must
#' also specify the stratum ID for each observation. Within each stratum, the initial priors on the regression
#' coefficients are independent normal priors, and the current and historical data models are assumed to have a common
#' scale parameter with a half-normal prior.
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
#' @param strata.list       a list of vectors specifying the stratum ID for each observation in the corresponding data set
#'                          in `data.list`. The first element in the list corresponds to the current data, and the rest
#'                          correspond to the historical data sets. Each vector should have the same length as the number
#'                          of rows in the respective data set in `data.list`, with values representing stratum labels
#'                          as positive integers (e.g., 1, 2, 3, ...).
#' @param a0.strata         A scalar or a vector of fixed power prior parameters (\eqn{a_0}'s) for each stratum, with values
#'                          between 0 and 1. If a scalar is provided, it will be replicated for all strata. If a vector is
#'                          provided, its length must match the total number of unique strata across all data sets. The first
#'                          element of `a0.strata` corresponds to stratum 1, the second to stratum 2, and so on.
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
#'  The function returns an object of class `draws_df` giving posterior samples, with an attribute called 'data' which includes
#'  the list of variables specified in the data block of the Stan program.
#'
#' @references
#'  Wang, C., Li, H., Chen, W.-C., Lu, N., Tiwari, R., Xu, Y., & Yue, L. Q. (2019). Propensity score-integrated power prior approach for incorporating real-world evidence in single-arm clinical studies. Journal of Biopharmaceutical Statistics, 29(5), 731–748.
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
#'     data_list = list(currdata = E1690, histdata = E1684)
#'     strata_list = list(rep(1:2, each = 25), rep(1:2, each = 50))
#'     # Alternatively, we can determine the strata based on propensity scores
#'     # using the psrwe package, which is available on GitHub
#'     aft.stratified.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment,
#'       data.list = data_list,
#'       strata.list = strata_list,
#'       a0.strata = c(0.3, 0.5),
#'       dist = "weibull",
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
aft.stratified.pp = function(
    formula,
    data.list,
    strata.list,
    a0.strata,
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
  ## get Stan data for stratified PP
  standat = get.aft.stan.data.stratified.pp(
    formula     = formula,
    data.list   = data.list,
    strata.list = strata.list,
    a0.strata   = a0.strata,
    dist        = dist,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    scale.mean  = scale.mean,
    scale.sd    = scale.sd,
    get.loglik  = get.loglik
  )

  aft_stratified_pp = instantiate::stan_package_model(
    name = "aft_stratified_pp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = aft_stratified_pp$sample(data = standat,
                                 iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                 ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X_obs
  K        = standat$K
  oldnames = c( paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"),
                paste0("scaleVec[", 1:K, "]") )
  newnames = c( paste0( colnames(X), '_stratum_', rep(1:K, each = p) ),
                paste0("scale_stratum_", 1:K) )

  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
