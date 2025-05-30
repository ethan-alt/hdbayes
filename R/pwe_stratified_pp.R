#' Posterior of stratified power prior (PP) with fixed \eqn{a_0}
#'
#' Sample from the posterior distribution of a piecewise exponential (PWE) model (i.e., a proportional hazards model
#' with a piecewise constant baseline hazard) using the power prior (PP) within predefined strata. If the strata and
#' power prior parameters (\eqn{a_0}'s) are determined based on propensity scores, this function can be used to
#' sample from the posterior of a PWE model under the propensity score-integrated power prior (PSIPP) by Wang et al.
#' (2019) <doi:10.1080/10543406.2019.1657133>.
#'
#' The power prior parameters (\eqn{a_0}'s) are treated as fixed and must be provided for each stratum. Users must
#' also specify the stratum ID for each observation. Within each stratum, the initial priors on the regression
#' coefficients are independent normal priors, and the current and historical data models are assumed to share the
#' baseline hazard parameters with half-normal priors.
#'
#' @include data_checks_pwe.R
#' @include get_stan_data_pwe.R
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
#' @param strata.list       a list of vectors specifying the stratum ID for each observation in the corresponding data set
#'                          in `data.list`. The first element in the list corresponds to the current data, and the rest
#'                          correspond to the historical data sets. Each vector should have the same length as the number
#'                          of rows in the respective data set in `data.list`, with values representing stratum labels
#'                          as positive integers (e.g., 1, 2, 3, ...).
#' @param breaks            a numeric vector specifying the time points that define the boundaries of the piecewise
#'                          intervals. The values should be in ascending order, with the final value being greater than
#'                          or equal to the maximum observed time.
#' @param a0.strata         A scalar or a vector of fixed power prior parameters (\eqn{a_0}'s) for each stratum, with values
#'                          between 0 and 1. If a scalar is provided, it will be replicated for all strata. If a vector is
#'                          provided, its length must match the total number of unique strata across all data sets. The first
#'                          element of `a0.strata` corresponds to stratum 1, the second to stratum 2, and so on.
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     pwe.stratified.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment,
#'       data.list = data_list,
#'       strata.list = strata_list,
#'       breaks = breaks,
#'       a0.strata = c(0.3, 0.5),
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'   }
#' }
pwe.stratified.pp = function(
    formula,
    data.list,
    strata.list,
    breaks,
    a0.strata,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for stratified PP
  standat = get.pwe.stan.data.stratified.pp(
    formula          = formula,
    data.list        = data.list,
    strata.list      = strata.list,
    breaks           = breaks,
    a0.strata        = a0.strata,
    beta.mean        = beta.mean,
    beta.sd          = beta.sd,
    base.hazard.mean = base.hazard.mean,
    base.hazard.sd   = base.hazard.sd,
    get.loglik       = get.loglik
  )

  pwe_stratified_pp = instantiate::stan_package_model(
    name = "pwe_stratified_pp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = pwe_stratified_pp$sample(data = standat,
                                 iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                 ...)

  ## rename parameters
  p        = standat$p
  X1       = standat$X1
  J        = standat$J
  K        = standat$K
  if( p > 0 ){
    oldnames = c(paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"),
                 paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]"))
    newnames = c(paste0( colnames(X1), '_stratum_', rep(1:K, each = p) ),
                 paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"))
  }else{
    oldnames = c(paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]"))
    newnames = c(paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"))
  }

  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
