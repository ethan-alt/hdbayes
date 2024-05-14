#' Posterior of normalized asymptotic power prior (NAPP)
#'
#' Sample from the posterior distribution of a GLM using the NAPP by Ibrahim et al. (2015) <doi:10.1002/sim.6728>.
#'
#' The NAPP assumes that the regression coefficients and logarithm of the dispersion parameter are a multivariate
#' normal distribution with mean equal to the maximum likelihood estimate of the historical data and covariance
#' matrix equal to \eqn{a_0^{-1}} multiplied by the inverse Fisher information matrix of the historical data,
#' where \eqn{a_0} is the power prior parameter (treated as random).
#'
#' @include data_checks.R
#' @include get_stan_data.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical datasets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns an object of class `draws_df` giving posterior samples.
#'
#' @references
#'  Ibrahim, J. G., Chen, M., Gwon, Y., and Chen, F. (2015). The power prior: Theory and applications. Statistics in Medicine, 34(28), 3724–3749.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.napp(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data_list,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.napp = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NAPP
  standat = get.stan.data.napp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    offset.list = offset.list,
    a0.shape1   = a0.shape1,
    a0.shape2   = a0.shape2
  )

  glm_napp   = instantiate::stan_package_model(
    name = "glm_napp",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_napp$sample(data = standat,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                        ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X
  K        = standat$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0s[', 1:K, ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:K))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  return(d)
}
