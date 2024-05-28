#' Posterior of commensurate prior (CP)
#'
#' Sample from the posterior distribution of a GLM using the CP by Hobbs et al. (2011) <doi:10.1111/j.1541-0420.2011.01564.x>.
#'
#' The CP assumes that the regression coefficients for the current data conditional on those for the historical
#' data are independent normal distributions with mean equal to the corresponding regression coefficients
#' for the historical data and variance equal to the inverse of the corresponding elements of a user-specified
#' vector (tau) of precision parameters. The number of current data regression coefficients is assumed to be the
#' same as that of historical data regression coefficients. The priors on the dispersion parameters (if applicable)
#' for the current and historical data sets are independent half-normal distributions.
#'
#' @include data_checks.R
#' @include get_stan_data.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param tau               a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the commensurate prior parameters. If a scalar is provided, tau will be a vector of repeated
#'                          elements of the given scalar. Each element of tau must be positive, corresponding to a normal
#'                          precision parameter.
#' @param beta0.mean        a scalar or a vector whose dimension is equal to the number of regression coefficients
#'                          giving the mean parameters for the prior on the historical data regression coefficients. If a
#'                          scalar is provided, same as for tau. Defaults to a vector of 0s.
#' @param beta0.sd          a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the prior on the historical data regression coefficients. If a scalar is
#'                          provided, same as for tau. Defaults to a vector of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the means for the half-normal priors on the dispersion parameters. If a scalar is
#'                          provided, same as for tau. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the sds for the half-normal priors on the dispersion parameters. If a scalar is
#'                          provided, same as for tau. Defaults to a vector of 10s.
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
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.commensurate(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson(), data.list = data_list,
#'     tau = rep(5, 4),    ## 4 parameters including intercept
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.commensurate = function(
    formula,
    family,
    data.list,
    tau,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for CP
  standat = get.stan.data.cp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    tau            = tau,
    offset.list    = offset.list,
    beta0.mean     = beta0.mean,
    beta0.sd       = beta0.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd
  )

  glm_commensurate = instantiate::stan_package_model(
    name = "glm_commensurate",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_commensurate$sample(data = standat,
                                iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                                ...)

  ## rename parameters
  p        = standat$p
  K        = standat$K
  X        = standat$X
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"))
  newnames = c(colnames(X), paste0( colnames(X), '_hist') )

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
  }
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  return(d)
}
