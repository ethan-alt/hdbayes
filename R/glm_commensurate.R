#' Posterior of commensurate prior (CP)
#'
#' Sample from the posterior distribution of a GLM using the commensurate prior (CP) by Hobbs et al. (2011) <doi:10.1111/j.1541-0420.2011.01564.x>.
#'
#' The commensurate prior (CP) assumes that the regression coefficients for the current data conditional on those for
#' the historical data are independent normal distributions with mean equal to the corresponding regression coefficients
#' for the historical data and variance equal to the inverse of the corresponding elements of a vector of precision
#' parameters (referred to as the commensurability parameter \eqn{\tau}). We regard \eqn{\tau} as random and elicit
#' a spike-and-slab prior, which is specified as a mixture of two half-normal priors, on \eqn{\tau}.
#'
#' The number of current data regression coefficients is assumed to be the same as that of historical data regression
#' coefficients. The priors on the dispersion parameters (if applicable) for the current and historical data sets are
#' independent half-normal distributions.
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
#' @param offset.list       a list of vectors giving the offsets for each data. The length of `offset.list` is equal to
#'                          the length of `data.list`. The length of each element of `offset.list` is equal to the number
#'                          of rows in the corresponding element of `data.list`. Defaults to a list of vectors of 0s.
#' @param beta0.mean        a scalar or a vector whose dimension is equal to the number of regression coefficients
#'                          giving the mean parameters for the prior on the historical data regression coefficients. If a
#'                          scalar is provided, `beta0.mean` will be a vector of repeated elements of the given scalar.
#'                          Defaults to a vector of 0s.
#' @param beta0.sd          a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the prior on the historical data regression coefficients. If a scalar is
#'                          provided, same as for `beta0.mean`. Defaults to a vector of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the location parameters for the half-normal priors on the dispersion parameters.
#'                          If a scalar is provided, same as for `beta0.mean`. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the scale parameters for the half-normal priors on the dispersion parameters. If a
#'                          scalar is provided, same as for `beta0.mean`. Defaults to a vector of 10s.
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
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   actg036 = actg036[1:50, ]
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   glm.commensurate(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson(), data.list = data_list,
#'     p.spike = 0.1,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.commensurate = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
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
  standat = get.stan.data.cp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    offset.list    = offset.list,
    beta0.mean     = beta0.mean,
    beta0.sd       = beta0.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    p.spike        = p.spike,
    spike.mean     = spike.mean,
    spike.sd       = spike.sd,
    slab.mean      = slab.mean,
    slab.sd        = slab.sd,
    get.loglik     = get.loglik
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
  X        = standat$X
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"))
  newnames = c(colnames(X), paste0( colnames(X), '_hist') )

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
  }
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  ## add model name as an attribute
  attr(x = d, which = 'model') = "glm_commensurate"
  return(d)
}
