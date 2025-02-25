#' Posterior of latent exchangeability prior (LEAP)
#'
#' Sample from the posterior distribution of a GLM using the latent exchangeability prior (LEAP) by Alt et al. (2024)
#' <doi:10.1093/biomtc/ujae083>.
#'
#' The latent exchangeability prior (LEAP) discounts the historical data by identifying the most relevant individuals
#' from the historical data. It is equivalent to a prior induced by the posterior of a finite mixture model for the
#' historical data set.
#'
#' @include data_checks.R
#' @include get_stan_data.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For LEAP implementation, all historical data sets will be
#'                          stacked into one historical data set.
#' @param K                 the desired number of classes to identify. Defaults to 2.
#' @param prob.conc         a scalar or a vector of length `K` giving the concentration parameters for Dirichlet prior.
#'                          If length == 2, a `Beta(prob.conc[1], prob.conc[2])` prior is used. If a scalar is provided,
#'                          `prob.conc` will be a vector of repeated elements of the given scalar. Defaults to a vector of 1s.
#' @param offset.list       a list of matrices giving the offset for current data followed by historical data. For each
#'                          matrix, the number of rows corresponds to observations and columns correspond to classes.
#'                          Defaults to a list of matrices of 0s. Note that the first element of `offset.list` (corresponding
#'                          to the offset for current data) should be a matrix of repeated columns if `offset.list` is not NULL.
#' @param beta.mean         a scalar or a `p x K` matrix of mean parameters for initial prior on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept). If a scalar is provided,
#'                          `beta.mean` will be a matrix of repeated elements of the given scalar. Defaults to a matrix of 0s.
#' @param beta.sd           a scalar or a `p x K` matrix of sd parameters for the initial prior on regression coefficients,
#'                          where `p` is the number of regression coefficients (including intercept). If a scalar is provided,
#'                          same as for `beta.mean`. Defaults to a matrix of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the location
#'                          parameters for the half-normal priors on the dispersion parameters. If a scalar is provided,
#'                          `disp.mean` will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the scale
#'                          parameters for the half-normal priors on the dispersion parameters. If a scalar is provided, same
#'                          as for `disp.mean`. Defaults to a vector of 10s.
#' @param gamma.lower       a scalar giving the lower bound for probability of subjects in historical data being exchangeable
#'                          with subjects in current data. Defaults to 0.
#' @param gamma.upper       a scalar giving the upper bound for probability of subjects in historical data being exchangeable
#'                          with subjects in current data. Defaults to 1.
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
#' data(actg019)
#' data(actg036)
#' # take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' if (instantiate::stan_cmdstan_exists()) {
#'   glm.leap(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = list(actg019, actg036),
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.leap = function(
    formula,
    family,
    data.list,
    K                 = 2,
    prob.conc         = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    gamma.lower       = 0,
    gamma.upper       = 1,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  if ( length(data.list) == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for LEAP
  standat = get.stan.data.leap(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    K              = K,
    prob.conc      = prob.conc,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    gamma.lower    = gamma.lower,
    gamma.upper    = gamma.upper
  )

  glm_leap = instantiate::stan_package_model(
    name = "glm_leap",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_leap$sample(data = standat,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)

  ## rename and reorder parameters so that regression coefficients appear at the top
  p        = standat$p
  X        = standat$X
  K        = standat$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  oldnames = c(oldnames, paste0( 'probs[', 1:K, ']' ))
  newnames = c(newnames, paste0( 'probs[', 1:K, ']' ))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
