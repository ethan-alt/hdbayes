#' Posterior of latent exchangeability prior (LEAP) with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the latent exchangeability prior (LEAP) by Alt et al. (
#' 2024) <doi:10.1093/biomtc/ujae083>, extended to handle missing covariates via a Gaussian copula model.
#'
#' The latent exchangeability prior (LEAP) discounts the historical data by identifying the most relevant individuals
#' from the historical data. It is equivalent to a prior induced by the posterior of a finite mixture model for the
#' historical data set.
#'
#' This function extends LEAP to handle missing covariate values using a Gaussian copula approach. Missing values
#' are imputed by sampling from the conditional distribution implied by the copula, with the correlation matrix
#' either estimated from complete cases or provided by the user.
#'
#' @include data_checks_mc.R
#' @include get_stan_data_mc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For LEAP implementation, all historical data sets will be
#'                          stacked into one historical data set. Covariates may contain missing values (NA), which
#'                          will be imputed using the Gaussian copula model.
#' @param K                 the desired number of classes to identify. Defaults to 2.
#' @param prob.conc         a scalar or a vector of length `K` giving the concentration parameters for Dirichlet prior.
#'                          If length == 2, a `Beta(prob.conc[1], prob.conc[2])` prior is used. If a scalar is provided,
#'                          `prob.conc` will be a vector of repeated elements of the given scalar. Defaults to a vector of 1s.
#' @param copula.corr       optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases across all datasets. Must be a q x q positive definite
#'                          matrix where q is the number of covariates with missing values.
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
#'  The function returns an object of class `draws_df` containing posterior samples, including the imputed covariate
#'  values. The object has two attributes:
#'
#'  \describe{
#'    \item{data}{a list of variables specified in the data block of the Stan program}
#'
#'    \item{model}{a character string indicating the model name}
#'  }
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
#' # introduce missing values in the historical data for demonstration
#' actg036$age[sample(1:50, 5)] = NA
#' if (instantiate::stan_cmdstan_exists()) {
#'   glm.leap.mc(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = list(actg019, actg036),
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.leap.mc = function(
    formula,
    family,
    data.list,
    K                 = 2,
    prob.conc         = NULL,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
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

  ## get Stan data for LEAP with missing covariates
  standat = get.stan.data.leap.mc(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    K              = K,
    prob.conc      = prob.conc,
    copula.corr    = copula.corr,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    gamma.lower    = gamma.lower,
    gamma.upper    = gamma.upper,
    get.loglik     = get.loglik
  )

  glm_leap_mc = instantiate::stan_package_model(
    name = "glm_leap_mc",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_leap_mc$sample(data = standat,
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
  ## add model name as an attribute
  attr(x = d, which = 'model') = "glm_leap_mc"
  return(d)
}
