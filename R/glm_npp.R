#' Posterior of normalized power prior (NPP)
#'
#' Sample from the posterior distribution of a GLM using the NPP by Duan et al. (2006) <doi:10.1002/env.752>.
#'
#' Before using this function, users must estimate the logarithm of the normalizing constant across a
#' range of different values for the power prior parameter (\eqn{a_0}), possibly smoothing techniques
#' over a find grid. The power prior parameters (\eqn{a_0}'s) are treated as random with independent
#' beta priors. The initial priors on the regression coefficients are independent normal priors. The
#' current and historical data sets are assumed to have a common dispersion parameter with a
#' half-normal prior (if applicable). For normal linear models, the exact normalizing constants for
#' NPP can be computed. See the implementation in [lm.npp()].
#'
#'
#' @include data_checks.R
#' @include glm_npp_lognc.R
#' @include get_stan_data.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for beta.mean. Defaults to a vector of 10s.
#' @param disp.mean         mean parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           sd parameter for the half-normal prior on dispersion parameter. Defaults to 10.
#' @param a0.lognc          a vector giving values of the power prior parameter for which the logarithm of the normalizing
#'                          constant has been evaluated.
#' @param lognc             an S by T matrix where S is the length of a0.lognc, T is the number of historical data sets, and
#'                          the j-th column, j = 1, ..., T, is a vector giving the logarithm of the normalizing constant (as
#'                          estimated by [glm.npp.lognc()] for a0.lognc using the j-th historical data set.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.lower          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          lower bounds for each element of the a0 vector. If a scalar is provided, a0.lower will be a
#'                          vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param a0.upper          a scalar or a vector whose dimension is equal to the number of historical data sets giving the
#'                          upper bounds for each element of the a0 vector. If a scalar is provided, same as for a0.lower.
#'                          Defaults to a vector of 1s.
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
#' @seealso [glm.npp.lognc()]
#'
#' @references
#'  Duan, Y., Ye, K., and Smith, E. P. (2005). Evaluating water quality using power priors to incorporate historical information. Environmetrics, 17(1), 95–106.
#'
#' @examples
#' \donttest{
#'   if(requireNamespace("parallel")){
#'     data(actg019)
#'     data(actg036)
#'     ## take subset for speed purposes
#'     actg019 = actg019[1:100, ]
#'     actg036 = actg036[1:50, ]
#'
#'     library(parallel)
#'     ncores    = 2
#'     data.list = list(data = actg019, histdata = actg036)
#'     formula   = cd4 ~ treatment + age + race
#'     family    = poisson()
#'     a0        = seq(0, 1, length.out = 11)
#'     if (instantiate::stan_cmdstan_exists()) {
#'       ## call created function
#'       ## wrapper to obtain log normalizing constant in parallel package
#'       logncfun = function(a0, ...){
#'         hdbayes::glm.npp.lognc(
#'           formula = formula, family = family, a0 = a0, histdata = data.list[[2]],
#'           ...
#'         )
#'       }
#'
#'       cl = makeCluster(ncores)
#'       clusterSetRNGStream(cl, 123)
#'       clusterExport(cl, varlist = c('formula', 'family', 'data.list'))
#'       a0.lognc = parLapply(
#'         cl = cl, X = a0, fun = logncfun, iter_warmup = 500,
#'         iter_sampling = 1000, chains = 1, refresh = 0
#'       )
#'       stopCluster(cl)
#'       a0.lognc = data.frame( do.call(rbind, a0.lognc) )
#'
#'       ## sample from normalized power prior
#'       glm.npp(
#'         formula = formula,
#'         family = family,
#'         data.list = data.list,
#'         a0.lognc = a0.lognc$a0,
#'         lognc = matrix(a0.lognc$lognc, ncol = 1),
#'         chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'         refresh = 0
#'       )
#'     }
#'   }
#' }
glm.npp = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NPP
  standat = get.stan.data.npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta.mean      = beta.mean,
    beta.sd        = beta.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper
  )

  glm_npp_posterior = instantiate::stan_package_model(
    name = "glm_npp_posterior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_npp_posterior$sample(data = standat,
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
  oldnames = c(oldnames, paste0('a0s[', 1:(K-1), ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:(K-1)))
  d        = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  return(d)
}
