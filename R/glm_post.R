#' Posterior of a normal/half-normal prior
#'
#' Sample from the posterior distribution of a GLM using a normal/half-normal prior.
#'
#' The priors on the regression coefficients are independent normal distributions. When the normal priors are elicited
#' with large variances, the prior is also referred to as the reference or vague prior. The dispersion parameter is
#' assumed to be independent of the regression coefficients with a half-normal prior (if applicable).
#'
#' @include data_checks.R
#' @include get_stan_data.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list consisting of one `data.frame` giving the current data. If `data.list` has more
#'                          than one `data.frame`, only the first element will be used as the current data.
#' @param offset.list       a list consisting of one vector giving the offset for the current data. The length of
#'                          the vector is equal to the number of rows in the current data. The vector has all values
#'                          set to 0 by default. If `offset.list` has more than one vector, same as for `data.list`.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the normal prior on regression coefficients. If a scalar is provided,
#'                          `beta.mean` will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the normal prior on regression coefficients. If a scalar is provided,
#'                          same as for `beta.mean`. Defaults to a vector of 10s.
#' @param disp.mean         location parameter for the half-normal prior on dispersion parameter. Defaults to 0. If
#'                          `disp.mean` is a vector with length > 1, only the first element will be used as `disp.mean`.
#' @param disp.sd           scale parameter for the half-normal prior on dispersion parameter. Defaults to 10. If
#'                          `disp.sd` is a vector with length > 1, same as for `disp.mean`.
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
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   ## take subset for speed purposes
#'   actg019 = actg019[1:100, ]
#'   data.list = list(currdata = actg019)
#'   glm.post(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     data.list = data.list,
#'     beta.sd   = 10,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.post = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for normal/half-normal prior
  standat = get.stan.data.post(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    offset.list = offset.list,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd
  )

  glm_post      = instantiate::stan_package_model(
    name = "glm_post",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_post$sample(data = standat,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                      ...)

  ## rename parameters
  p        = standat$p
  X        = standat$X
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  d = rename.params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standat
  return(d)
}
