#'
#' Final step for sampling from the posterior distribution of a GLM using the Robust Meta-Analytic Predictive (MAP)
#' Prior Of Schmidli. See `glm_rmap_bhm.R` for the first step and `glm_rmap_bhm_approx.R` for the second step.
#'
#' This function samples from the posterior distribution of a GLM using the Robust MAP prior. The first component of
#' the Robust MAP prior is a prior induced by the Bayesian Hierarchical Model (BHM). We approximate this component by
#' a mixture of multivariate normal distributions where the parameters are obtained from the outputs of the
#' `glm.rmap.bhm.approx()` function. The second component is a vague (noninformative) multivariate normal prior. We assume
#' that the covariance matrix of the vague prior is a diagonal matrix.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param curr.data         a `data.frame` giving the current data.
#' @param include.intercept logical; if TRUE, an intercept will be included in the model.
#' @param curr.offset       a vector whose dimension is equal to the rows of the current dataset giving an offset for
#'                          the current data. Defaults to a vector of 0s.
#' @param probs             a vector of mixing proportions in the mixture approximation to the prior induced by the BHM.
#'                          Obtained from the outputs of the `glm.rmap.bhm.approx()` function.
#' @param means             a matrix with the jth column being the mean vector for the jth component in the mixture
#'                          approximation to the prior induced by the BHM. Obtained from the outputs of the
#'                          `glm.rmap.bhm.approx()` function.
#' @param covs              a 3-dimensional array giving the covariance matrices for the mixture approximation to the prior
#'                          induced by the BHM. Obtained from the outputs of the `glm.rmap.bhm.approx()` function.
#'                          the means for the half-normal hyperpriors on the sd hyperparameters of regression coefficients.
#' @param w                 a scalar between 0 and 1 giving how much weight to put on the historical data.
#' @param norm.vague.mean   a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the means for the vague normal prior on regression coefficients. If a scalar is provided,
#'                          norm.vague.mean will be a vector of repeated elements of the given scalar. Defaults to a
#'                          vector of 0s.
#' @param norm.vague.sd     a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sds for the vague normal prior on regression coefficients. If a scalar is provided, same as
#'                          for norm.vague.mean. Defaults to a vector of 10s.
#' @param curr.disp.mean    a scalar giving the mean for the half-normal hyperprior on the dispersion parameter for the current
#'                          data. Defaults to a vector of 0s.
#' @param curr.disp.sd      a scalar giving the sd for the half-normal hyperprior on the dispersion parameter for the current
#'                          data. Defaults to a vector of 10s.
#' @param local.location    a file path giving the desired location of the local copies of all the .stan model files in the
#'                          package. Defaults to the path created by `rappdirs::user_cache_dir("hdbayes")`.
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. iter_warmup, iter_sampling, chains).
#'
#' @return                  an object of class `draws_df` giving posterior samples
#'
#' @examples
#' data(actg019) ## current data
#' data(actg036) ## historical data
#' ## take subset for speed purposes
#' actg019 = actg019[1:150, ]
#' actg036 = actg036[1:100, ]
#' hist_data_list = list(actg036)
#' samples_bhm = glm.rmap.bhm(
#'   formula = cd4 ~ treatment + age + race,
#'   family = poisson('log'),
#'   hist.data.list = hist_data_list,
#'   chains = 1, iter_warmup = 1000, iter_sampling = 2000
#' )$beta_pred
#' res_approx = glm.rmap.bhm.approx(
#'   samples.bhm = samples_bhm,
#'   G = 1:5, verbose = FALSE
#' )
#' glm.rmap(
#'   formula = cd4 ~ treatment + age + race,
#'   family = poisson('log'),
#'   curr.data = actg019,
#'   probs = res_approx$probs,
#'   means = res_approx$means,
#'   covs = res_approx$covs,
#'   chains = 1, iter_warmup = 1000, iter_sampling = 2000
#' )
glm.rmap = function(
    formula,
    family,
    curr.data,
    probs,
    means,
    covs,
    include.intercept = TRUE,
    curr.offset       = NULL,
    w                 = 0.1,
    norm.vague.mean   = NULL,
    norm.vague.sd     = NULL,
    curr.disp.mean    = NULL,
    curr.disp.sd      = NULL,
    local.location    = NULL,
    ...
) {
  ## perform data checks
  if ( !is.null(curr.offset) ){
    data.checks(formula, family, list(curr.data), list(curr.offset))
  }else {
    data.checks(formula, family, list(curr.data), curr.offset)
  }

  y        = curr.data[, all.vars(formula)[1]]
  n        = length(y)
  X        = model.matrix(formula, curr.data)
  p        = ncol(X)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]
  G        = length(probs)

  if ( is.null(curr.offset) )
    curr.offset = rep(0, n)

  ## Default half-normal hyperprior on dispersion parameter (if exist) is N^{+}(0, 10^2)
  if ( is.null(curr.disp.mean) )
    curr.disp.mean = 0
  if ( is.null(curr.disp.sd) )
    curr.disp.sd = 10

  ## Default vague prior is indepndent N(0, 10^2)
  if ( !is.null(norm.vague.mean) ){
    if ( !( is.vector(norm.vague.mean) & (length(norm.vague.mean) %in% c(1, p)) ) )
      stop("norm.vague.mean must be a scalar or a vector of length ", p, " if norm.vague.mean is not NULL")
  }
  norm.vague.mean = to.vector(param = norm.vague.mean, default.value = 0, len = p)
  if ( !is.null(norm.vague.sd) ){
    if ( !( is.vector(norm.vague.sd) & (length(norm.vague.sd) %in% c(1, p)) ) )
      stop("norm.vague.sd must be a scalar or a vector of length ", p, " if norm.vague.sd is not NULL")
  }
  norm.vague.sd = to.vector(param = norm.vague.sd, default.value = 10, len = p)


  standat = list(
    'n1'              = n,
    'p'               = p,
    'y1'              = y,
    'X1'              = X,
    'G'               = G,
    'probs'           = probs,
    'means'           = means,
    'covars'          = covs,
    'norm_vague_mean' = norm.vague.mean,
    'norm_vague_sd'   = norm.vague.sd,
    'disp_mean1'      = curr.disp.mean,
    'disp_sd1'        = curr.disp.sd,
    'w'               = w,
    'dist'            = dist,
    'link'            = link,
    'offs1'           = curr.offset
  )

  ## copy all the .stan model files to the specified local location
  if( is.null(local.location) )
    local.location <- rappdirs::user_cache_dir(appname = "hdbayes")

  if (length(list.files(local.location, pattern = ".stan")) >= 1) {
    cli::cli_alert_info("Using cached Stan models")
  } else {
    cli::cli_alert_info("Copying Stan models to cache")
    staninside::copy_models(pkgname = "hdbayes",
                            local_location = local.location)
    cli::cli_alert_success("Models copied!")
  }

  model_name      = "glm_robustmap"
  model_file_path = file.path(local.location, paste0(model_name, ".stan"))
  glm_robustmap   = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = glm_robustmap$sample(data = standat, ...)
  d   = fit$draws(format = 'draws_df')

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
