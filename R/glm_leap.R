#'
#' Posterior of LEAP
#'
#' ...
#'
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data.list   a list of `data.frame` giving the current data followed by historical data
#' @param K           the desired number of classes to identify
#' @param beta.mean   a `p x K` matrix of mean parameter for initial prior on regression coefficients (including intercept). Defaults to a matrix of zeros.
#' @param beta.cov    a list of `K` `data.frame` each size `p x p` giving covariance parameter for initial prior on regression coefficients (including intercept). Each defaults to a diagonal covariance matrix where each variance is equal to 100.
#' @param disp.shape  shape parameter for inverse-gamma prior on dispersion parameter
#' @param disp.scale  scale parameter for inverse-gamma prior on dispersion parameter
#' @param offset.list a list of `data.frame` giving the offset for current data followed by historical data. For each `data.frame`, the number of rows correpond to observations and columns correspond to classes. Defaults to matrices of 0s
#' @param local.location    a file path giving the desired location of the local copies of all the .stan model files in the
#'                          package. Defaults to the path created by `rappdirs::user_cache_dir("hdbayes")`.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          [cmdstanr::sample()].
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in [cmdstanr::sample()].
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in [cmdstanr::sample()].
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. seed, refresh, init).
#'
#' @return                  an object of class `draws_df` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' # take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' glm.leap(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(),
#'   data.list = list(actg019, actg036),
#'   K = 2,
#'   chains = 1, iter_warmup = 500, iter_sampling = 1000
#' )
#'
glm.leap = function(
    formula,
    family,
    data.list,
    K           = 2,
    beta.mean   = NULL,
    beta.cov    = NULL,
    disp.shape  = 2.1,
    disp.scale  = 1.1,
    offset.list = NULL,
    local.location    = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks(formula, family, data.list, offset.list)

  ## get model information
  data=data.list[[1]]
  histdata=data.list[[2]] # not yet support multiple histdata; only read 2nd dataset
  offset=offset.list[[1]]
  offset0=offset.list[[2]]

  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]
  n  = length(y)
  n0 = length(y0)
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)
  p  = ncol(X)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  ## Default offset is matrix of 0s
  if ( is.null(offset) )
    offset = matrix(rep(0, n*K), ncol=K)
  if ( is.null(offset0) )
    offset0 = matrix(rep(0, n0*K), ncol=K)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) )
    beta.mean = matrix(rep(0, ncol(X)*K), ncol=K)
  if ( is.null(beta.cov) )
    beta.cov  = replicate(K, diag(100, ncol(X)), simplify=F)

  standat = list(
    'n'           = n,
    'n0'          = n0,
    'p'           = p,
    'K'           = K,
    'y'           = y,
    'X'           = X,
    'y0'          = y0,
    'X0'          = X0,
    'mean_beta'   = beta.mean,
    'cov_beta'    = do.call(rbind, beta.cov),
    'disp_shape'  = disp.shape,
    'disp_scale'  = disp.scale,
    'conc'        = c(0.95, 0.95),
    'gamma_lower' = 0,
    'gamma_upper' = 1,
    'dist'        = dist,
    'link'        = link,
    'offs'        = offset,
    'offs0'       = offset0
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

  model_name      = "glm_leap"
  model_file_path = file.path(local.location, paste0(model_name, ".stan"))
  glm_leap        = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = glm_leap$sample(data = standat,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)
  d   = fit$draws(format = 'draws_df')

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
