#'
#' Commensurate prior
#'
#' Sample from the posterior distribution of a GLM using the Commensurate
#' prior of Hobbs et al. This prior assumes that the regression coefficients
#' for the current data set conditional on those for the historical data set are
#' multivariate normal with mean equal to the regression coefficients of the historical data
#' and covariance equal to the inverse of a user-specified precision parameter (tau).
#' Dispersion parameters (if applicable) between the
#' current and historical data sets are treated as independent.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rests
#'                          are the historical datasets.
#' @param include.intercept logical; if TRUE, an intercept will be included in the model. Defaults to TRUE.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param tau               a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the commensurate prior parameters. If a scalar is provided, tau will be a vector of repeated
#'                          elements of the given scalar. Each element of tau must be positive, corresponding to a normal
#'                          precision parameter.
#' @param beta0.mean        a scalar or a vector whose dimension is equal to the number of regression coefficients
#'                          giving the mean parameters for the initial prior on regression coefficients. If a scalar
#'                          is provided, same as for tau. Defaults to a vector of 0s.
#' @param beta0.sd          a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for tau. Defaults to a vector of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of datasets (including the current
#'                          data) giving the means for the half-normal hyperpriors on the dispersion parameters. If a
#'                          scalar is provided, same as for tau. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of datasets (including the current
#'                          data) giving the sds for the half-normal hyperpriors on the dispersion parameters. If a scalar
#'                          is provided, same as for tau. Defaults to a vector of 10s.
#' @param local.location    a file path giving the desired location of the local copies of all the .stan model files in
#'                          the package. Defaults to the path created by `rappdirs::user_cache_dir("hdbayes")`.
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. iter_warmup, iter_sampling, chains).
#'
#' @return                  an object of class `draws_df` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' ## take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' data_list = list(currdata = actg019, histdata = actg036)
#' glm.commensurate(
#'   formula = cd4 ~ treatment + age + race,
#'   family = poisson(), data.list = data_list,
#'   tau = rep(5, 4),    ## 4 parameters including intercept
#'   chains = 1, iter_warmup = 500, iter_sampling = 1000
#' )
glm.commensurate = function(
    formula,
    family,
    data.list,
    tau,
    include.intercept = TRUE,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    local.location    = NULL,
    ...
) {
  data.checks(formula, family, data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list,
                            include.intercept = include.intercept)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  ## Default offset for each dataset is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  ## check tau
  if ( !( is.vector(tau) & (length(tau) %in% c(1, p)) ) )
    stop("tau must be a scalar or a vector of length ", p)
  tau = to.vector(param = tau, len = p)
  ## Check if each element of tau is positive
  if ( any(is.na(tau) ) )
    stop('tau must be a vector of non-missing values')
  if ( any(tau <= 0) )
    stop("Each element of tau must be positive")

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta0.mean) ){
    if ( !( is.vector(beta0.mean) & (length(beta0.mean) %in% c(1, p)) ) )
      stop("beta0.mean must be a scalar or a vector of length ", p, " if beta0.mean is not NULL")
  }
  beta0.mean = to.vector(param = beta0.mean, default.value = 0, len = p)
  if ( !is.null(beta0.sd) ){
    if ( !( is.vector(beta0.sd) & (length(beta0.sd) %in% c(1, p)) ) )
      stop("beta0.sd must be a scalar or a vector of length ", p, " if beta0.sd is not NULL")
  }
  beta0.sd = to.vector(param = beta0.sd, default.value = 10, len = p)

  ## Default half-normal hyperprior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) %in% c(1, K)) ) )
      stop("disp.mean must be a scalar or a vector of length ", K, " if disp.mean is not NULL")
  }
  disp.mean = to.vector(param = disp.mean, default.value = 0, len = K)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) %in% c(1, K)) ) )
      stop("disp.sd must be a scalar or a vector of length ", K, " if disp.sd is not NULL")
  }
  disp.sd = to.vector(param = disp.sd, default.value = 10, len = K)


  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'beta0_mean'      = beta0.mean,
    'beta0_sd'        = beta0.sd,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'tau'             = tau,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
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

  model_name       = "glm_commensurate"
  model_file_path  = file.path(local.location, paste0(model_name, ".stan"))
  glm_commensurate = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = glm_commensurate$sample(data = standat, ...)
  d   = fit$draws(format = 'draws_df')


  ## rename parameters
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"))
  newnames = c(colnames(X), paste0( colnames(X), '_hist') )

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
    newnames = c(newnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
  }
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
