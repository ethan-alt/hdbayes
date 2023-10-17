#'
#' Posterior of normalized power prior
#'
#' Sample from the posterior distribution of a normal linear model
#' using the normalized power prior (NPP)
#'
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical datasets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients.
#'                          Conditional on the variance parameter sigmasq for the outcome, beta.sd * sqrt(sigmasq) gives
#'                          the sd for the initial prior on regression coefficients. If a scalar is provided, same as for
#'                          beta.mean. Defaults to a vector of 10s.
#' @param sigmasq.shape     shape parameter for inverse-gamma prior on variance parameter. Defaults to 2.1.
#' @param sigmasq.scale     scale parameter for inverse-gamma prior on variance parameter. Defaults to 1.1.
#' @param a0.shape1         first shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param a0.shape2         second shape parameter for the i.i.d. beta prior on a0 vector. When \code{a0.shape1 == 1} and
#'                          \code{a0.shape2 == 1}, a uniform prior is used.
#' @param local.location    a file path giving the desired location of the local copies of all the .stan model files in the
#'                          package. Defaults to the path created by `rappdirs::user_cache_dir("hdbayes")`.
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. iter_warmup, iter_sampling, chains).
#'
#' @return                  an object of class `draws_df` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' lm.npp(
#'    formula = cd4 ~ treatment + age + race,
#'    data.list = list(data = actg019, histdata = actg036), chains = 1,
#'    iter_warmup = 500, iter_sampling = 1000, refresh = 0
#' )
#'
#'
lm.npp = function(
    formula,
    data.list,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    sigmasq.shape     = 2.1,
    sigmasq.scale     = 1.1,
    a0.shape1         = 1,
    a0.shape2         = 1,
    local.location    = NULL,
    ...
) {
  data.checks(formula, gaussian('identity'), data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  num.obs      = 1 + end.index - start.index
  p            = ncol(X)
  N            = length(y)
  K            = length(data.list) - 1 # number of historical datasets

  ## Default offset for each dataset is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  ## outcome after adjusting for offset
  y_offset = y - offset

  ## Default prior on regression coefficients given sigmasq is N(0, sigmasq * 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = to.vector(param = beta.sd, default.value = 10, len = p)
  beta.cov = diag(beta.sd^2, nrow = p, ncol = p)

  hist.mle  = matrix(NA, nrow = p, ncol = K)
  hist.prec = array(NA, c(K, p, p)) # each element is X0'X0 for a historical dataset
  sumy0sq   = vector(length = K) # each element is sum(y0^2) for a historical dataset (y0 has been adjusted for offsets)
  for(k in 1:K) {
    histdata       = data.list[[1+k]]
    offs0          = offset[ start.index[1+k]:end.index[1+k] ]
    histdata$offs0 = offs0
    sumy0sq[k]     = sum( y_offset[ start.index[1+k]:end.index[1+k] ]^2 )
    X0             = model.matrix(formula, histdata)
    hist.prec[k, , ] = crossprod(X0)
    fit.lm =
      suppressWarnings(
        lm(formula = formula, data = histdata, offset = offs0)
      )
    hist.mle[, k] = as.numeric(fit.lm$coefficients)
  }


  standat = list(
    'K'               = K, # number of historical datasets
    'n'               = num.obs[1], # current data sample size
    'n0s'             = num.obs[-1],
    'p'               = p,
    'y'               = y_offset[start.index[1]:end.index[1]], # current data outcome after adjusting for offsets
    'X'               = X[start.index[1]:end.index[1], ], # current data design matrix
    "hist_prec"       = hist.prec,
    "hist_mle"        = hist.mle,
    "sumy0sq"         = sumy0sq,
    'mean_beta'       = beta.mean,
    'cov_beta'        = beta.cov,
    'sigmasq_shape'   = sigmasq.shape,
    'sigmasq_scale'   = sigmasq.scale,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2
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

  model_name      = "lm_npp"
  model_file_path = file.path(local.location, paste0(model_name, ".stan"))
  lm_npp          = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = lm_npp$sample(data = standat, ...)
  d   = fit$draws(format = 'draws_df')

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  oldnames = c(oldnames, paste0('a0s[', 1:K, ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:K))
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
