#' get Stan data for PP
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.pp = function(
    formula,
    data.list,
    dist              = "weibull",
    a0                = 0.5,
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE
) {
  data.checks.aft(formula, data.list, dist)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  t             = data[, time.name]
  y             = log(t)
  eventind      = as.integer( data[, eventind.name] )
  X             = stats::model.matrix(formula, data)
  p             = ncol(X)

  ## historical data
  if( length(data.list) == 1 ){
    t0        = 0
    y0        = log(t0)
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    t0            = histdata[, time.name]
    y0            = log(t0)
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)
  }

  ## Default prior on regression coefficients is N(0, 10^2)
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

  ## Default half-normal prior on scale parameter is N^{+}(0, 10^2)
  if ( !is.null(scale.mean) ){
    if ( !( is.vector(scale.mean) & (length(scale.mean) == 1) ) )
      stop("scale.mean must be a scalar if scale.mean is not NULL")
  }
  scale.mean = to.vector(param = scale.mean, default.value = 0, len = 1)
  if ( !is.null(scale.sd) ){
    if ( !( is.vector(scale.sd) & (length(scale.sd) == 1) ) )
      stop("scale.sd must be a scalar if scale.sd is not NULL")
  }
  scale.sd = to.vector(param = scale.sd, default.value = 10, len = 1)

  ## check a0 value
  a0 = as.numeric(a0)
  if ( any(a0 < 0 | a0 > 1 ) )
    stop("a0 must be a scalar between 0 and 1")

  standat = list(
    'dist'            = dist.to.integer(dist),
    'n'               = length(eventind),
    'n_obs'           = sum(eventind),
    'n_cen'           = sum(1 - eventind),
    'n0_obs'          = sum(eventind0),
    'n0_cen'          = sum(1 - eventind0),
    'p'               = p,
    'y_obs'           = y[which(eventind == 1)],
    'y_cen'           = y[which(eventind == 0)],
    'X_obs'           = X[which(eventind == 1), ],
    'X_cen'           = X[which(eventind == 0), ],
    'y0_obs'          = y0[which(eventind0 == 1)],
    'y0_cen'          = y0[which(eventind0 == 0)],
    'X0_obs'          = X0[which(eventind0 == 1), ],
    'X0_cen'          = X0[which(eventind0 == 0), ],
    'a0'              = a0,
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' get Stan data for normal/half-normal prior
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.post = function(
    formula,
    data.list,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE
) {
  if( length(data.list) > 1 ){
    data.list   = list(data.list[[1]])
  }

  standat_pp = get.aft.stan.data.pp(
    formula     = formula,
    data.list   = data.list,
    dist        = dist,
    a0          = 0,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    scale.mean  = scale.mean,
    scale.sd    = scale.sd,
    get.loglik  = get.loglik
  )

  standat = standat_pp[c("dist", 'n', 'n_obs', 'n_cen', 'p', 'y_obs', 'y_cen', 'X_obs', 'X_cen',
                         'beta_mean', 'beta_sd', 'scale_mean', 'scale_sd', 'get_loglik')]
  return(standat)
}


#' get Stan data for LEAP
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.leap = function(
    formula,
    data.list,
    dist              = "weibull",
    K                 = 2,
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    probs.conc        = NULL,
    get.loglik        = FALSE
) {
  data.checks.aft(formula, data.list, dist)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  t             = data[, time.name]
  y             = log(t)
  eventind      = as.integer( data[, eventind.name] )
  X             = stats::model.matrix(formula, data)
  p             = ncol(X)

  ## historical data
  if( length(data.list) == 1 ){
    t0        = 0
    y0        = log(t0)
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    t0            = histdata[, time.name]
    y0            = log(t0)
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)
  }

  ## Default probs.conc is a vector of 1s
  if ( !is.null(probs.conc) ){
    if ( !( is.vector(probs.conc) & (length(probs.conc) %in% c(1, K)) ) )
      stop("probs.conc must be a scalar or a vector of length ", K, " if probs.conc is not NULL")
  }
  probs.conc = to.vector(param = probs.conc, default.value = 1, len = K)

  ## Default prior on regression coefficients is N(0, 10^2)
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

  ## Default half-normal prior on scale parameter is N^{+}(0, 10^2)
  if ( !is.null(scale.mean) ){
    if ( !( is.vector(scale.mean) & (length(scale.mean) == 1) ) )
      stop("scale.mean must be a scalar if scale.mean is not NULL")
  }
  scale.mean = to.vector(param = scale.mean, default.value = 0, len = 1)
  if ( !is.null(scale.sd) ){
    if ( !( is.vector(scale.sd) & (length(scale.sd) == 1) ) )
      stop("scale.sd must be a scalar if scale.sd is not NULL")
  }
  scale.sd = to.vector(param = scale.sd, default.value = 10, len = 1)

  standat = list(
    'dist'            = dist.to.integer(dist),
    'n'               = length(eventind),
    'K'               = K,
    'n_obs'           = sum(eventind),
    'n_cen'           = sum(1 - eventind),
    'n0_obs'          = sum(eventind0),
    'n0_cen'          = sum(1 - eventind0),
    'p'               = p,
    'y_obs'           = y[which(eventind == 1)],
    'y_cen'           = y[which(eventind == 0)],
    'X_obs'           = X[which(eventind == 1), ],
    'X_cen'           = X[which(eventind == 0), ],
    'y0_obs'          = y0[which(eventind0 == 1)],
    'y0_cen'          = y0[which(eventind0 == 0)],
    'X0_obs'          = X0[which(eventind0 == 1), ],
    'X0_cen'          = X0[which(eventind0 == 0), ],
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'probs_conc'      = probs.conc,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' get Stan data for BHM
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.bhm = function(
    formula,
    data.list,
    dist              = "weibull",
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE
) {
  data.checks.aft(formula, data.list, dist)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  t             = data[, time.name]
  y             = log(t)
  eventind      = as.integer( data[, eventind.name] )
  X             = stats::model.matrix(formula, data)
  p             = ncol(X)

  ## historical data
  if( length(data.list) == 1 ){
    t0        = 0
    y0        = log(t0)
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    t0            = histdata[, time.name]
    y0            = log(t0)
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)
  }

  ## Default normal hyperprior on mean of regression coefficients is N(0, 10^2)
  if ( !is.null(meta.mean.mean) ){
    if ( !( is.vector(meta.mean.mean) & (length(meta.mean.mean) %in% c(1, p)) ) )
      stop("meta.mean.mean must be a scalar or a vector of length ", p, " if meta.mean.mean is not NULL")
  }
  meta.mean.mean = to.vector(param = meta.mean.mean, default.value = 0, len = p)
  if ( !is.null(meta.mean.sd) ){
    if ( !( is.vector(meta.mean.sd) & (length(meta.mean.sd) %in% c(1, p)) ) )
      stop("meta.mean.sd must be a scalar or a vector of length ", p, " if meta.mean.sd is not NULL")
  }
  meta.mean.sd = to.vector(param = meta.mean.sd, default.value = 10, len = p)

  ## Default half-normal hyperprior on sd of regression coefficients is N^{+}(0, 1)
  if ( !is.null(meta.sd.mean) ){
    if ( !( is.vector(meta.sd.mean) & (length(meta.sd.mean) %in% c(1, p)) ) )
      stop("meta.sd.mean must be a scalar or a vector of length ", p, " if meta.sd.mean is not NULL")
  }
  meta.sd.mean = to.vector(param = meta.sd.mean, default.value = 0, len = p)
  if ( !is.null(meta.sd.sd) ){
    if ( !( is.vector(meta.sd.sd) & (length(meta.sd.sd) %in% c(1, p)) ) )
      stop("meta.sd.sd must be a scalar or a vector of length ", p, " if meta.sd.sd is not NULL")
  }
  meta.sd.sd = to.vector(param = meta.sd.sd, default.value = 1, len = p)

  ## Default half-normal prior on scale parameter is N^{+}(0, 10^2)
  if ( !is.null(scale.mean) ){
    if ( !( is.vector(scale.mean) & (length(scale.mean) == 1) ) )
      stop("scale.mean must be a scalar if scale.mean is not NULL")
  }
  scale.mean = to.vector(param = scale.mean, default.value = 0, len = 1)
  if ( !is.null(scale.sd) ){
    if ( !( is.vector(scale.sd) & (length(scale.sd) == 1) ) )
      stop("scale.sd must be a scalar if scale.sd is not NULL")
  }
  scale.sd = to.vector(param = scale.sd, default.value = 10, len = 1)

  standat = list(
    'dist'            = dist.to.integer(dist),
    'n'               = length(eventind),
    'n_obs'           = sum(eventind),
    'n_cen'           = sum(1 - eventind),
    'n0_obs'          = sum(eventind0),
    'n0_cen'          = sum(1 - eventind0),
    'p'               = p,
    'y_obs'           = y[which(eventind == 1)],
    'y_cen'           = y[which(eventind == 0)],
    'X_obs'           = X[which(eventind == 1), ],
    'X_cen'           = X[which(eventind == 0), ],
    'y0_obs'          = y0[which(eventind0 == 1)],
    'y0_cen'          = y0[which(eventind0 == 0)],
    'X0_obs'          = X0[which(eventind0 == 1), ],
    'X0_cen'          = X0[which(eventind0 == 0), ],
    'meta_mean_mean'  = meta.mean.mean,
    'meta_mean_sd'    = meta.mean.sd,
    'meta_sd_mean'    = meta.sd.mean,
    'meta_sd_sd'      = meta.sd.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}
