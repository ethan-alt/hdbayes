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
    prob.conc         = NULL,
    gamma.lower       = 0,
    gamma.upper       = 1,
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

  ## Default prob.conc is a vector of 1s
  if ( !is.null(prob.conc) ){
    if ( !( is.vector(prob.conc) & (length(prob.conc) %in% c(1, K)) ) )
      stop("prob.conc must be a scalar or a vector of length ", K, " if prob.conc is not NULL")
  }
  prob.conc = to.vector(param = prob.conc, default.value = 1, len = K)

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

  ## gamma.upper should be smaller than or equal to 1
  if ( gamma.upper > 1 )
    gamma.upper = 1
  ## gamma.lower should be larger than or equal to 0
  if ( gamma.lower < 0 )
    gamma.lower = 0

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
    'prob_conc'       = prob.conc,
    'gamma_lower'     = gamma.lower,
    'gamma_upper'     = gamma.upper,
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

#' get Stan data for commensurate prior (CP)
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.cp = function(
    formula,
    data.list,
    dist              = "weibull",
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
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

  ## check p.spike
  p.spike = as.numeric(p.spike)
  if ( p.spike < 0 | p.spike > 1 )
    stop("p.spike must be a scalar between 0 and 1")

  ## Default half-normal prior (spike component) on commensurability parameter is N^{+}(200, 0.1)
  spike.mean = as.numeric(spike.mean)
  spike.sd   = as.numeric(spike.sd)

  ## Default half-normal prior (slab component) on commensurability parameter is N^{+}(0, 5)
  slab.mean = as.numeric(slab.mean)
  slab.sd   = as.numeric(slab.sd)

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
    'beta0_mean'      = beta0.mean,
    'beta0_sd'        = beta0.sd,
    'p_spike'         = p.spike,
    'mu_spike'        = spike.mean,
    'sigma_spike'     = spike.sd,
    'mu_slab'         = slab.mean,
    'sigma_slab'      = slab.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}

#' get Stan data for NPP
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.npp = function(
    formula,
    data.list,
    a0.lognc,
    lognc,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = 0,
    a0.upper          = 1,
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

  ## check a0.lognc and lognc
  if ( length(a0.lognc) != length(lognc) )
    stop('the length of lognc must be the same as that of a0.lognc')
  if ( any(is.na(a0.lognc) ) )
    stop('a0.lognc must not have missing values')
  if ( any(is.na(lognc)) )
    stop('lognc must not have missing values')
  if ( any(a0.lognc < 0) || any(a0.lognc > 1) )
    stop('each element of a0.lognc should be between 0 and 1')

  ## a0.upper should be smaller than or equal to 1
  if ( a0.upper > 1 )
    a0.upper = 1
  ## a0.lower should be larger than or equal to 0
  if ( a0.lower < 0 )
    a0.lower = 0

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
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    's'               = length(a0.lognc),
    'a0_lognc'        = a0.lognc,
    'lognc'           = lognc,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}

#' get Stan data for stratified PP
#'
#' @include data_checks_aft.R
#'
#' @noRd
get.aft.stan.data.stratified.pp = function(
    formula,
    data.list,
    strata.list,
    a0.strata,
    dist              = "weibull",
    beta.mean         = NULL,
    beta.sd           = NULL,
    scale.mean        = NULL,
    scale.sd          = NULL,
    get.loglik        = FALSE
) {
  data.checks.aft(formula, data.list, dist,
                  strata.list = strata.list, is.stratified.pp = TRUE)

  ## current data
  data          = data.list[[1]]
  stratum.curr  = strata.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  ## re-order data by event indicator and stratum
  data          = data[order(data[, eventind.name], stratum.curr), ]
  t             = data[, time.name]
  y             = log(t)
  eventind      = as.integer( data[, eventind.name] )
  X             = stats::model.matrix(formula, data)
  p             = ncol(X)

  ## historical data
  histdata      = do.call(rbind, data.list[-1])
  stratum.hist  = do.call(c, strata.list[-1])
  ## re-order histdata by event indicator and stratum
  histdata      = histdata[order(histdata[, eventind.name], stratum.hist), ]
  t0            = histdata[, time.name]
  y0            = log(t0)
  eventind0     = as.integer( histdata[, eventind.name] )
  X0            = stats::model.matrix(formula, histdata)

  ## get the number of strata
  K = as.integer( max( stratum.curr, stratum.hist ) )

  ## get strata assignment for current and historical data
  stratumID.obs  = stratum.curr[which(eventind == 1)]
  stratumID.cen  = stratum.curr[which(eventind == 0)]
  stratumID0.obs = stratum.hist[which(eventind0 == 1)]
  stratumID0.cen = stratum.hist[which(eventind0 == 0)]

  ## check a0.strata values
  if ( !( is.vector(a0.strata) & (length(a0.strata) %in% c(1, K)) ) )
    stop("a0.strata must be a scalar or a vector of length ", K)
  a0.strata = to.vector(param = a0.strata, len = K)
  if ( any(a0.strata < 0 | a0.strata > 1 ) )
    stop("Each element of a0.strata must be a scalar between 0 and 1")

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
    'K'               = K,
    'stratumID_obs'   = stratumID.obs,
    'stratumID_cen'   = stratumID.cen,
    'stratumID0_obs'  = stratumID0.obs,
    'stratumID0_cen'  = stratumID0.cen,
    'a0s'             = a0.strata,
    'beta_mean'       = beta.mean,
    'beta_sd'         = beta.sd,
    'scale_mean'      = scale.mean,
    'scale_sd'        = scale.sd,
    'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}
