#' get Stan data for PP
#'
#' @include data_checks_pwe.R
#'
#' @noRd
get.pwe.stan.data.pp = function(
    formula,
    data.list,
    breaks,
    a0                = 0.5,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE
) {
  data.checks.pwe(formula, data.list, breaks)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  y1            = data[, time.name]
  eventind      = as.integer( data[, eventind.name] )
  X1            = stats::model.matrix(formula, data)

  ## make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X1) )
    X1 = X1[, -1, drop = F]

  p = ncol(X1)
  J = length(breaks) - 1  ## number of intervals

  ## historical data
  if( length(data.list) == 1 ){
    y0        = 0
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    y0            = histdata[, time.name]
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)

    if ( '(Intercept)' %in% colnames(X0) )
      X0 = X0[, -1, drop = F]
  }

  ## create index giving interval into which obs failed / was censored
  intindx = sapply(y1, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })
  intindx0 = sapply(y0, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })

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

  ## Default half-normal priors on baseline hazards are N^{+}(0, 10^2)
  if ( !is.null(base.hazard.mean) ){
    if ( !( is.vector(base.hazard.mean) & (length(base.hazard.mean) %in% c(1, J)) ) )
      stop("base.hazard.mean must be a scalar or a vector of length ", J, " if base.hazard.mean is not NULL")
  }
  base.hazard.mean = to.vector(param = base.hazard.mean, default.value = 0, len = J)
  if ( !is.null(base.hazard.sd) ){
    if ( !( is.vector(base.hazard.sd) & (length(base.hazard.sd) %in% c(1, J)) ) )
      stop("base.hazard.sd must be a scalar or a vector of length ", J, " if base.hazard.sd is not NULL")
  }
  base.hazard.sd = to.vector(param = base.hazard.sd, default.value = 10, len = J)

  ## check a0 value
  a0 = as.numeric(a0)
  if ( any(a0 < 0 | a0 > 1 ) )
    stop("a0 must be a scalar between 0 and 1")

  standat = list(
    'n1'                = nrow(X1)
    , 'n0'              = nrow(X0)
    , 'J'               = J
    , 'p'               = p
    , 'y1'              = y1
    , 'y0'              = y0
    , 'X1'              = X1
    , 'X0'              = X0
    , 'intindx'         = intindx
    , 'intindx0'        = intindx0
    , 'death_ind'       = eventind
    , 'death_ind0'      = eventind0
    , 'breaks'          = breaks
    , 'a0'              = a0
    , 'beta_mean'       = beta.mean
    , 'beta_sd'         = beta.sd
    , 'hazard_mean'     = base.hazard.mean
    , 'hazard_sd'       = base.hazard.sd
    , 'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' get Stan data for normal/half-normal prior
#'
#' @include data_checks_pwe.R
#'
#' @noRd
get.pwe.stan.data.post = function(
    formula,
    data.list,
    breaks,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE
) {
  if( length(data.list) > 1 ){
    data.list   = list(data.list[[1]])
  }

  standat_pp = get.pwe.stan.data.pp(
    formula          = formula,
    data.list        = data.list,
    breaks           = breaks,
    a0               = 0,
    beta.mean        = beta.mean,
    beta.sd          = beta.sd,
    base.hazard.mean = base.hazard.mean,
    base.hazard.sd   = base.hazard.sd,
    get.loglik       = get.loglik
  )
  standat = standat_pp[c('n1', 'J', 'p', 'y1', 'X1', 'intindx', 'death_ind', 'breaks',
                         'beta_mean', 'beta_sd', 'hazard_mean', 'hazard_sd', 'get_loglik')]
  return(standat)
}


#' get Stan data for LEAP
#'
#' @include data_checks_pwe.R
#'
#' @noRd
get.pwe.stan.data.leap = function(
    formula,
    data.list,
    breaks,
    K                 = 2,
    prob.conc         = NULL,
    gamma.lower       = 0,
    gamma.upper       = 1,
    beta.mean         = NULL,
    beta.sd           = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE
) {
  data.checks.pwe(formula, data.list, breaks)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  y1            = data[, time.name]
  eventind      = as.integer( data[, eventind.name] )
  X1            = stats::model.matrix(formula, data)

  ## make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X1) )
    X1 = X1[, -1, drop = F]

  p = ncol(X1)
  J = length(breaks) - 1  ## number of intervals

  ## historical data
  if( length(data.list) == 1 ){
    y0        = 0
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    y0            = histdata[, time.name]
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)

    if ( '(Intercept)' %in% colnames(X0) )
      X0 = X0[, -1, drop = F]
  }

  ## create index giving interval into which obs failed / was censored
  intindx = sapply(y1, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })
  intindx0 = sapply(y0, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })

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

  ## Default half-normal priors on baseline hazards are N^{+}(0, 10^2)
  if ( !is.null(base.hazard.mean) ){
    if ( !( is.vector(base.hazard.mean) & (length(base.hazard.mean) %in% c(1, J)) ) )
      stop("base.hazard.mean must be a scalar or a vector of length ", J, " if base.hazard.mean is not NULL")
  }
  base.hazard.mean = to.vector(param = base.hazard.mean, default.value = 0, len = J)
  if ( !is.null(base.hazard.sd) ){
    if ( !( is.vector(base.hazard.sd) & (length(base.hazard.sd) %in% c(1, J)) ) )
      stop("base.hazard.sd must be a scalar or a vector of length ", J, " if base.hazard.sd is not NULL")
  }
  base.hazard.sd = to.vector(param = base.hazard.sd, default.value = 10, len = J)

  ## gamma.upper should be smaller than or equal to 1
  if ( gamma.upper > 1 )
    gamma.upper = 1
  ## gamma.lower should be larger than or equal to 0
  if ( gamma.lower < 0 )
    gamma.lower = 0

  standat = list(
    'n1'                = nrow(X1)
    , 'n0'              = nrow(X0)
    , 'J'               = J
    , 'p'               = p
    , 'K'               = K
    , 'y1'              = y1
    , 'y0'              = y0
    , 'X1'              = X1
    , 'X0'              = X0
    , 'intindx'         = intindx
    , 'intindx0'        = intindx0
    , 'death_ind'       = eventind
    , 'death_ind0'      = eventind0
    , 'breaks'          = breaks
    , 'conc'            = prob.conc
    , 'gamma_lower'     = gamma.lower
    , 'gamma_upper'     = gamma.upper
    , 'beta_mean'       = beta.mean
    , 'beta_sd'         = beta.sd
    , 'hazard_mean'     = base.hazard.mean
    , 'hazard_sd'       = base.hazard.sd
    , 'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' get Stan data for BHM
#'
#' @include data_checks_pwe.R
#'
#' @noRd
get.pwe.stan.data.bhm = function(
    formula,
    data.list,
    breaks,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE
) {
  data.checks.pwe(formula, data.list, breaks)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  y1            = data[, time.name]
  eventind      = as.integer( data[, eventind.name] )
  X1            = stats::model.matrix(formula, data)

  ## make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X1) )
    X1 = X1[, -1, drop = F]

  p = ncol(X1)
  J = length(breaks) - 1  ## number of intervals

  ## historical data
  if( length(data.list) == 1 ){
    y0        = 0
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    y0            = histdata[, time.name]
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)

    if ( '(Intercept)' %in% colnames(X0) )
      X0 = X0[, -1, drop = F]
  }

  ## create index giving interval into which obs failed / was censored
  intindx = sapply(y1, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })
  intindx0 = sapply(y0, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })

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

  ## Default half-normal priors on baseline hazards are N^{+}(0, 10^2)
  if ( !is.null(base.hazard.mean) ){
    if ( !( is.vector(base.hazard.mean) & (length(base.hazard.mean) %in% c(1, J)) ) )
      stop("base.hazard.mean must be a scalar or a vector of length ", J, " if base.hazard.mean is not NULL")
  }
  base.hazard.mean = to.vector(param = base.hazard.mean, default.value = 0, len = J)
  if ( !is.null(base.hazard.sd) ){
    if ( !( is.vector(base.hazard.sd) & (length(base.hazard.sd) %in% c(1, J)) ) )
      stop("base.hazard.sd must be a scalar or a vector of length ", J, " if base.hazard.sd is not NULL")
  }
  base.hazard.sd = to.vector(param = base.hazard.sd, default.value = 10, len = J)

  standat = list(
    'n1'                = nrow(X1)
    , 'n0'              = nrow(X0)
    , 'J'               = J
    , 'p'               = p
    , 'y1'              = y1
    , 'y0'              = y0
    , 'X1'              = X1
    , 'X0'              = X0
    , 'intindx'         = intindx
    , 'intindx0'        = intindx0
    , 'death_ind'       = eventind
    , 'death_ind0'      = eventind0
    , 'breaks'          = breaks
    , 'meta_mean_mean'  = meta.mean.mean
    , 'meta_mean_sd'    = meta.mean.sd
    , 'meta_sd_mean'    = meta.sd.mean
    , 'meta_sd_sd'      = meta.sd.sd
    , 'hazard_mean'     = base.hazard.mean
    , 'hazard_sd'       = base.hazard.sd
    , 'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}


#' get Stan data for commensurate prior (CP)
#'
#' @include data_checks_pwe.R
#'
#' @noRd
get.pwe.stan.data.cp = function(
    formula,
    data.list,
    breaks,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
    base.hazard.mean  = NULL,
    base.hazard.sd    = NULL,
    get.loglik        = FALSE
) {
  data.checks.pwe(formula, data.list, breaks)

  ## current data
  data          = data.list[[1]]
  ## extract names and variables for response, censoring, etc.
  time.name     = all.vars(formula)[1]
  eventind.name = all.vars(formula)[2]
  y1            = data[, time.name]
  eventind      = as.integer( data[, eventind.name] )
  X1            = stats::model.matrix(formula, data)

  ## make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X1) )
    X1 = X1[, -1, drop = F]

  p = ncol(X1)
  J = length(breaks) - 1  ## number of intervals

  ## historical data
  if( length(data.list) == 1 ){
    y0        = 0
    eventind0 = 0
    X0        = matrix(0, 1, p)
  } else{
    histdata      = do.call(rbind, data.list[-1])
    y0            = histdata[, time.name]
    eventind0     = as.integer( histdata[, eventind.name] )
    X0            = stats::model.matrix(formula, histdata)

    if ( '(Intercept)' %in% colnames(X0) )
      X0 = X0[, -1, drop = F]
  }

  ## create index giving interval into which obs failed / was censored
  intindx = sapply(y1, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })
  intindx0 = sapply(y0, function(t){
    if( t == 0 ){
      return(1)
    }else{
      return(findInterval(t, breaks, left.open = TRUE))
    }
  })

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

  ## Default half-normal priors on baseline hazards are N^{+}(0, 10^2)
  if ( !is.null(base.hazard.mean) ){
    if ( !( is.vector(base.hazard.mean) & (length(base.hazard.mean) %in% c(1, J)) ) )
      stop("base.hazard.mean must be a scalar or a vector of length ", J, " if base.hazard.mean is not NULL")
  }
  base.hazard.mean = to.vector(param = base.hazard.mean, default.value = 0, len = J)
  if ( !is.null(base.hazard.sd) ){
    if ( !( is.vector(base.hazard.sd) & (length(base.hazard.sd) %in% c(1, J)) ) )
      stop("base.hazard.sd must be a scalar or a vector of length ", J, " if base.hazard.sd is not NULL")
  }
  base.hazard.sd = to.vector(param = base.hazard.sd, default.value = 10, len = J)

  standat = list(
    'n1'                = nrow(X1)
    , 'n0'              = nrow(X0)
    , 'J'               = J
    , 'p'               = p
    , 'y1'              = y1
    , 'y0'              = y0
    , 'X1'              = X1
    , 'X0'              = X0
    , 'intindx'         = intindx
    , 'intindx0'        = intindx0
    , 'death_ind'       = eventind
    , 'death_ind0'      = eventind0
    , 'breaks'          = breaks
    , 'beta0_mean'      = beta0.mean
    , 'beta0_sd'        = beta0.sd
    , 'p_spike'         = p.spike
    , 'mu_spike'        = spike.mean
    , 'sigma_spike'     = spike.sd
    , 'mu_slab'         = slab.mean
    , 'sigma_slab'      = slab.sd
    , 'hazard_mean'     = base.hazard.mean
    , 'hazard_sd'       = base.hazard.sd
    , 'get_loglik'      = as.integer(get.loglik)
  )
  return(standat)
}
