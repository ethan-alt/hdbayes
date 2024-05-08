#' get Stan data for PP
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.pp = function(
    formula,
    family,
    data.list,
    a0.vals,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL
) {
  data.checks(formula, family, data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list)
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

  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
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

  ## check a0 values
  if ( !( is.vector(a0.vals) & (length(a0.vals) %in% c(1, K-1)) ) )
    stop("a0.vals must be a scalar or a vector of length ", K-1)
  a0.vals = to.vector(param = a0.vals, len = K-1)
  if ( any(a0.vals < 0 | a0.vals > 1 ) )
    stop("Each element of a0.vals must be a scalar between 0 and 1")
  a0.vals = c(1, a0.vals) # first element = 1 for current data

  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) == 1) ) )
      stop("disp.mean must be a scalar if disp.mean is not NULL")
  }
  disp.mean = to.vector(param = disp.mean, default.value = 0, len = 1)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) == 1) ) )
      stop("disp.sd must be a scalar if disp.sd is not NULL")
  }
  disp.sd = to.vector(param = disp.sd, default.value = 10, len = 1)

  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta.mean,
    'sd_beta'         = beta.sd,
    'a0_vals'         = a0.vals,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standat)
}


#' get Stan data for BHM
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.bhm = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL
) {
  data.checks(formula, family, data.list, offset.list)

  res          = stack.data(formula = formula,
                            data.list = data.list)
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

  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
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
    'meta_mean_mean'  = meta.mean.mean,
    'meta_mean_sd'    = meta.mean.sd,
    'meta_sd_mean'    = meta.sd.mean,
    'meta_sd_sd'      = meta.sd.sd,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standat)
}
