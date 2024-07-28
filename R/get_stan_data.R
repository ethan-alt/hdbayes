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


#' get Stan data for NAPP
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.napp = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1
) {
  data.checks(formula, family, data.list, offset.list)

  K        = length(data.list) - 1 # number of historical datasets
  if ( K == 0 ){
    stop("data.list should include at least one historical data set")
  }
  data     = data.list[[1]] # current data
  y        = data[, all.vars(formula)[1]]
  n        = length(y)
  X        = stats::model.matrix(formula, data)
  p        = ncol(X)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]
  ind.disp = ifelse(dist > 2, 1, 0)

  res          = stack.data(formula = formula,
                            data.list = data.list)
  N            = length(res$y) # total number of observations
  start.index  = res$start.index
  end.index    = res$end.index
  rm(res)

  ## Default offset for each dataset is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  offset.curr = offset[1:n] # offset for current data

  ## Fit enriched GLM to get hessian
  theta.means = matrix(NA, nrow = p + ind.disp, ncol = K)
  theta.covars = array(NA, c(K, p + ind.disp, p + ind.disp))
  for(k in 1:K) {
    histdata = data.list[[1+k]]
    offs0 = offset[ start.index[1+k]:end.index[1+k] ]
    histdata$offs0 = offs0
    fit.glm =
      suppressWarnings(
        stats::glm(formula = formula, family = family, data = histdata, offset = offs0)
      )
    fit.glm     = enrichwith::enrich(fit.glm)
    theta.mean  = fit.glm$coefficients
    theta.covar = fit.glm$expected_information_mle

    if ( !( family$family %in% c('binomial', 'poisson') ) ) {
      ## theta = (beta, log(dispersion) )
      theta.mean = c(theta.mean, log(fit.glm$dispersion_mle))

      ## jacobian adjustment to fisher information for dispersion --> log dispersion
      theta.covar[nrow(theta.covar), nrow(theta.covar)] =
        fit.glm$dispersion_mle^2 * theta.covar[nrow(theta.covar), nrow(theta.covar)]
    }
    theta.covar         = chol2inv(chol(theta.covar))
    theta.means[, k]    = theta.mean
    theta.covars[k, , ] = theta.covar
  }

  standat = list(
    'K'            = K,
    'n'            = n,
    'p'            = p,
    'ind_disp'     = ind.disp,
    'y'            = y,
    'X'            = X,
    'theta_hats'   = theta.means,
    'theta_covars' = theta.covars,
    'a0_shape1'    = a0.shape1,
    'a0_shape2'    = a0.shape2,
    'dist'         = dist,
    'link'         = link,
    'offs'         = offset.curr
  )
  return(standat)
}


#' get Stan data for commensurate prior (CP)
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.cp = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5
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
    'p_spike'         = p.spike,
    'mu_spike'        = spike.mean,
    'sigma_spike'     = spike.sd,
    'mu_slab'         = slab.mean,
    'sigma_slab'      = slab.sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standat)
}


#' get Stan data for LEAP
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.leap = function(
    formula,
    family,
    data.list,
    K                 = 2,
    prob.conc         = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    all.hist          = FALSE ## indicator for whether data.list consists of historical data sets only
) {
  data.checks.leap(formula, family, data.list, K, offset.list)

  ## get model information
  if (all.hist){
    hist.data.list   = data.list
    hist.offset.list = offset.list
  }else {
    data             = data.list[[1]]
    offset           = offset.list[[1]]
    y                = data[, all.vars(formula)[1]]
    n                = length(y)
    X                = model.matrix(formula, data)
    hist.data.list   = data.list[-1]
    hist.offset.list = offset.list[-1]

    ## Default offsets are matrices of 0s
    if ( is.null(offset) )
      offset = matrix(rep(0, n*K), ncol=K)
  }

  ## stack all historical data sets into one historical data set
  res.hist = stack.data(formula, hist.data.list)
  y0       = res.hist$y
  n0       = length(y0)
  X0       = res.hist$X
  p        = ncol(X0)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  ## Default prob.conc is a vector of 1s
  if ( !is.null(prob.conc) ){
    if ( !( is.vector(prob.conc) & (length(prob.conc) %in% c(1, K)) ) )
      stop("prob.conc must be a scalar or a vector of length ", K, " if prob.conc is not NULL")
  }
  prob.conc = to.vector(param = prob.conc, default.value = 1, len = K)

  ## Default offsets are matrices of 0s
  if ( is.null(hist.offset.list) ){
    offset0 = matrix(rep(0, n0*K), ncol=K)
  }else {
    offset0 = do.call(rbind, hist.offset.list)
  }

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) ){
    beta.mean = matrix(rep(0, p*K), ncol=K)
  }else {
    if ( length(beta.mean) == 1 ){
      beta.mean = matrix(rep(as.numeric(beta.mean), p*K), ncol=K)
    }else {
      if ( !( is.matrix(beta.mean) ) )
        stop("beta.mean must be a matrix")
      if ( nrow(beta.mean) != p )
        stop("beta.mean must have ", p, " row(s) if it is not NULL")
      if ( ncol(beta.mean) != K )
        stop("beta.mean must have ", K, " column(s) if it is not NULL")
      beta.mean = matrix(as.numeric(beta.mean), nrow = p, ncol = K)
    }
  }

  if ( is.null(beta.sd) ){
    beta.sd = matrix(rep(10, p*K), ncol=K)
  }else {
    if ( length(beta.sd) == 1 ){
      beta.sd = matrix(rep(as.numeric(beta.sd), p*K), ncol=K)
    }else {
      if ( !( is.matrix(beta.sd) ) )
        stop("beta.sd must be a matrix")
      if ( nrow(beta.sd) != p )
        stop("beta.sd must have ", p, " row(s) if it is not NULL")
      if ( ncol(beta.sd) != K )
        stop("beta.sd must have ", K, " column(s) if it is not NULL")
      beta.sd = matrix(as.numeric(beta.sd), nrow = p, ncol = K)
    }
  }

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

  if (all.hist){
    standat = list(
      'n0'          = n0,
      'p'           = p,
      'K'           = K,
      'y0'          = y0,
      'X0'          = X0,
      'mean_beta'   = beta.mean,
      'sd_beta'     = beta.sd,
      'disp_mean'   = disp.mean,
      'disp_sd'     = disp.sd,
      'conc'        = prob.conc,
      'gamma_lower' = 0,
      'gamma_upper' = 1,
      'dist'        = dist,
      'link'        = link,
      'offs0'       = offset0
    )
  }else {
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
      'sd_beta'     = beta.sd,
      'disp_mean'   = disp.mean,
      'disp_sd'     = disp.sd,
      'conc'        = prob.conc,
      'gamma_lower' = 0,
      'gamma_upper' = 1,
      'dist'        = dist,
      'link'        = link,
      'offs'        = offset,
      'offs0'       = offset0
    )
  }
  return(standat)
}

#' get Stan data for NPP
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.npp = function(
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
    a0.upper          = NULL
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

  ## check a0.lognc and lognc
  if ( length(a0.lognc) != nrow(lognc) )
    stop('the number of rows in lognc must be the same as the length of a0.lognc')
  if ( ncol(lognc) != (K - 1) )
    stop('the number of columns in lognc must be the same as the number of historical data sets')
  if ( any(is.na(a0.lognc) ) )
    stop('a0.lognc must not have missing values')
  if ( any(is.na(lognc)) )
    stop('lognc must not have missing values')
  if ( any(a0.lognc < 0) || any(a0.lognc > 1) )
    stop('each element of a0.lognc should be between 0 and 1')

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

  ## Default lower bound for each a0 is 0; default upper bound for each a0 is 1
  if ( !is.null(a0.lower) ){
    if ( !( is.vector(a0.lower) & (length(a0.lower) %in% c(1, K - 1)) ) )
      stop("a0.lower must be a scalar or a vector of length ", K - 1, " if a0.lower is not NULL")
  }
  a0.lower = to.vector(param = a0.lower, default.value = 0, len = K - 1)
  if ( !is.null(a0.upper) ){
    if ( !( is.vector(a0.upper) & (length(a0.upper) %in% c(1, K - 1)) ) )
      stop("a0.upper must be a scalar or a vector of length ", K - 1, " if a0.upper is not NULL")
  }
  a0.upper = to.vector(param = a0.upper, default.value = 1, len = K - 1)

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
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    's'               = length(a0.lognc),
    'a0_lognc'        = a0.lognc,
    'lognc'           = lognc,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standat)
}


#' get Stan data for normal/half-normal prior
#'
#' @include data_checks.R
#'
#' @noRd
get.stan.data.post = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL
) {
  if( length(data.list) > 1 ){
    data.list   = list(data.list[[1]])
  }
  if ( length(offset.list) > 1 ) {
    offset.list = list(offset.list[[1]])
  }
  if ( length(disp.mean) > 1 ) {
    disp.mean = disp.mean[1]
  }
  if ( length(disp.sd) > 1 ) {
    disp.sd = disp.sd[1]
  }

  standat_pp = get.stan.data.pp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    a0.vals     = 0,
    offset.list = offset.list,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd
  )

  standat = list(
    'n'               = standat_pp$N,
    'p'               = standat_pp$p,
    'y'               = standat_pp$y,
    'X'               = standat_pp$X,
    'mean_beta'       = standat_pp$mean_beta,
    'sd_beta'         = standat_pp$sd_beta,
    'disp_mean'       = standat_pp$disp_mean,
    'disp_sd'         = standat_pp$disp_sd,
    'dist'            = standat_pp$dist,
    'link'            = standat_pp$link,
    'offs'            = standat_pp$offs
  )
  return(standat)
}
