#' Get Stan data for GLM models with missing covariates (Gaussian copula)
#'
#' @include data_checks_mc.R
#' @include missing_data_utils.R
#'
#' @name get_stan_data_mc
#' @keywords internal
NULL


#' Get Stan data for BHM with missing covariates
#'
#' @param formula         a two-sided formula giving the relationship between the response variable and covariates.
#' @param family          an object of class `family`.
#' @param data.list       a list of `data.frame`s.
#' @param copula.corr     optional user-provided correlation matrix for the copula.
#' @param offset.list     a list of vectors giving the offsets for each data.
#' @param meta.mean.mean  hyperprior mean for regression coefficient means.
#' @param meta.mean.sd    hyperprior sd for regression coefficient means.
#' @param meta.sd.mean    hyperprior mean for regression coefficient sds.
#' @param meta.sd.sd      hyperprior sd for regression coefficient sds.
#' @param disp.mean       prior mean for dispersion parameters.
#' @param disp.sd         prior sd for dispersion parameters.
#' @param get.loglik      whether to generate log-likelihood matrix.
#' @noRd
get.stan.data.bhm.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    get.loglik        = FALSE
) {
  # Validate data (allowing missing covariates)
  data.checks.mc(formula, family, data.list, offset.list)

  # Prepare copula data
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    stop("No missing covariates detected. Use glm.bhm() instead of glm.bhm.mc()")
  }

  # Stack data with missing value handling
  stacked = stack.data.mc(formula, data.list, copula_info)

  y            = stacked$y
  X            = stacked$X
  start.index  = stacked$start_index
  end.index    = stacked$end_index
  p            = stacked$p
  N            = stacked$N
  K            = stacked$K
  miss_idx     = stacked$miss_idx
  n_miss       = stacked$n_miss
  q            = stacked$q

  fam.indx     = get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  # Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  } else {
    offset = unlist(offset.list)
  }

  # Default normal hyperprior on mean of regression coefficients is N(0, 10^2)
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

  # Default half-normal hyperprior on sd of regression coefficients is N^{+}(0, 1)
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

  # Default half-normal hyperprior on dispersion parameters (if exist) is N^{+}(0, 10^2)
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

  # Prepare marginal parameters for Stan
  marginal_stan = prepare.marginal.stan(
    copula_info$marginal_types,
    copula_info$marginal_params,
    copula_info$all_missing_covs
  )

  # Build Stan data list
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
    'offs'            = offset,
    'get_loglik'      = as.integer(get.loglik),
    # Missing data parameters
    'n_miss'          = n_miss,
    'q'               = q,
    'miss_idx'        = miss_idx,
    'copula_corr'     = copula_info$copula_corr,
    'marginal_type'   = marginal_stan$marginal_type,
    'marginal_mean'   = marginal_stan$marginal_mean,
    'marginal_sd'     = marginal_stan$marginal_sd,
    'n_levels'        = marginal_stan$n_levels,
    'cutpoints'       = marginal_stan$cutpoints,
    'max_levels'      = marginal_stan$max_levels
  )

  return(standat)
}


#' Get Stan data for PP (Power Prior) with missing covariates
#'
#' @noRd
get.stan.data.pp.mc = function(
    formula,
    family,
    data.list,
    a0.vals,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    get.loglik        = FALSE
) {
  # Validate data (allowing missing covariates)
  data.checks.mc(formula, family, data.list, offset.list)

  # Prepare copula data
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    stop("No missing covariates detected. Use glm.pp() instead of glm.pp.mc()")
  }

  # Stack data with missing value handling
  stacked = stack.data.mc(formula, data.list, copula_info)

  y            = stacked$y
  X            = stacked$X
  start.index  = stacked$start_index
  end.index    = stacked$end_index
  p            = stacked$p
  N            = stacked$N
  K            = stacked$K
  miss_idx     = stacked$miss_idx
  n_miss       = stacked$n_miss
  q            = stacked$q

  fam.indx     = get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  # Default offset
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  } else {
    offset = unlist(offset.list)
  }

  # Default prior on regression coefficients is N(0, 10^2)
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

  # Check a0 values
  if ( !( is.vector(a0.vals) & (length(a0.vals) %in% c(1, K-1)) ) )
    stop("a0.vals must be a scalar or a vector of length ", K-1)
  a0.vals = to.vector(param = a0.vals, len = K-1)
  if ( any(a0.vals < 0 | a0.vals > 1 ) )
    stop("Each element of a0.vals must be a scalar between 0 and 1")
  a0.vals = c(1, a0.vals) # first element = 1 for current data

  # Default half-normal prior on dispersion parameters
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

  # Prepare marginal parameters for Stan
  marginal_stan = prepare.marginal.stan(
    copula_info$marginal_types,
    copula_info$marginal_params,
    copula_info$all_missing_covs
  )

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
    'offs'            = offset,
    'get_loglik'      = as.integer(get.loglik),
    # Missing data parameters
    'n_miss'          = n_miss,
    'q'               = q,
    'miss_idx'        = miss_idx,
    'copula_corr'     = copula_info$copula_corr,
    'marginal_type'   = marginal_stan$marginal_type,
    'marginal_mean'   = marginal_stan$marginal_mean,
    'marginal_sd'     = marginal_stan$marginal_sd,
    'n_levels'        = marginal_stan$n_levels,
    'cutpoints'       = marginal_stan$cutpoints,
    'max_levels'      = marginal_stan$max_levels
  )

  return(standat)
}


#' Get Stan data for normal/half-normal prior (reference prior) with missing covariates
#'
#' @noRd
get.stan.data.post.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    get.loglik        = FALSE
) {
  # Use only current data (first element)
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

  standat_pp = get.stan.data.pp.mc(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    a0.vals     = 0,
    copula.corr = copula.corr,
    offset.list = offset.list,
    beta.mean   = beta.mean,
    beta.sd     = beta.sd,
    disp.mean   = disp.mean,
    disp.sd     = disp.sd,
    get.loglik  = get.loglik
  )

  # Restructure for posterior-only model
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
    'offs'            = standat_pp$offs,
    'get_loglik'      = standat_pp$get_loglik,
    # Missing data parameters
    'n_miss'          = standat_pp$n_miss,
    'q'               = standat_pp$q,
    'miss_idx'        = standat_pp$miss_idx,
    'copula_corr'     = standat_pp$copula_corr,
    'marginal_type'   = standat_pp$marginal_type,
    'marginal_mean'   = standat_pp$marginal_mean,
    'marginal_sd'     = standat_pp$marginal_sd,
    'n_levels'        = standat_pp$n_levels,
    'cutpoints'       = standat_pp$cutpoints,
    'max_levels'      = standat_pp$max_levels
  )

  return(standat)
}


#' Get Stan data for commensurate prior (CP) with missing covariates
#'
#' @noRd
get.stan.data.cp.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    p.spike           = 0.1,
    spike.mean        = 200,
    spike.sd          = 0.1,
    slab.mean         = 0,
    slab.sd           = 5,
    get.loglik        = FALSE
) {
  # Validate data (allowing missing covariates)
  data.checks.mc(formula, family, data.list, offset.list)

  # Prepare copula data
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    stop("No missing covariates detected. Use glm.commensurate() instead of glm.commensurate.mc()")
  }

  # Stack data with missing value handling
  stacked = stack.data.mc(formula, data.list, copula_info)

  y            = stacked$y
  X            = stacked$X
  start.index  = stacked$start_index
  end.index    = stacked$end_index
  p            = stacked$p
  N            = stacked$N
  K            = stacked$K
  miss_idx     = stacked$miss_idx
  n_miss       = stacked$n_miss
  q            = stacked$q

  fam.indx     = get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  # Default offset
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  } else {
    offset = unlist(offset.list)
  }

  # Default prior on regression coefficients is N(0, 10^2)
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

  # Default half-normal hyperprior on dispersion parameters
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

  # Check p.spike
  p.spike = as.numeric(p.spike)
  if ( p.spike < 0 | p.spike > 1 )
    stop("p.spike must be a scalar between 0 and 1")

  spike.mean = as.numeric(spike.mean)
  spike.sd   = as.numeric(spike.sd)
  slab.mean = as.numeric(slab.mean)
  slab.sd   = as.numeric(slab.sd)

  # Prepare marginal parameters for Stan
  marginal_stan = prepare.marginal.stan(
    copula_info$marginal_types,
    copula_info$marginal_params,
    copula_info$all_missing_covs
  )

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
    'offs'            = offset,
    'get_loglik'      = as.integer(get.loglik),
    # Missing data parameters
    'n_miss'          = n_miss,
    'q'               = q,
    'miss_idx'        = miss_idx,
    'copula_corr'     = copula_info$copula_corr,
    'marginal_type'   = marginal_stan$marginal_type,
    'marginal_mean'   = marginal_stan$marginal_mean,
    'marginal_sd'     = marginal_stan$marginal_sd,
    'n_levels'        = marginal_stan$n_levels,
    'cutpoints'       = marginal_stan$cutpoints,
    'max_levels'      = marginal_stan$max_levels
  )

  return(standat)
}


#' Get Stan data for LEAP (Latent Exchangeability Prior) with missing covariates
#'
#' @noRd
get.stan.data.leap.mc = function(
    formula,
    family,
    data.list,
    K                 = 2,
    prob.conc         = NULL,
    copula.corr       = NULL,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    gamma.lower       = 0,
    gamma.upper       = 1,
    get.loglik        = FALSE
) {
  # Validate data (allowing missing covariates)
  data.checks.mc.leap(formula, family, data.list, K, offset.list)

  # Prepare copula data from all datasets
  copula_info = prepare.copula.data(formula, data.list, copula.corr)

  if (is.null(copula_info)) {
    stop("No missing covariates detected. Use glm.leap() instead of glm.leap.mc()")
  }

  # Current data
  data         = data.list[[1]]
  y            = data[, all.vars(formula)[1]]
  n            = length(y)

  # Historical data - stack all historical datasets
  hist.data.list   = data.list[-1]
  hist.offset.list = if (!is.null(offset.list)) offset.list[-1] else NULL

  # Build design matrices with missing value handling
  # For current data - fill in missing with means for model.matrix
  data_filled = data
  for (cov in copula_info$all_missing_covs) {
    if (any(is.na(data[[cov]]))) {
      data_filled[[cov]][is.na(data[[cov]])] = mean(data[[cov]], na.rm = TRUE)
    }
  }
  X = stats::model.matrix(formula, data_filled)
  p = ncol(X)

  # Track missing indices for current data
  miss_pattern_curr = copula_info$miss_patterns[[1]]
  miss_idx = build.miss.idx.for.design.matrix(miss_pattern_curr, X, copula_info$all_missing_covs)
  n_miss = nrow(miss_idx)

  # For historical data - stack and fill
  hist_data_stacked = do.call(rbind, hist.data.list)
  for (cov in copula_info$all_missing_covs) {
    if (any(is.na(hist_data_stacked[[cov]]))) {
      hist_data_stacked[[cov]][is.na(hist_data_stacked[[cov]])] = mean(hist_data_stacked[[cov]], na.rm = TRUE)
    }
  }
  X0 = stats::model.matrix(formula, hist_data_stacked)
  y0 = hist_data_stacked[, all.vars(formula)[1]]
  n0 = length(y0)

  # Track missing indices for historical data
  miss_idx0_list = list()
  row_offset = 0
  for (k in seq_along(hist.data.list)) {
    miss_pattern_k = copula_info$miss_patterns[[k + 1]]  # +1 because first is current
    if (miss_pattern_k$n_miss > 0) {
      miss_idx_k = build.miss.idx.for.design.matrix(miss_pattern_k, X0, copula_info$all_missing_covs, row_offset)
      if (nrow(miss_idx_k) > 0) {
        miss_idx0_list[[length(miss_idx0_list) + 1]] = miss_idx_k
      }
    }
    row_offset = row_offset + nrow(hist.data.list[[k]])
  }

  if (length(miss_idx0_list) > 0) {
    miss_idx0 = do.call(rbind, miss_idx0_list)
  } else {
    miss_idx0 = matrix(integer(0), ncol = 3)
  }
  n_miss0 = nrow(miss_idx0)

  # Ensure miss_idx has at least one row for Stan (will be ignored if n_miss=0)
  if (n_miss == 0) {
    miss_idx = matrix(c(1L, 1L, 1L), nrow = 1)
  }
  if (n_miss0 == 0) {
    miss_idx0 = matrix(c(1L, 1L, 1L), nrow = 1)
  }

  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  q = length(copula_info$all_missing_covs)

  # Default prob.conc is a vector of 1s
  if (!is.null(prob.conc)) {
    if (!(is.vector(prob.conc) & (length(prob.conc) %in% c(1, K))))
      stop("prob.conc must be a scalar or a vector of length ", K, " if prob.conc is not NULL")
  }
  prob.conc = to.vector(param = prob.conc, default.value = 1, len = K)

  # Default offsets are matrices of 0s
  offset = offset.list[[1]]
  if (is.null(offset))
    offset = matrix(rep(0, n * K), ncol = K)
  if (is.null(hist.offset.list)) {
    offset0 = matrix(rep(0, n0 * K), ncol = K)
  } else {
    offset0 = do.call(rbind, hist.offset.list)
  }

  # Default prior on regression coefficients is N(0, 10^2)
  if (is.null(beta.mean)) {
    beta.mean = matrix(rep(0, p * K), ncol = K)
  } else {
    if (length(beta.mean) == 1) {
      beta.mean = matrix(rep(as.numeric(beta.mean), p * K), ncol = K)
    } else {
      if (!is.matrix(beta.mean))
        stop("beta.mean must be a matrix")
      if (nrow(beta.mean) != p)
        stop("beta.mean must have ", p, " row(s) if it is not NULL")
      if (ncol(beta.mean) != K)
        stop("beta.mean must have ", K, " column(s) if it is not NULL")
      beta.mean = matrix(as.numeric(beta.mean), nrow = p, ncol = K)
    }
  }

  if (is.null(beta.sd)) {
    beta.sd = matrix(rep(10, p * K), ncol = K)
  } else {
    if (length(beta.sd) == 1) {
      beta.sd = matrix(rep(as.numeric(beta.sd), p * K), ncol = K)
    } else {
      if (!is.matrix(beta.sd))
        stop("beta.sd must be a matrix")
      if (nrow(beta.sd) != p)
        stop("beta.sd must have ", p, " row(s) if it is not NULL")
      if (ncol(beta.sd) != K)
        stop("beta.sd must have ", K, " column(s) if it is not NULL")
      beta.sd = matrix(as.numeric(beta.sd), nrow = p, ncol = K)
    }
  }

  # Default half-normal prior on dispersion parameters
  if (!is.null(disp.mean)) {
    if (!(is.vector(disp.mean) & (length(disp.mean) %in% c(1, K))))
      stop("disp.mean must be a scalar or a vector of length ", K, " if disp.mean is not NULL")
  }
  disp.mean = to.vector(param = disp.mean, default.value = 0, len = K)
  if (!is.null(disp.sd)) {
    if (!(is.vector(disp.sd) & (length(disp.sd) %in% c(1, K))))
      stop("disp.sd must be a scalar or a vector of length ", K, " if disp.sd is not NULL")
  }
  disp.sd = to.vector(param = disp.sd, default.value = 10, len = K)

  # Bounds on gamma
  if (gamma.upper > 1) gamma.upper = 1
  if (gamma.lower < 0) gamma.lower = 0

  # Prepare marginal parameters for Stan
  marginal_stan = prepare.marginal.stan(
    copula_info$marginal_types,
    copula_info$marginal_params,
    copula_info$all_missing_covs
  )

  # Ensure q >= 1 for Stan array dimensions
  q_stan = if (q == 0) 1 else q

  standat = list(
    'n'             = n,
    'n0'            = n0,
    'p'             = p,
    'K'             = K,
    'y'             = y,
    'X'             = X,
    'y0'            = y0,
    'X0'            = X0,
    'mean_beta'     = beta.mean,
    'sd_beta'       = beta.sd,
    'disp_mean'     = disp.mean,
    'disp_sd'       = disp.sd,
    'conc'          = prob.conc,
    'gamma_lower'   = gamma.lower,
    'gamma_upper'   = gamma.upper,
    'dist'          = dist,
    'link'          = link,
    'offs'          = offset,
    'offs0'         = offset0,
    'get_loglik'    = as.integer(get.loglik),
    # Missing data parameters
    'n_miss'        = n_miss,
    'n_miss0'       = n_miss0,
    'q'             = q_stan,
    'miss_idx'      = miss_idx,
    'miss_idx0'     = miss_idx0,
    'copula_corr'   = copula_info$copula_corr,
    'marginal_type' = marginal_stan$marginal_type,
    'marginal_mean' = marginal_stan$marginal_mean,
    'marginal_sd'   = marginal_stan$marginal_sd,
    'n_levels'      = marginal_stan$n_levels,
    'cutpoints'     = marginal_stan$cutpoints,
    'max_levels'    = marginal_stan$max_levels
  )

  return(standat)
}


#' Helper function to build missing index matrix for a design matrix
#'
#' @param miss_pattern output from detect.missing.pattern
#' @param X design matrix
#' @param all_missing_covs character vector of all missing covariate names
#' @param row_offset offset to add to row indices (for stacked data)
#' @return matrix with columns (row, dm_col, cov_idx)
#' @noRd
build.miss.idx.for.design.matrix = function(miss_pattern, X, all_missing_covs, row_offset = 0) {
  if (miss_pattern$n_miss == 0) {
    return(matrix(integer(0), ncol = 3))
  }

  dm_colnames = colnames(X)
  miss_idx_list = list()

  for (m in seq_len(nrow(miss_pattern$miss_idx))) {
    row_in_data = miss_pattern$miss_idx[m, "row"]
    cov_idx = miss_pattern$miss_idx[m, "col"]
    cov_name = miss_pattern$missing_covs[cov_idx]

    # Find which design matrix columns correspond to this covariate
    dm_cols = which(dm_colnames == cov_name)
    if (length(dm_cols) == 0) {
      dm_cols = grep(cov_name, dm_colnames, fixed = TRUE)
    }
    if (length(dm_cols) == 0) {
      dm_cols = grep(paste0("^", cov_name), dm_colnames)
    }

    for (dm_col in dm_cols) {
      miss_idx_list[[length(miss_idx_list) + 1]] = c(
        row_offset + row_in_data,
        dm_col,
        which(all_missing_covs == cov_name)
      )
    }
  }

  if (length(miss_idx_list) > 0) {
    miss_idx = do.call(rbind, miss_idx_list)
    colnames(miss_idx) = c("row", "dm_col", "cov_idx")
  } else {
    miss_idx = matrix(integer(0), ncol = 3)
  }

  return(miss_idx)
}
