
#'
#' Posterior of normalized power prior
#'
#' Sample from the posterior distribution of a GLM using the normalized power prior
#' (NPP). Before using this function, users must estimate the logarithm of the
#' normalizing constant across a range of power prior parameters (a0), possibly
#' smoothing techniques over a find grid.
#'
#'
#' @include data_checks.R
#' @include glm_npp_lognc.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rests
#'                          are the historical datasets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                          the length of data.list. The length of each element of offset.list is equal to the number
#'                          of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          beta.mean will be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. If a scalar is provided,
#'                          same as for beta.mean. Defaults to a vector of 10s.
#' @param disp.mean         mean parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           sd parameter for the half-normal prior on dispersion parameter. Defaults to 10.
#' @param a0.lognc          a vector giving values of the power prior parameter for which the logarithm of the normalizing
#'                          constant has been evaluated
#' @param lognc             an S by T matrix where S is the length of a0.lognc, T is the number of historical datasets, and
#'                          the j-th column, j = 1, ..., T, is a vector giving the logarithm of the normalizing constant (as
#'                          estimated by \code{\link[hdbayes]{glm.npp.lognc}}) for a0.lognc using the j-th historical dataset.
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
#' \dontrun{
#'   data(actg019)
#'   data(actg036)
#'   library(parallel)
#'   ncores = 5
#'
#'   data.list = list(data = actg019, histdata = actg036)
#'   formula   = cd4 ~ treatment + age + race
#'   family    = poisson()
#'   a0        = seq(0, 1, length.out = 100)
#'   a0.lognc  = list()
#'   a0.lognc.hdbayes = data.frame( a0 = a0 )
#'   ## call created function
#'   for (i in 2:length(data.list)) {
#'     histdata = data.list[[i]]
#'     ## wrapper to obtain log normalizing constant in parallel package
#'     logncfun = function(a0, ...){
#'       glm.npp.lognc(
#'         formula = formula, family = family, a0 = a0, histdata = histdata, ...
#'       )
#'     }
#'
#'     cl = makeCluster(ncores)
#'     clusterSetRNGStream(cl, 123)
#'     clusterExport(cl, varlist = c('formula', 'family', 'histdata', 'glm.npp.lognc',
#'                                   'glm_lp', 'get_lp2mean', 'normal_glm_lp', 'bernoulli_glm_lp',
#'                                   'poisson_glm_lp', 'gamma_glm_lp', 'invgauss_glm_lp'))
#'     a0.lognc[[i-1]] = parLapply(
#'       cl = cl, X = a0, fun = logncfun, iter_warmup = 2000,
#'       iter_sampling = 5000, chains = 1, refresh = 0
#'     )
#'     stopCluster(cl)
#'     a0.lognc[[i-1]] = data.frame( do.call(rbind, a0.lognc[[i-1]]) )
#'     a0.lognc.hdbayes = cbind(a0.lognc.hdbayes, a0.lognc[[i-1]]$lognc)
#'     colnames(a0.lognc.hdbayes)[i] = paste0("lognc_hist", i-1)
#'     }
#'   a0.lognc = a0.lognc.hdbayes$a0
#'   lognc    = as.matrix(a0.lognc.hdbayes[, -1, drop = F])
#'
#'   ## sample from normalized power prior
#'   glm.npp(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson(),
#'     data.list = data.list,
#'     a0.lognc = a0.lognc,
#'     lognc = lognc,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000,
#'     refresh = 0
#'   )
#' }
#'
glm.npp = function(
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
    local.location    = NULL,
    ...
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

  ## Default offset for each dataset is a vector of 0s
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
    stop('the number of columns in lognc must be the same as the number of historical datasets')
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

  model_name        = "glm_npp_posterior"
  model_file_path   = file.path(local.location, paste0(model_name, ".stan"))
  glm_npp_posterior = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = glm_npp_posterior$sample(data = standat, ...)
  d   = fit$draws(format = 'draws_df')


  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0_vals[', 1:(K-1), ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:(K-1)))
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
