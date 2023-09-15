
#'
#' Normalized asymptotic power prior (NAPP)
#'
#' Sample from the posterior distribution of a GLM using the normalized
#' asymptotic power prior. The regression coefficients and logarithm of the
#' dispersion parameter are a multivariate normal distribution with mean
#' equal to the maximum likelihood estimate of the historical data and
#' covariance matrix equal to \eqn{a0^{-1}} multiplied by the inverse Fisher
#' information matrix of the historical data, where a0 is the power prior
#' parameter (treated as random).
#'
#' @include data_checks.R
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
#' ## take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' data.list = list(data = actg019, histdata = actg036)
#' glm.napp(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(), data.list = data.list,
#'   chains = 1, iter_warmup = 500, iter_sampling = 1000
#' )
glm.napp = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    local.location    = NULL,
    ...
) {
  data.checks(formula, family, data.list, offset.list)

  K        = length(data.list) - 1 # number of historical datasets
  data     = data.list[[1]] # current data
  y        = data[, all.vars(formula)[1]]
  n        = length(y)
  X        = model.matrix(formula, data)
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

  model_name      = "glm_napp"
  model_file_path = file.path(local.location, paste0(model_name, ".stan"))
  glm_napp        = cmdstanr::cmdstan_model(model_file_path)

  ## fit model in cmdstanr
  fit = glm_napp$sample(data = standat, ...)
  d   = fit$draws(format = 'draws_df')


  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0s[', 1:K, ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:K))
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames

  return(d)
}
