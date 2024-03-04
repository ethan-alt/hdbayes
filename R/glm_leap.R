#' Posterior of Latent Exchangeability Prior (LEAP)
#'
#' Sample from the posterior distribution of a GLM using the LEAP by Alt et al. (2023).
#'
#' The LEAP discounts the historical data by identifying the most relevant individuals from the historical data.
#' It is equivalent to a prior induced by the posterior of a finite mixture model for the historical data set.
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of two `data.frame`s giving the current data followed by one historical data set.
#' @param K                 the desired number of classes to identify. Defaults to 2.
#' @param prob.conc         a scalar or a vector of length `K` giving the concentration parameters for Dirichlet prior.
#'                          If length == 2, a `beta(prob.conc[1], prob.conc[2])` prior is used. If a scalar is provided,
#'                          prob.conc will be a vector of repeated elements of the given scalar. Defaults to a vector of 1s.
#' @param beta.mean         a `p x K` matrix of mean parameters for initial prior on regression coefficients (including
#'                          intercept). Defaults to a matrix of 0s.
#' @param beta.sd           a `p x K` matrix of sd parameters for the initial prior on regression coefficients (including
#'                          intercept). Defaults to a matrix of 10s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the means for
#'                          the half-normal hyperpriors on the dispersion parameters. If a scalar is provided, disp.mean will
#'                          be a vector of repeated elements of the given scalar. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of classes (`K`) giving the sds for
#'                          the half-normal hyperpriors on the dispersion parameters. If a scalar is provided, same as for
#'                          disp.mean. Defaults to a vector of 10s.
#' @param offset.list       a list of matrices giving the offset for current data followed by historical data. For each
#'                          matrix, the number of rows corresponds to observations and columns correspond to classes.
#'                          Defaults to a list of matrices of 0s.
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          [cmdstanr::sample()].
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in [cmdstanr::sample()].
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in [cmdstanr::sample()].
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns an object of class `draws_df` giving posterior samples.
#'
#' @references
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2023). LEAP: The latent exchangeability prior for borrowing information from historical data. arXiv preprint.
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' # take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' if (instantiate::stan_cmdstan_exists()) {
#'   glm.leap(
#'     formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'     family = binomial('logit'),
#'     data.list = list(actg019, actg036),
#'     K = 2,
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#' }
glm.leap = function(
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
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  data.checks(formula, family, data.list, offset.list)

  ## get model information
  data     = data.list[[1]]
  ## restricted to one historical data set
  histdata = data.list[[2]]
  offset   = offset.list[[1]]
  offset0  = offset.list[[2]]

  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]
  n  = length(y)
  n0 = length(y0)
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)
  p  = ncol(X)
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
  if ( is.null(offset) )
    offset = matrix(rep(0, n*K), ncol=K)
  if ( is.null(offset0) )
    offset0 = matrix(rep(0, n0*K), ncol=K)

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( is.null(beta.mean) )
    beta.mean = matrix(rep(0, p*K), ncol=K)
  if ( is.null(beta.sd) )
    beta.sd   = matrix(rep(10, p*K), ncol=K)

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

  glm_leap        = instantiate::stan_package_model(
    name = "glm_leap",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit = glm_leap$sample(data = standat,
                       iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                       ...)
  pars = fit$metadata()$model_params

  ## rename parameters
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)

  ## reorder parameters so that regression coefficients appear at the top
  pars = c(pars[1], pars[pars %in% oldnames], (pars[!pars %in% oldnames])[-1])
  d   = fit$draws(format = 'draws_df', variables = pars)
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}
