#' Estimate the logarithm of the normalizing constant for normalized power prior (NPP) for one data set
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the normalizing
#' constant for the NPP for a fixed value of the power prior parameter \eqn{a_0 \in (0, 1)} for one data
#' set. The initial priors are independent normal priors on the regression coefficients and a half-normal
#' prior on the dispersion parameter (if applicable).
#'
#' @include data_checks.R
#' @include expfam_loglik.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param histdata          a `data.frame` giving the historical data.
#' @param a0                the power prior parameter (a scalar between 0 and 1).
#' @param offset0           vector whose dimension is equal to the rows of the historical data set giving an offset for the
#'                          historical data. Defaults to a vector of 0s.
#' @param beta.mean         a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the mean parameters for the normal initial prior on regression coefficients given the dispersion
#'                          parameter. If a scalar is provided, beta.mean will be a vector of repeated elements of the given
#'                          scalar. Defaults to a vector of 0s.
#' @param beta.sd           a scalar or a vector whose dimension is equal to the number of regression coefficients giving
#'                          the sd parameters for the initial prior on regression coefficients. The sd used is
#'                          \code{sqrt(dispersion) * beta.sd}. If a scalar is provided, same as for beta.mean. Defaults to
#'                          a vector of 10s.
#' @param disp.mean         mean parameter for the half-normal prior on dispersion parameter. Defaults to 0.
#' @param disp.sd           sd parameter for the half-normal prior on dispersion parameter. Defaults to 10.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                         [cmdstanr::sample()].
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in [cmdstanr::sample()].
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in [cmdstanr::sample()].
#' @param ...               arguments passed to [cmdstanr::sample()] (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns a vector giving the value of a0, the estimated logarithm of the normalizing constant, the minimum
#'  estimated bulk effective sample size of the MCMC sampling, and the maximum Rhat.
#'
#' @references
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg036 = actg036[1:50, ]
#'   glm.npp.lognc(
#'     cd4 ~ treatment + age + race,
#'     family = poisson(), histdata = actg036, a0 = 0.5,
#'     chains = 1, iter_warmup = 500, iter_sampling = 5000
#'   )
#' }
glm.npp.lognc = function(
    formula,
    family,
    histdata,
    a0,
    offset0           = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  if( !( is.null(offset0) ) ){
    data.checks(formula, family, list(histdata), list(offset0))
  }else{
    data.checks(formula, family, list(histdata), NULL)
  }

  y0 = histdata[, all.vars(formula)[1]]
  n0 = length(y0)
  X0 = model.matrix(formula, histdata)
  p  = ncol(X0)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  if (length(a0) != 1)
    stop('a0 must be a scalar')
  if ( a0 < 0 | a0 > 1 )
    stop("a0 must be between 0 and 1")

  ## Default offset is vector of 0s
  if ( is.null(offset0) )
    offset0 = rep(0, n0)

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
    'n0'          = n0,
    'p'           = p,
    'y0'          = y0,
    'X0'          = X0,
    'beta_mean'   = beta.mean,
    'beta_sd'     = beta.sd,
    'disp_mean'   = disp.mean,
    'disp_sd'     = disp.sd,
    'a0'          = a0,
    'dist'        = dist,
    'link'        = link,
    'offs0'       = offset0
  )

  glm_npp_prior = instantiate::stan_package_model(
    name = "glm_npp_prior",
    package = "hdbayes"
  )

  ## fit model in cmdstanr
  fit  = glm_npp_prior$sample(data = standat,
                              iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains,
                              ...)
  d    = fit$draws(format = 'draws_matrix')
  summ = posterior::summarise_draws(d)

  ## estimate log normalizing constant
  log_density = function(pars, data){
    beta       = pars[paste0("beta[", 1:data$p,"]")]
    prior_lp   = sum( dnorm(beta, mean = data$beta_mean, sd = data$beta_sd, log = T) )
    dist       = data$dist
    link       = data$link
    dispersion = 1.0
    if ( dist > 2 ) {
      dispersion = pars[["dispersion[1]"]]
      prior_lp   = prior_lp +
        dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) -
        pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T)
    }
    hist_lp = glm_lp(data$y0, beta, data$X0, dist, link, data$offs0, dispersion)
    return(data$a0 * hist_lp + prior_lp)
  }

  post_samples = as.matrix(d[, -1, drop=F])
  lb           = rep(-Inf, p)
  ub           = rep(Inf, p)
  if( dist > 2 ) {
    lb = c(lb, 0)
    ub = c(ub, Inf)
  }
  names(ub) = colnames(post_samples)
  names(lb) = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = post_samples,
        'log_posterior' = log_density,
        'data' = standat,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## Return vector of a0, lognc, min_n_eff, max_Rhat
  res = c(
    'a0'           = a0,
    'lognc'        = bs$logml,
    'min_ess_bulk' = min(summ[, 'ess_bulk']),
    'max_Rhat'     = max(summ[, 'rhat'])
  )

  if ( res['min_ess_bulk'] < 1000 )
    warning(
      paste0(
        'The minimum bulk effective sample size of the MCMC sampling is ',
        res['min_ess_bulk'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have at least 1000. Try increasing the number of iterations.'
      )
    )
  if ( res['max_Rhat'] > 1.10 )
    warning(
      paste0(
        'The maximum Rhat of the MCMC sampling is ',
        res['max_Rhat'],
        ' for a0 = ',
        round(a0, 4),
        '. It is recommended to have a maximum Rhat of no more than 1.1. Try increasing the number of iterations.'
      )
    )
  return(res)
}
