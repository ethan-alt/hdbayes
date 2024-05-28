#' Log marginal likelihood of a GLM under commensurate prior (CP)
#'
#' Uses Markov chain Monte Carlo (MCMC) and bridge sampling to estimate the logarithm of the marginal
#' likelihood of a GLM under CP.
#'
#' This function shares the same arguments as [glm.commensurate()], while introducing two additional
#' parameters: `post.samples` and `bridge.args`. `post.samples` provides posterior samples under CP
#' (e.g., the output from [glm.commensurate()]), whereas `bridge.args` specifies arguments to pass onto
#' [bridgesampling::bridge_sampler()] (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`).
#'
#' It is important to ensure that the values assigned to the shared arguments (excluding those relevant
#' for MCMC sampling) in this function and [glm.commensurate()] match with those used in generating
#' `post.samples`. The arguments pertinent to MCMC sampling are utilized to compute the normalizing
#' constants for CP.
#'
#' @include get_stan_data.R
#' @include data_checks.R
#' @include expfam_loglik.R
#' @include glm_commensurate_lognc.R
#'
#' @export
#'
#' @inheritParams glm.commensurate
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under CP, such as the output from [glm.commensurate()]. Each row corresponds
#'                          to the posterior samples obtained from one iteration of MCMC. The column names of `post.samples`
#'                          should include the names of covariates for regression coefficients, such as "(Intercept)",
#'                          and "dispersion" for the dispersion parameter, if applicable.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{model}{"Commensurate"}
#'
#'    \item{logml}{the estimated logarithm of the marginal likelihood}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using all data sets (including current and historical data)}
#'
#'    \item{bs.hist}{an object of class `bridge` or `bridge_list` giving the output from
#'    [bridgesampling::bridge_sampler()] using historical data sets (for computing the log normalizing
#'    constant for CP)}
#'
#'  }
#'
#' @references
#'  Hobbs, B. P., Carlin, B. P., Mandrekar, S. J., and Sargent, D. J. (2011). Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. Biometrics, 67(3), 1047–1056.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg019)
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg019   = actg019[1:100, ]
#'   actg036   = actg036[1:50, ]
#'   formula   = cd4 ~ treatment + age + race
#'   family    = poisson()
#'   data_list = list(currdata = actg019, histdata = actg036)
#'   d.cp      = glm.commensurate(
#'     formula = formula,
#'     family = family,
#'     data.list = data_list,
#'     tau = rep(5, 4),
#'     chains = 1, iter_warmup = 500, iter_sampling = 1000
#'   )
#'   glm.logml.commensurate(
#'     formula = formula, family = family,
#'     data.list = data_list,
#'     tau = rep(5, 4),
#'     post.samples = d.cp,
#'     bridge.args = list(silent = TRUE),
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )
#' }
glm.logml.commensurate = function(
    post.samples,
    formula,
    family,
    data.list,
    tau,
    offset.list       = NULL,
    beta0.mean        = NULL,
    beta0.sd          = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  K = length(data.list)
  if ( K == 1 ){
    stop("data.list should include at least one historical data set")
  }

  ## get Stan data for CP
  standat = get.stan.data.cp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    tau            = tau,
    offset.list    = offset.list,
    beta0.mean     = beta0.mean,
    beta0.sd       = beta0.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd
  )

  ## check the format of post.samples
  post.samples.checks(post.samples, colnames(standat$X), family)

  d        = as.matrix(post.samples)
  ## rename parameters
  p        = standat$p
  X        = standat$X
  oldnames = c(paste0("beta[", 1:p, "]"), paste0("beta0[", 1:p, "]"))
  newnames = c(colnames(X), paste0( colnames(X), '_hist') )
  colnames(d)[colnames(d) %in% newnames] = oldnames
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion', paste0( 'dispersion', '_hist_', 1:(K-1) ))
  }
  d = d[, oldnames, drop=F]

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    N          = data$N
    beta       = pars[paste0("beta[", 1:p,"]")]
    beta0      = pars[paste0("beta0[", 1:p,"]")]
    prior_lp   = sum( dnorm(beta0, mean = data$beta0_mean, sd = data$beta0_sd, log = T) ) +
      sum( dnorm(beta, mean = beta0, sd = 1/sqrt(data$tau), log = T) )
    dist       = data$dist
    link       = data$link
    start.idx  = data$start_idx
    end.idx    = data$end_idx
    dispersion = rep(1.0, K)
    if ( dist > 2 ){
      ## prior on dispersion
      dispersion = pars[c("dispersion", paste0( "dispersion", "_hist_", 1:(K-1) ))]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) -
        sum( pnorm(0, mean = data$disp_mean, sd = data$disp_sd, lower.tail = F, log.p = T) )
      ## historical data likelihood
      prior_lp    = prior_lp + sum( sapply(2:K, function(k){
        y         = data$y[ start.idx[k]:end.idx[k] ]
        X         = data$X[ start.idx[k]:end.idx[k], ]
        offs      = data$offs[ start.idx[k]:end.idx[k] ]
        glm_lp(y, beta0, X, dist, link, offs, dispersion[k])
      }) )
    }else {
      ## historical data likelihood
      prior_lp    = prior_lp +
        glm_lp(data$y[ start.idx[2]:N ], beta0,
               data$X[ start.idx[2]:N, ], dist, link,
               data$offs[ start.idx[2]:N ], 1.0)
    }
    ## current data likelihood
    y         = data$y[ start.idx[1]:end.idx[1] ]
    X         = data$X[ start.idx[1]:end.idx[1], ]
    offs      = data$offs[ start.idx[1]:end.idx[1] ]
    data_lp   = glm_lp(y, beta, X, dist, link, offs, dispersion[1])
    return(data_lp + prior_lp)
  }

  lb = rep(-Inf, p*2)
  ub = rep(Inf, p*2)
  if( standat$dist > 2 ) {
    lb = c(lb, rep(0, K))
    ub = c(ub, rep(Inf, K))
  }
  names(ub) = colnames(d)
  names(lb) = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = d,
        'log_posterior' = log_density,
        'data' = standat,
        'lb' = lb,
        'ub' = ub),
      bridge.args
    )
  )

  ## estimate log normalizing constant for CP using historical data sets
  if ( !is.null(disp.mean) && (length(disp.mean) != 1) ){
    disp.mean = disp.mean[-1]
  }
  if ( !is.null(disp.sd) && (length(disp.sd) != 1) ){
    disp.sd = disp.sd[-1]
  }
  res.hist = glm.commensurate.lognc(
    formula           = formula,
    family            = family,
    hist.data.list    = data.list[-1],
    tau               = tau,
    hist.offset.list  = offset.list[-1],
    beta0.mean        = beta0.mean,
    beta0.sd          = beta0.sd,
    hist.disp.mean    = disp.mean,
    hist.disp.sd      = disp.sd,
    bridge.args       = bridge.args,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    chains            = chains,
    ...
  )

  ## Return a list of model name, estimated log marginal likelihood, and output from bridgesampling::bridge_sampler
  res = list(
    'model'   = "Commensurate",
    'logml'   = bs$logml - res.hist$lognc,
    'bs'      = bs,
    'bs.hist' = res.hist$bs
  )
}
