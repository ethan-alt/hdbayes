#' Estimate the logarithm of the normalizing constant for latent exchangeability prior (LEAP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the latent exchangeability
#' prior (LEAP) using historical data sets.
#'
#' @include expfam_loglik.R
#' @include mixture_loglik.R
#'
#' @noRd
#'
#' @param post.samples      samples from the latent exchangeability prior (LEAP), with an attribute called 'data' which
#'                          includes the list of variables specified in the data block of the Stan program.
#' @param bridge.args       a `list` giving arguments (other than `samples`, `log_posterior`, `data`, `lb`, and `ub`)
#'                          to pass onto [bridgesampling::bridge_sampler()].
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{lognc}{the estimated logarithm of the normalizing constant}
#'
#'    \item{bs}{an object of class `bridge` or `bridge_list` giving the output from [bridgesampling::bridge_sampler()]}
#'  }
#'
#' @references
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2023). LEAP: The latent exchangeability prior for borrowing information from historical data. arXiv preprint.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036)
#'   ## take subset for speed purposes
#'   actg036   = actg036[1:50, ]
#'   formula   = outcome ~ scale(age) + race + treatment + scale(cd4)
#'   family    = binomial('logit')
#'   hist.data.list = list(histdata = actg036)
#'   hist.standat = hdbayes:::get.stan.data.leap(
#'     formula = formula,
#'     family = family,
#'     data.list = hist.data.list,
#'     K = 2,
#'     all.hist = T
#'   )
#'   glm_leap_prior = instantiate::stan_package_model(
#'     name = "glm_leap_prior",
#'     package = "hdbayes"
#'   )
#'   fit = glm_leap_prior$sample(
#'     data = hist.standat,
#'     iter_warmup = 500, iter_sampling = 1000, chains = 1
#'   )
#'   d.leap.hist  = fit$draws(format = 'draws_df')
#'   attr(x = d.leap.hist, which = 'data') = hist.standat
#'   glm.leap.lognc(
#'     post.samples = d.leap.hist
#'   )
#' }
glm.leap.lognc = function(
    post.samples,
    bridge.args       = NULL
) {
  ## get Stan data for LEAP
  stan.data = attr(post.samples, 'data')

  p        = stan.data$p
  K        = stan.data$K
  oldnames = paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  if ( stan.data$dist > 2 ) {
    oldnames = c(oldnames, paste0( 'dispersion[', 1:K, ']' ))
  }
  oldnames = c(oldnames, "gamma")
  if ( K > 2 ){
    oldnames = c(oldnames, paste0("delta_raw[", 1:(K-2), "]"))
  }
  d = suppressWarnings(
    as.matrix(post.samples[, oldnames, drop=F])
  )

  ## compute log normalizing constants for half-normal priors
  stan.data$lognc_disp  = sum( pnorm(0, mean = stan.data$disp_mean, sd = stan.data$disp_sd, lower.tail = F, log.p = T) )

  ## compute log normalizing constants for gamma
  gamma_shape1    = stan.data$conc[1]
  gamma_shape2    = sum(stan.data$conc[2:K])

  stan.data$lognc_gamma = 0
  if( stan.data$gamma_lower != 0 || stan.data$gamma_upper != 1 ) {
    stan.data$lognc_gamma = log( pbeta(stan.data$gamma_upper, shape1 = gamma_shape1, shape2 = gamma_shape2) -
                                         pbeta(stan.data$gamma_lower, shape1 = gamma_shape1, shape2 = gamma_shape2) )
  }

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    betaMat    = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    prior_lp   = sum( dnorm(betaMat, mean = as.numeric(data$mean_beta),
                            sd = as.numeric(data$sd_beta), log = T) )
    betaMat    = matrix(betaMat, nrow = p, ncol = K)

    ## prior on gamma
    conc         = data$conc
    gamma_shape1 = conc[1]
    gamma_shape2 = sum(conc[2:K])
    gamma        = pars[["gamma"]]
    probs        = c(gamma, 1 - gamma)

    if ( gamma_shape1 != 1 || gamma_shape2 != 1 ){
      prior_lp = prior_lp + dbeta(gamma, gamma_shape1, gamma_shape2, log = T)
    }
    prior_lp     = prior_lp - data$lognc_gamma

    if( K > 2 ){
      delta_raw = pars[paste0("delta_raw[", 1:(K-2), "]")]
      delta_raw = c(delta_raw, 1 - sum(delta_raw))
      prior_lp  = prior_lp + dirichlet_lp(delta_raw, conc[2:K])
      probs     = c(gamma, (1 - gamma) * delta_raw)
    }

    dist         = data$dist
    link         = data$link
    if ( dist > 2 ) {
      dispersion = pars[paste0("dispersion[", 1:K,"]")]
      prior_lp   = prior_lp +
        sum( dnorm(dispersion, mean = data$disp_mean, sd = data$disp_sd, log = T) ) - data$lognc_disp
    }else {
      dispersion = rep(1.0, K)
    }
    hist_lp      = glm_mixture_lp(data$y0, betaMat, dispersion, probs, data$X0, dist, link, data$offs0)
    return(hist_lp + prior_lp)
  }

  lb = rep(-Inf, p*K)
  ub = rep(Inf, p*K)
  if( stan.data$dist > 2 ) {
    lb = c(lb, rep(0, K) )
    ub = c(ub, rep(Inf, K) )
  }
  lb = c(lb, stan.data$gamma_lower)
  ub = c(ub, stan.data$gamma_upper)
  if( K > 2 ){
    lb = c(lb, rep(0, K-2))
    ub = c(ub, rep(1, K-2))
  }
  names(ub) = colnames(d)
  names(lb) = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples"       = d,
        'log_posterior' = log_density,
        'data'          = stan.data,
        'lb'            = lb,
        'ub'            = ub
        ),
      bridge.args
    )
  )

  ## Return a list of lognc and output from bridgesampling::bridge_sampler
  res = list(
    'lognc'        = bs$logml,
    'bs'           = bs
  )
  return(res)
}
