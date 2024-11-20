#' Estimate the logarithm of the normalizing constant for latent exchangeability prior (LEAP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the latent exchangeability
#' prior (LEAP) using historical data set.
#'
#' @include aft_loglik.R
#' @include mixture_aft_loglik.R
#' @include mixture_loglik.R
#' @include expfam_loglik.R
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
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1684)
#'     data(E1690)
#'     ## take subset for speed purposes
#'     E1684 = E1684[1:100, ]
#'     E1690 = E1690[1:50, ]
#'     ## replace 0 failure times with 0.50 days
#'     E1684$failtime[E1684$failtime == 0] = 0.50/365.25
#'     E1690$failtime[E1690$failtime == 0] = 0.50/365.25
#'     E1684$cage = as.numeric(scale(E1684$age))
#'     E1690$cage = as.numeric(scale(E1690$age))
#'     data_list = list(currdata = E1690, histdata = E1684)
#'     d.leap = aft.leap(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       K= 2,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     aft.leap.lognc(
#'       post.samples = d.leap,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
aft.leap.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  ## get Stan data for LEAP
  stan.data = attr(post.samples, 'data')

  p        = stan.data$p
  K        = stan.data$K
  oldnames = paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")
  oldnames = c(oldnames, paste0( 'scaleVec[', 1:K, ']' ))
  oldnames = c(oldnames, "logit_gamma")
  if ( K > 2 ){
    oldnames = c(oldnames, paste0("delta_raw[", 1:(K-2), "]"))
  }
  d = suppressWarnings(
    as.matrix(post.samples[, oldnames, drop=F])
  )

  ## compute log normalizing constants for half-normal priors
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)

  ## compute log normalizing constants for gamma
  gamma_shape1    = stan.data$prob_conc[1]
  gamma_shape2    = sum(stan.data$prob_conc[2:K])

  stan.data$lognc_logit_gamma = 0
  if( stan.data$gamma_lower != 0 || stan.data$gamma_upper != 1 ) {
    stan.data$lognc_logit_gamma = log( pbeta(stan.data$gamma_upper, shape1 = gamma_shape1, shape2 = gamma_shape2) -
                                         pbeta(stan.data$gamma_lower, shape1 = gamma_shape1, shape2 = gamma_shape2) )
  }

  stan.data$is_prior = is.prior

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    betaMat    = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    betaMat    = matrix(betaMat, nrow = p, ncol = K)
    scaleVec   = as.numeric( pars[paste0('scaleVec[', 1:K, ']')] )

    prior_lp   = 0
    for( k in 1:K ){
      prior_lp = prior_lp + sum( dnorm(betaMat[, k], mean = as.numeric(data$beta_mean),
                                       sd = as.numeric(data$beta_sd), log = T) )
      prior_lp = prior_lp + as.numeric( dnorm(scaleVec[k], mean = data$scale_mean, sd = data$scale_sd, log = T) ) - data$scale_prior_lognc
    }

    ## prior on logit(gamma)
    conc         = data$prob_conc
    gamma_shape1 = conc[1]
    gamma_shape2 = sum(conc[2:K])
    logit_gamma  = pars[["logit_gamma"]]
    log1m_gamma  = -log1p_exp(logit_gamma) # log(1 - gamma)
    log_probs    = c(logit_gamma, 0) + log1m_gamma

    prior_lp     = prior_lp + logit_beta_lp(logit_gamma, gamma_shape1, gamma_shape2) -
      data$lognc_logit_gamma

    if( K > 2 ){
      delta_raw = as.numeric(pars[paste0("delta_raw[", 1:(K-2), "]")])
      delta_raw = c(delta_raw, 1 - sum(delta_raw))
      prior_lp  = prior_lp + dirichlet_lp(delta_raw, conc[2:K])
      log_probs = c(logit_gamma, log(delta_raw)) + log1m_gamma
    }

    Eta0_obs    = data$X0_obs %*% betaMat
    Eta0_cen    = data$X0_cen %*% betaMat
    y0_obs      = data$y0_obs
    y0_cen      = data$y0_cen
    dist        = data$dist
    data_lp     = sum ( sapply(1:data$n0_obs, function(i){
      aft_model_obs_mixture_lp(y0_obs[i], Eta0_obs[i, ], scaleVec, log_probs, dist)
    }) )
    data_lp     = data_lp + sum(
      sapply(1:data$n0_cen, function(i){
        aft_model_cen_mixture_lp(y0_cen[i], Eta0_cen[i, ], scaleVec, log_probs, dist)
      })
    )

    if( !data$is_prior ){
      eta_obs = data$X_obs %*% betaMat[, 1]
      eta_cen = data$X_cen %*% betaMat[, 1]
      data_lp = data_lp + sum( aft_model_lp(data$y_obs, data$y_cen, eta_obs, eta_cen, scaleVec[1], dist) )
    }

    return(data_lp + prior_lp)
  }

  lb = c(rep(-Inf, p*K), rep(0, K), binomial('logit')$linkfun(stan.data$gamma_lower))
  ub = c(rep(Inf, (p+1)*K), binomial('logit')$linkfun(stan.data$gamma_upper))
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
