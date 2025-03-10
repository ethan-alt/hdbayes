#' Estimate the logarithm of the normalizing constant for latent exchangeability prior (LEAP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the latent exchangeability
#' prior (LEAP) using historical data set.
#'
#' @include pwe_loglik.R
#' @include mixture_loglik.R
#' @include expfam_loglik.R
#'
#' @noRd
#'
#' @param post.samples      samples from the latent exchangeability prior (LEAP), with an attribute called 'data' which
#'                          includes the list of variables specified in the data block of the Stan program.
#' @param is.prior          whether the samples are from the LEAP (using historical data set only). Defaults to FALSE.
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
#'  Alt, E. M., Chang, X., Jiang, X., Liu, Q., Mo, M., Xia, H. M., and Ibrahim, J. G. (2024). LEAP: The latent exchangeability prior for borrowing information from historical data. Biometrics, 80(3).
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     d.leap = curepwe.leap(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       K= 2,
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     curepwe.leap.lognc(
#'       post.samples = d.leap,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
curepwe.leap.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  ## get Stan data for LEAP
  stan.data = attr(post.samples, 'data')

  p        = stan.data$p
  J        = stan.data$J
  K        = stan.data$K
  oldnames = c(paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"),
               paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]"),
               "logit_gamma")
  if ( is.prior ){
    oldnames = c(oldnames, "logit_p_cured0")
  }else{
    oldnames = c(oldnames, "logit_p_cured", "logit_p_cured0")
  }
  if ( K > 2 ){
    oldnames = c(oldnames, paste0("delta_raw[", 1:(K-2), "]"))
  }
  d = suppressWarnings(
    as.matrix(post.samples[, oldnames, drop=F])
  )

  ## compute log normalizing constants (lognc) for half-normal prior on baseline hazards
  stan.data$lognc_hazard = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )

  ## compute log normalizing constants for gamma
  gamma_shape1    = stan.data$conc[1]
  gamma_shape2    = sum(stan.data$conc[2:K])

  stan.data$lognc_logit_gamma = 0
  if( stan.data$gamma_lower != 0 || stan.data$gamma_upper != 1 ) {
    stan.data$lognc_logit_gamma = log( pbeta(stan.data$gamma_upper, shape1 = gamma_shape1, shape2 = gamma_shape2) -
                                         pbeta(stan.data$gamma_lower, shape1 = gamma_shape1, shape2 = gamma_shape2) )
  }

  stan.data$is_prior = is.prior

  ## estimate log normalizing constant
  log_density = function(pars, data){
    p              = data$p
    J              = data$J
    K              = data$K
    betaMat        = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    betaMat        = matrix(betaMat, nrow = p, ncol = K)
    lambdaMat      = pars[paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]")]
    lambdaMat      = matrix(lambdaMat, nrow = J, ncol = K)
    logit_p_cured0 = as.numeric( pars[["logit_p_cured0"]] )
    log1m_p_cured0 = -log1p_exp(logit_p_cured0) # log(1 - p_cured0)
    log_p_cured0   = logit_p_cured0 + log1m_p_cured0 # log(p_cured0)

    prior_lp   = 0
    for( k in 1:K ){
      prior_lp = prior_lp + sum( dnorm(betaMat[, k], mean = as.numeric(data$beta_mean),
                                       sd = as.numeric(data$beta_sd), log = T) )
      prior_lp = prior_lp + sum( dnorm(lambdaMat[, k], mean = as.numeric(data$hazard_mean),
                                       sd = as.numeric(data$hazard_sd), log = T) ) - data$lognc_hazard
    }

    ## prior on logit_p_cured0
    prior_lp = prior_lp + dnorm(logit_p_cured0, mean = data$logit_p_cured_mean, sd = data$logit_p_cured_sd, log = T)

    ## prior on logit(gamma)
    conc         = data$conc
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

    eta0Mat  = data$X0 %*% betaMat

    # log likelihood contribution for the non-cured population
    contribs0 = sapply(1:K, function(k){
      log_probs[k] + pwe_lpdf(data$y0, eta0Mat[, k], lambdaMat[, k], data$breaks, data$intindx0, data$J, data$death_ind0)
    })
    contribs0   = as.matrix(contribs0, ncol = K)
    noncured_lp = apply(contribs0, 1, log_sum_exp)

    data_lp     = cbind(log_p_cured0 + log(1 - data$death_ind0),
                        log1m_p_cured0 + noncured_lp)
    data_lp     = sum( apply(data_lp, 1, log_sum_exp) )

    if( !data$is_prior ){
      eta           = data$X1 %*% betaMat[, 1]
      logit_p_cured = as.numeric( pars[["logit_p_cured"]] ) # logit of p_cured
      log1m_p_cured = -log1p_exp(logit_p_cured) # log(1 - p_cured)
      log_p_cured   = logit_p_cured + log1m_p_cured # log(p_cured)
      contribs      = cbind(log_p_cured + log(1 - data$death_ind),
                            log1m_p_cured + pwe_lpdf(data$y1, eta, lambdaMat[, 1], data$breaks, data$intindx, data$J, data$death_ind))
      data_lp       = data_lp + sum( apply(contribs, 1, log_sum_exp) )
    }

    return(data_lp + prior_lp)
  }

  lb = c(rep(-Inf, p*K), rep(0, J*K), binomial('logit')$linkfun(stan.data$gamma_lower), -Inf)
  ub = c(rep(Inf, (p+J)*K), binomial('logit')$linkfun(stan.data$gamma_upper), Inf)
  if ( !is.prior ){
    lb = c(lb, -Inf)
    ub = c(ub, Inf)
  }
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
