#' Estimate the logarithm of the normalizing constant for stratified power prior (PP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the stratified power
#' prior (PP) using all data sets or using historical data sets only. Note that the power prior parameters
#' (\eqn{a_0}'s) are treated as fixed.
#'
#' @include aft_loglik.R
#' @include mixture_aft_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of an AFT model under the stratified power prior (PP) or samples from
#'                          the stratified PP, with an attribute called 'data' which includes the list of variables
#'                          specified in the data block of the Stan program.
#' @param is.prior          whether the samples are from the stratified PP (using historical data sets only).
#'                          Defaults to FALSE.
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
#'  Wang, C., Li, H., Chen, W.-C., Lu, N., Tiwari, R., Xu, Y., & Yue, L. Q. (2019). Propensity score-integrated power prior approach for incorporating real-world evidence in single-arm clinical studies. Journal of Biopharmaceutical Statistics, 29(5), 731–748.
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
#'     data_list = list(currdata = E1690, histdata = E1684)
#'     strata_list = list(rep(1:2, each = 25), rep(1:2, each = 50))
#'     # Alternatively, we can determine the strata based on propensity scores
#'     # using the psrwe package, which is available on GitHub
#'     d.stratified.pp = aft.stratified.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       strata.list = strata_list,
#'       a0.strata = c(0.3, 0.5),
#'       dist = "weibull",
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     aft.stratified.pp.lognc(
#'       post.samples = d.stratified.pp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
aft.stratified.pp.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X         = stan.data$X0_obs
  K         = stan.data$K
  oldnames  = c( paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"),
                 paste0("scaleVec[", 1:K, "]") )
  newnames  = c( paste0( colnames(X), '_stratum_', rep(1:K, each = p) ),
                 paste0("scale_stratum_", 1:K) )
  colnames(d)[colnames(d) %in% newnames] = oldnames
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on scale
  stan.data$scale_prior_lognc = pnorm(0, mean = stan.data$scale_mean, sd = stan.data$scale_sd, lower.tail = F, log.p = T)
  stan.data$is_prior          = is.prior

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta       = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta       = matrix(beta, nrow = p, ncol = K)
    scale      = as.numeric( pars[paste0("scaleVec[", 1:K, "]")] )
    prior_lp   = sum( sapply(1:K, function(k){
      sum(dnorm(beta[, k], mean = data$beta_mean, sd = data$beta_sd, log = T)) +
        dnorm(scale[k], mean = data$scale_mean, sd = data$scale_sd, log = T) - data$scale_prior_lognc
    })
    )

    Eta0_obs       = data$X0_obs %*% beta
    Eta0_cen       = data$X0_cen %*% beta
    stratumID0_obs = data$stratumID0_obs
    stratumID0_cen = data$stratumID0_cen
    y0_obs         = data$y0_obs
    y0_cen         = data$y0_cen
    a0s            = data$a0s

    eta0_obs = sapply(1:length(stratumID0_obs), function(i){
      Eta0_obs[i, stratumID0_obs[i]]
    })
    eta0_cen = sapply(1:length(stratumID0_cen), function(i){
      Eta0_cen[i, stratumID0_cen[i]]
    })
    data_lp = sum( a0s[stratumID0_obs] * aft_model_obs_lpdf(y0_obs, eta0_obs, scale[stratumID0_obs], data$dist) ) +
      sum( a0s[stratumID0_cen] * aft_model_cen_lpdf(y0_cen, eta0_cen, scale[stratumID0_cen], data$dist) )

    if( !data$is_prior ){
      Eta_obs        = data$X_obs %*% beta
      Eta_cen        = data$X_cen %*% beta
      stratumID_obs  = data$stratumID_obs
      stratumID_cen  = data$stratumID_cen
      y_obs          = data$y_obs
      y_cen          = data$y_cen

      eta_obs        = sapply(1:length(stratumID_obs), function(i){
        Eta_obs[i, stratumID_obs[i]]
      })
      eta_cen        = sapply(1:length(stratumID_cen), function(i){
        Eta_cen[i, stratumID_cen[i]]
      })
      data_lp        = data_lp + sum( aft_model_obs_lpdf(y_obs, eta_obs, scale[stratumID_obs], data$dist) ) +
        sum( aft_model_cen_lpdf(y_cen, eta_cen, scale[stratumID_cen], data$dist) )
    }
    return(data_lp + prior_lp)
  }

  lb           = c(rep(-Inf, p*K), rep(0, K))
  ub           = rep(Inf, length(lb))
  names(ub)    = colnames(d)
  names(lb)    = names(ub)

  bs = do.call(
    what = bridgesampling::bridge_sampler,
    args = append(
      list(
        "samples" = d,
        'log_posterior' = log_density,
        'data' = stan.data,
        'lb' = lb,
        'ub' = ub),
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
