#' Estimate the logarithm of the normalizing constant for stratified power prior (PP)
#'
#' Uses bridge sampling to estimate the logarithm of the normalizing constant for the stratified power
#' prior (PP) using all data sets or using historical data set only. Note that the power prior parameters
#' (\eqn{a_0}'s) are treated as fixed.
#'
#' @include pwe_loglik.R
#'
#' @noRd
#'
#' @param post.samples      posterior samples of a PWE model under the stratified power prior (PP) or samples from
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
#'     nbreaks = 3
#'     probs   = 1:nbreaks / nbreaks
#'     breaks  = as.numeric(
#'       quantile(E1690[E1690$failcens==1, ]$failtime, probs = probs)
#'     )
#'     breaks  = c(0, breaks)
#'     breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
#'     d.stratified.pp = pwe.stratified.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment,
#'       data.list = data_list,
#'       strata.list = strata_list,
#'       breaks = breaks,
#'       a0.strata = c(0.3, 0.5),
#'       chains = 1, iter_warmup = 500, iter_sampling = 1000
#'     )
#'     pwe.stratified.pp.lognc(
#'       post.samples = d.stratified.pp,
#'       is.prior = FALSE,
#'       bridge.args = list(silent = TRUE)
#'     )
#'   }
#' }
pwe.stratified.pp.lognc = function(
    post.samples,
    is.prior          = FALSE,
    bridge.args       = NULL
) {
  stan.data = attr(post.samples, 'data')
  d         = as.matrix(post.samples)
  ## rename parameters
  p         = stan.data$p
  X0        = stan.data$X0
  J         = stan.data$J
  K         = stan.data$K
  oldnames = c(paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]"),
               paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]"))
  newnames = c(paste0( colnames(X0), '_stratum_', rep(1:K, each = p) ),
               paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"))
  colnames(d)[colnames(d) %in% newnames] = oldnames
  d = d[, oldnames, drop=F]

  ## compute log normalizing constants (lognc) for half-normal prior on baseline hazards
  stan.data$lognc_hazard = sum( pnorm(0, mean = stan.data$hazard_mean, sd = stan.data$hazard_sd, lower.tail = F, log.p = T) )
  stan.data$is_prior     = is.prior

  ## log of the unnormalized posterior density function
  log_density = function(pars, data){
    p          = data$p
    K          = data$K
    beta       = pars[paste0("betaMat[", rep(1:p, K), ',', rep(1:K, each = p), "]")]
    beta       = matrix(beta, nrow = p, ncol = K)
    lambda     = pars[paste0("lambdaMat[", rep(1:J, K), ',', rep(1:K, each = J), "]")]
    lambda     = matrix(lambda, nrow = J, ncol = K)
    prior_lp   = sum( sapply(1:K, function(k){
      sum(dnorm(beta[, k], mean = data$beta_mean, sd = data$beta_sd, log = T)) +
        sum(dnorm(lambda[, k], mean = data$hazard_mean, sd = data$hazard_sd, log = T)) - data$lognc_hazard
    }) )

    Eta0       = data$X0 %*% beta
    stratumID0 = data$stratumID0
    a0s        = data$a0s
    y0         = data$y0
    data_lp    = sapply(1:data$n0, function(i){
      a0s[ stratumID0[i] ] * pwe_lpdf(y0[i], Eta0[i, stratumID0[i]], lambda[, stratumID0[i]],
                                      data$breaks, data$intindx0[i], data$J, data$death_ind0[i])
    })
    data_lp    = sum(data_lp)

    if( !data$is_prior ){
      Eta       = data$X1 %*% beta
      stratumID = data$stratumID
      y1        = data$y1
      data_lp = data_lp + sum( sapply(1:data$n1, function(i){
        pwe_lpdf(y1[i], Eta[i, stratumID[i]], lambda[, stratumID[i]],
                 data$breaks, data$intindx[i], data$J, data$death_ind[i])
      }))
    }
    return(data_lp + prior_lp)
  }

  lb           = c(rep(-Inf, p*K), rep(0, J*K))
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
