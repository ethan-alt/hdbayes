#' Posterior of robust meta-analytic predictive prior (RMAP) with missing covariates
#'
#' Sample from the posterior distribution of a GLM using the robust meta-analytic predictive prior (RMAP)
#' by Schmidli et al. (2014) <doi:10.1111/biom.12242>, extended to handle missing covariates via a
#' Gaussian copula model.
#'
#' The robust meta-analytic predictive prior (RMAP) is a two-part mixture prior consisting of a meta-analytic
#' predictive (MAP) prior (the prior induced by Bayesian hierarchical model (BHM)) and a vague (i.e.,
#' non-informative) prior (specifically, the normal/half-normal prior with large variances).
#'
#' This function extends RMAP to handle missing covariate values using a Gaussian copula approach. Missing
#' values are imputed by sampling from the conditional distribution implied by the copula, with the
#' correlation matrix either estimated from complete cases or provided by the user.
#'
#' @include glm_bhm_mc.R
#' @include glm_logml_map.R
#' @include glm_post_mc.R
#' @include glm_logml_post.R
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. Covariates may contain missing values (NA), which will be
#'                          imputed using the Gaussian copula model.
#' @param copula.corr       optional correlation matrix for the Gaussian copula. If NULL (default), the correlation
#'                          is estimated from complete cases across all datasets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of `offset.list` is equal to
#'                          the length of `data.list`. The length of each element of `offset.list` is equal to the number
#'                          of rows in the corresponding element of `data.list`. Defaults to a list of vectors of 0s.
#' @param w                 a scalar between 0 and 1 giving how much weight to put on the historical data. Defaults to 0.1.
#' @param meta.mean.mean    same as `meta.mean.mean` in [glm.bhm.mc()]. It is a scalar or a vector whose dimension is equal
#'                          to the number of regression coefficients giving the means for the normal hyperpriors on the
#'                          mean hyperparameters of regression coefficients in Bayesian hierarchical model (BHM).
#' @param meta.mean.sd      same as `meta.mean.sd` in [glm.bhm.mc()]. It is a scalar or a vector giving the sds for the
#'                          normal hyperpriors on the mean hyperparameters of regression coefficients in BHM.
#' @param meta.sd.mean      same as `meta.sd.mean` in [glm.bhm.mc()]. It is a scalar or a vector giving the means for the
#'                          half-normal hyperpriors on the sd hyperparameters of regression coefficients in BHM.
#' @param meta.sd.sd        same as `meta.sd.sd` in [glm.bhm.mc()]. It is a scalar or a vector giving the sds for the
#'                          half-normal hyperpriors on the sd hyperparameters of regression coefficients in BHM.
#' @param disp.mean         a scalar or a vector giving the location parameters for the half-normal priors on the
#'                          dispersion parameters. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector giving the scale parameters for the half-normal priors on the
#'                          dispersion parameters. Defaults to a vector of 10s.
#' @param norm.vague.mean   same as `beta.mean` in [glm.post.mc()]. It is a scalar or a vector giving the means for
#'                          the vague normal prior on regression coefficients. Defaults to a vector of 0s.
#' @param norm.vague.sd     same as `beta.sd` in [glm.post.mc()]. It is a scalar or a vector giving the sds for the
#'                          vague normal prior on regression coefficients. Defaults to a vector of 10s.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000.
#' @param chains            number of Markov chains to run. Defaults to 4.
#' @param ...               arguments passed to `sample()` method in cmdstanr package.
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{post.samples}{an object of class `draws_df` giving posterior samples under the robust meta-analytic predictive prior (RMAP)}
#'
#'    \item{post.weight.bhm}{a scalar between 0 and 1 giving the updated mixture weight for posterior density under the MAP prior}
#'
#'    \item{post.samples.bhm}{an object of class `draws_df` giving posterior samples under the Bayesian hierarchical model (BHM)}
#'
#'    \item{post.samples.vague}{an object of class `draws_df` giving posterior samples under the vague/non-informative prior}
#'
#'    \item{bs.map}{output from computing log marginal likelihood of the MAP prior}
#'
#'    \item{bs.vague}{output from computing log marginal likelihood of the vague prior}
#'  }
#'
#' @references
#'  Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., and Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in clinical trials with historical control information. Biometrics, 70(4), 1023-1032.
#'
#' @examples
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data(actg019)
#'     data(actg036)
#'     ## take subset for speed purposes
#'     actg019 = actg019[1:150, ]
#'     actg036 = actg036[1:100, ]
#'     ## introduce missing values
#'     actg036$age[sample(1:100, 10)] = NA
#'     data.list = list(actg019, actg036)
#'     glm.rmap.mc(
#'       formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'       family = binomial('logit'),
#'       data.list = data.list,
#'       w = 0.1,
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'   }
#' }
glm.rmap.mc = function(
    formula,
    family,
    data.list,
    copula.corr       = NULL,
    offset.list       = NULL,
    w                 = 0.1,
    meta.mean.mean    = NULL,
    meta.mean.sd      = NULL,
    meta.sd.mean      = NULL,
    meta.sd.sd        = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    norm.vague.mean   = NULL,
    norm.vague.sd     = NULL,
    bridge.args       = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get posterior samples under BHM (equivalently, under the MAP prior)
  d.bhm = glm.bhm.mc(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    copula.corr    = copula.corr,
    offset.list    = offset.list,
    meta.mean.mean = meta.mean.mean,
    meta.mean.sd   = meta.mean.sd,
    meta.sd.mean   = meta.sd.mean,
    meta.sd.sd     = meta.sd.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    iter_warmup    = iter_warmup,
    iter_sampling  = iter_sampling,
    chains         = chains,
    ...
  )
  ## compute log marginal likelihood of the MAP prior
  bs.map = glm.logml.map(
    post.samples   = d.bhm,
    bridge.args    = bridge.args,
    iter_warmup    = iter_warmup,
    iter_sampling  = iter_sampling,
    chains         = chains,
    ...
  )
  logml.map = bs.map$logml

  ## get posterior samples under the vague prior
  d.vague = glm.post.mc(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    copula.corr    = copula.corr,
    offset.list    = offset.list,
    beta.mean      = norm.vague.mean,
    beta.sd        = norm.vague.sd,
    disp.mean      = disp.mean,
    disp.sd        = disp.sd,
    iter_warmup    = iter_warmup,
    iter_sampling  = iter_sampling,
    chains         = chains,
    ...
  )
  ## compute log marginal likelihood under the vague prior
  bs.vague = glm.logml.post(
    post.samples   = d.vague,
    bridge.args    = bridge.args
  )
  logml.vague = bs.vague$logml

  ## compute updated mixture weight for posterior density under the MAP prior using marginal likelihoods and w
  post.wt.log = c( log(w) + logml.map, log(1-w) + logml.vague )
  post.wt     = exp( post.wt.log[1] - log_sum_exp( post.wt.log ) )

  m        = nrow(d.vague) ## number of posterior samples under RMAP
  varnames = colnames(d.vague)
  d.bhm2   = d.bhm[, varnames]
  ## merge chains of d.bhm into a single chain
  d.bhm2   = posterior::merge_chains(d.bhm2)
  ## merge chains of d.vague into a single chain
  d.vague2 = posterior::merge_chains(d.vague)

  ## obtain posterior samples under RMAP
  d        = lapply(1:m, function(i){
    x = rbinom(1, 1, prob = post.wt)
    if( x == 1 ){
      d.bhm2[sample(1:m, 1), ]
    }else{
      d.vague2[sample(1:m, 1), ]
    }
  })
  d        = do.call(rbind, d)

  res      = list(
    'post.samples'       = d,
    'post.weight.bhm'    = post.wt,
    'post.samples.bhm'   = d.bhm,
    'post.samples.vague' = d.vague,
    'bs.map'             = bs.map,
    'bs.vague'           = bs.vague
  )
  return(res)
}
