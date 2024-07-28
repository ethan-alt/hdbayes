#' Posterior of robust meta-analytic predictive prior (RMAP)
#'
#' Sample from the posterior distribution of a GLM using the robust meta-analytic predictive prior (RMAP)
#' by Schmidli et al. (2014) <doi:10.1111/biom.12242>.
#'
#' The robust meta-analytic predictive prior (RMAP) is a two-part mixture prior consisting of a meta-analytic
#' predictive (MAP) prior (the prior induced by Bayesian hierarchical model (BHM)) and a vague (i.e.,
#' non-informative) prior (specifically, the normal/half-normal prior with large variances). Although Schmidli et al.
#' (2014) recommends to use a finite mixture of conjugate priors to approximate the BHM, it can be difficult and
#' time-consuming to come up with an appropriate approximation.
#'
#' Instead, the approach taken by hdbayes is to use the marginal likelihood of the MAP and vague priors.
#' Specifically, note that the posterior distribution of a GLM under RMAP is also a two-part mixture distribution.
#' The updated mixture weight for posterior density under the MAP prior is
#' \deqn{\widetilde{w} = \frac{w Z_I(D, D_0)}{w Z_I(D, D_0) + (1-w) Z_V(D)},}
#' where \eqn{w} is the prior mixture weight for the MAP prior in RMAP, \eqn{Z_I(D, D_0)} is the marginal likelihood
#' of the MAP prior, and \eqn{Z_V(D)} is the marginal likelihood of the vague prior.
#'
#' @include glm_bhm.R
#' @include glm_logml_map.R
#' @include glm_post.R
#' @include glm_logml_post.R
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets.
#' @param offset.list       a list of vectors giving the offsets for each data. The length of `offset.list` is equal to
#'                          the length of `data.list`. The length of each element of `offset.list` is equal to the number
#'                          of rows in the corresponding element of `data.list`. Defaults to a list of vectors of 0s.
#' @param w                 a scalar between 0 and 1 giving how much weight to put on the historical data. Defaults to 0.1.
#' @param meta.mean.mean    same as `meta.mean.mean` in [glm.bhm()]. It is a scalar or a vector whose dimension is equal
#'                          to the number of regression coefficients giving the means for the normal hyperpriors on the
#'                          mean hyperparameters of regression coefficients in Bayesian hierarchical model (BHM). If a
#'                          scalar is provided, `meta.mean.mean` will be a vector of repeated elements of the given scalar.
#'                          Defaults to a vector of 0s.
#' @param meta.mean.sd      same as `meta.mean.sd` in [glm.bhm()]. It is a scalar or a vector whose dimension is equal
#'                          to the number of regression coefficients giving the sds for the normal hyperpriors on the
#'                          mean hyperparameters of regression coefficients in BHM. If a scalar is provided, same as for
#'                          `meta.mean.mean`. Defaults to a vector of 10s.
#' @param meta.sd.mean      same as `meta.sd.mean` in [glm.bhm()]. It is a scalar or a vector whose dimension is equal
#'                          to the number of regression coefficients giving the means for the half-normal hyperpriors
#'                          on the sd hyperparameters of regression coefficients in BHM. If a scalar is provided, same
#'                          as for `meta.mean.mean`. Defaults to a vector of 0s.
#' @param meta.sd.sd        same as `meta.sd.sd` in [glm.bhm()]. It is a scalar or a vector whose dimension is equal to
#'                          the number of regression coefficients giving the sds for the half-normal hyperpriors on the
#'                          sd hyperparameters of regression coefficients in BHM. If a scalar is provided, same as for
#'                          `meta.mean.mean`. Defaults to a vector of 1s.
#' @param disp.mean         a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the location parameters for the half-normal priors on the dispersion parameters. If
#'                          a scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 0s.
#' @param disp.sd           a scalar or a vector whose dimension is equal to the number of data sets (including the current
#'                          data) giving the scale parameters for the half-normal priors on the dispersion parameters. If a
#'                          scalar is provided, same as for `meta.mean.mean`. Defaults to a vector of 10s.
#' @param norm.vague.mean   same as `beta.mean` in [glm.post()]. It is a scalar or a vector whose dimension is equal to the
#'                          number of regression coefficients giving the means for the vague normal prior on regression
#'                          coefficients. If a scalar is provided, `norm.vague.mean` will be a vector of repeated elements
#'                          of the given scalar. Defaults to a vector of 0s.
#' @param norm.vague.sd     same as `beta.sd` in [glm.post()]. It is a scalar or a vector whose dimension is equal to the
#'                          number of regression coefficients giving the sds for the vague normal prior on regression
#'                          coefficients. If a scalar is provided, same as for `norm.vague.mean`. Defaults to a vector of 10s.
#' @param bridge.args       a `list` giving arguments (other than samples, log_posterior, data, lb, ub) to pass
#'                          onto [bridgesampling::bridge_sampler()].
#' @param iter_warmup       number of warmup iterations to run per chain. Defaults to 1000. See the argument `iter_warmup` in
#'                          `sample()` method in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Defaults to 1000. See the argument `iter_sampling`
#'                          in `sample()` method in cmdstanr package.
#' @param chains            number of Markov chains to run. Defaults to 4. See the argument `chains` in `sample()` method in
#'                          cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g. seed, refresh, init).
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{post.samples}{an object of class `draws_df` giving posterior samples under the robust meta-analytic predictive prior (RMAP)}
#'
#'    \item{post.samples.bhm}{an object of class `draws_df` giving posterior samples under the Bayesian hierarchical model (BHM),
#'    obtained from using [glm.bhm()]}
#'
#'    \item{post.samples.vague}{an object of class `draws_df` giving posterior samples under the vague/non-informative prior, obtained
#'    from using [glm.post()]}
#'
#'    \item{bs.map}{output from computing log marginal likelihood of the prior induced by the BHM (referred to as the meta-analytic predictive
#'    (MAP) prior) via [glm.logml.map()] function}
#'
#'    \item{bs.vague}{output from computing log marginal likelihood of the vague prior via [glm.logml.post()] function}
#'  }
#'
#' @references
#'  Schmidli, H., Gsteiger, S., Roychoudhury, S., O’Hagan, A., Spiegelhalter, D., and Neuenschwander, B. (2014). Robust meta‐analytic‐predictive priors in clinical trials with historical control information. Biometrics, 70(4), 1023–1032.
#'
#'  Gronau, Q. F., Singmann, H., and Wagenmakers, E.-J. (2020). bridgesampling: An r package for estimating normalizing constants. Journal of Statistical Software, 92(10).
#'
#' @examples
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data(actg019) ## current data
#'     data(actg036) ## historical data
#'     ## take subset for speed purposes
#'     actg019 = actg019[1:150, ]
#'     actg036 = actg036[1:100, ]
#'     data.list = list(actg019, actg036)
#'     glm.rmap(
#'       formula = outcome ~ scale(age) + race + treatment + scale(cd4),
#'       family = binomial('logit'),
#'       data.list = data.list,
#'       w = 0.1,
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'   }
#' }
glm.rmap = function(
    formula,
    family,
    data.list,
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
  d.bhm = glm.bhm(
    formula        = formula,
    family         = family,
    data.list      = data.list,
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
  d.vague = glm.post(
    formula        = formula,
    family         = family,
    data.list      = data.list,
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
  ## updated mixture weight = w * Z_I / (w * Z_I + (1-w) * Z_V ), where
  ## Z_I and Z_V are marginal likelihoods of the MAP and vague priors, respectively
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
    'post.samples.bhm'   = d.bhm,
    'post.samples.vague' = d.vague,
    'bs.map'             = bs.map,
    'bs.vague'           = bs.vague
  )
  return(res)
}
