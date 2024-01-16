#'
#' Second step for sampling from the posterior distribution of a GLM using the Robust Meta-Analytic Predictive (MAP)
#' Prior by Schmidli et al. See `glm_rmap_bhm.R` for the first step.
#'
#' This function approximates the distribution of the output "beta_pred" from the `glm.rmap.bhm()` function (i.e.,
#' samples from the prior induced by the Bayesian Hierarchical Model (BHM)) by a mixture of multivariate normal
#' distributions. We use the `mclust` package by Scrucca et al. to implement the mixture approximation.
#'
#' @export
#'
#' @param samples.bhm       a matrix of the samples of regression coefficients from the prior induced by the BHM, output
#'                          from the `glm.rmap.bhm()` function.
#' @param G                 an integer vector giving the numbers of mixtures components (in the mixture approximation to
#'                          the BHM) for which the BIC is to be calculated. Defaults to 1:9. See the argument `G` in
#'                          \code{\link[mclust:Mclust]{?mclust::Mclust}}. 
#' @param ...               arguments passed to [mclust::Mclust()]. See <https://mclust-org.github.io/mclust/reference/Mclust.html>.
#'
#' @return                  a list giving the parameters estimated from the optimal (according to BIC) mixture model
#'                          including the mixing proportions, the estimated mean and covariance matrix for each mixture
#'                          component. An object of class 'Mclust' from the mclust::Mclust() function is also included.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   data(actg036) ## historical data
#'   ## take subset for speed purposes
#'   actg036 = actg036[1:100, ]
#'   hist_data_list = list(actg036)
#'   samples_bhm = glm.rmap.bhm(
#'     formula = cd4 ~ treatment + age + race,
#'     family = poisson('log'),
#'     hist.data.list = hist_data_list,
#'     chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'   )$beta_pred
#'   glm.rmap.bhm.approx(
#'     samples.bhm = samples_bhm,
#'     G = 1:5, verbose = FALSE
#'   )
#' }
glm.rmap.bhm.approx = function(
  samples.bhm,
  G = NULL,
  ...
) {
  ## mixture approximation using the mclust package
  bhm_approx = mclust::Mclust(samples.bhm, G = G, ...)
  G          = bhm_approx$G
  params     = bhm_approx$parameters
  probs      = params$pro
  means      = params$mean
  covs       = params$variance$sigma
  p          = ncol(samples.bhm)

  ## ensure that the mean matrix has dimension = c(p, G)
  means = matrix(means, nrow = p, ncol = G)
  ## reshape the covariance array to have dimension = c(G, p, p)
  if ( (length(covs) == 1) & (G > 1) ) {
    covs = rep(covs, G)
  }
  if ( is.null(dim(covs)) ){
    covs = array(covs, dim = c(p, p, G))
  }
  covs = aperm(covs, perm = c(3, 1, 2))

  res = list(
    probs        = probs
    , means      = means
    , covs       = covs
    , bhm_approx = bhm_approx
  )
  return(res)
}
