#' @keywords internal
"_PACKAGE"
#'
#' The 'hdbayes' package.
#'
#' @description Bayesian analysis of generalized linear models using historical data
#'
#' @name hdbayes-package
#' @aliases hdbayes
#' @family help
#' @importFrom instantiate stan_package_model
#' @importFrom callr r
#' @importFrom fs dir_copy
#' @importFrom formula.tools is.two.sided lhs.vars
#' @importFrom stats model.matrix family glm binomial dgamma dnorm pnorm lm gaussian dbeta pbeta rbinom dlogis plogis
#' @importFrom posterior variables merge_chains
#' @importFrom enrichwith enrich
#' @importFrom bridgesampling bridge_sampler
#' @importFrom mvtnorm dmvnorm
NULL
