#' Compute model averaging weights
#'
#' Compute model averaging weights for a set of Bayesian models using Bayesian model averaging (BMA), pseudo-BMA,
#' pseudo-BMA+ (pseudo-BMA with the Bayesian bootstrap), or stacking. This function takes a list of model fit objects,
#' each containing posterior samples from a generalized linear model (GLM) or survival model, and returns normalized
#' weights that can be used for model comparison or combining posterior samples using functions like [sample.ensemble()].
#'
#' The input `fit.list` should be a list of outputs from model fitting functions in the `hdbayes` package, such as [glm.pp()]
#' (for generalized linear models), [aft.pp()] (for accelerated failure time models), [pwe.pp()] (for piecewise exponential (PWE)
#' models), or [curepwe.pp()] (for mixture cure rate models with a PWE component for the non-cured population). To compute
#' pseudo-BMA, pseudo-BMA+, or stacking weights, each fit must include pointwise log-likelihood values. To ensure this, the
#' fitting function must be called with `get.loglik = TRUE`.
#'
#' The arguments related to Markov chain Monte Carlo (MCMC) sampling are utilized to compute the logarithm of the normalizing
#' constant for BMA, if applicable.
#'
#' @include mixture_loglik.R
#'
#' @export
#'
#' @param fit.list          a list of model fit objects returned by functions in the `hdbayes` package. Each fit contains
#'                          posterior samples from a generalized linear model (GLM) (e.g., via [glm.pp()]), an accelerated
#'                          failure time (AFT) model (e.g., via [aft.pp()]), a piecewise exponential (PWE) model (e.g., via
#'                          [pwe.pp()]), or a mixture cure rate model with a PWE component for the non-cured population (e.g.,
#'                          via [curepwe.pp()]). Each fit also includes two attributes: `data`, a list of variables specified
#'                          in the data block of the Stan program, and `model`, a character string indicating the model name.
#'                          To compute pseudo-BMA, pseudo-BMA+, or stacking weights, the fitting function must be called with
#'                          `get.loglik = TRUE`.
#' @param type              a character string specifying the ensemble method used to compute model weights. Options are "bma"
#'                          (Bayesian model averaging (BMA)), "pseudobma" (pseudo-BMA without the Bayesian bootstrap), "pseudobma+"
#'                          (pseudo-BMA with the Bayesian bootstrap), and "stacking".
#' @param prior.prob        a numeric vector of prior model probabilities, used only when `type = "bma"`. Must be non-negative
#'                          and sum to 1. If set to NULL, a uniform prior is used (i.e., all models are equally likely). Defaults
#'                          to NULL.
#' @param bridge.args       a `list` of optional arguments (excluding `samples`, `log_posterior`, `data`, `lb`, and `ub`) to be
#'                          passed to [bridgesampling::bridge_sampler()]. These arguments are used when estimating the log marginal
#'                          likelihood, which is required if `type = "bma"`.
#' @param loo.args          a `list` of optional arguments (excluding `x`) to be passed to [loo::loo()] when computing pseudo-BMA,
#'                          pseudo-BMA+, or stacking weights.
#' @param loo.wts.args      a `list` of optional arguments (excluding `x`, `method`, and `BB`) to be passed to [loo::loo_model_weights()]
#'                          when computing pseudo-BMA, pseudo-BMA+, or stacking weights.
#' @param iter_warmup       number of warmup iterations to run per chain. Used only when computing the log marginal likelihood
#'                          (i.e., when `type = "bma"`). Defaults to 1000. See the argument `iter_warmup` in `sample()` method
#'                          in cmdstanr package.
#' @param iter_sampling     number of post-warmup iterations to run per chain. Used only when computing the log marginal likelihood
#'                          (i.e., when `type = "bma"`). Defaults to 1000. See the argument `iter_sampling` in `sample()` method
#'                          in cmdstanr package.
#' @param chains            number of Markov chains to run. Used only when computing the log marginal likelihood (i.e., when
#'                          `type = "bma"`). Defaults to 4. See the argument `chains` in `sample()` method in cmdstanr package.
#' @param ...               arguments passed to `sample()` method in cmdstanr package (e.g., `seed`, `refresh`, `init`). These are
#'                          used only when computing the log marginal likelihood (i.e., when `type = "bma"`).
#'
#' @return
#'  The function returns a `list` with the following objects
#'
#'  \describe{
#'    \item{weights}{a numeric vector of normalized model weights corresponding to the models in `fit.list`. The names of the
#'    weights are made unique based on the model identifiers.}
#'
#'    \item{type}{a character string indicating the method used to compute the model weights (e.g., "bma", "pseudobma", "pseudobma+",
#'    or "stacking")}
#'
#'    \item{res.logml}{a list of log marginal likelihood estimation results, returned only when `type = "bma"`}
#'
#'    \item{loo.list}{a list of outputs from [loo::loo()], returned only when `type` is "pseudobma", "pseudobma+", or "stacking"}
#'  }
#'
#' @seealso [sample.ensemble()]
#'
#' @references
#'  Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). Using stacking to average Bayesian predictive distributions. Bayesian Analysis, 13(3), 917–1007.
#'
#'  Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413–1432.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   if(requireNamespace("survival")){
#'     library(survival)
#'     data(E1684)
#'     data(E1690)
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
#'     fit.pwe.pp = pwe.pp(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       a0 = 0.5,
#'       get.loglik = TRUE,
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'     fit.pwe.post = pwe.post(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       breaks = breaks,
#'       get.loglik = TRUE,
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'     fit.aft.post = aft.post(
#'       formula = survival::Surv(failtime, failcens) ~ treatment + sex + cage + node_bin,
#'       data.list = data_list,
#'       dist = "weibull",
#'       beta.sd = 10,
#'       get.loglik = TRUE,
#'       chains = 1, iter_warmup = 1000, iter_sampling = 2000
#'     )
#'     compute.ensemble.weights(
#'       fit.list = list(fit.pwe.post, fit.pwe.pp, fit.aft.post),
#'       type = "pseudobma+",
#'       loo.args = list(save_psis = FALSE),
#'       loo.wts.args = list(optim_method="BFGS")
#'     )
#'   }
#' }
compute.ensemble.weights = function(
    fit.list,
    type              = c("bma", "pseudobma", "pseudobma+", "stacking"),
    prior.prob        = NULL,
    bridge.args       = NULL,
    loo.args          = NULL,
    loo.wts.args      = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  type   = match.arg(type)
  models = sapply(fit.list, function(d){
    attr(d, "model")
  })
  if (any(models == "lm_npp")) {
    stop("The model lm_npp is not supported for this function.")
  }
  K      = length(models)
  dots   = list(...)

  if( type == "bma" ){
    res.logml = lapply(seq_along(fit.list), function(i){
      model_name = models[i]
      parts      = strsplit(model_name, "_")[[1]]
      prefix     = parts[1]
      suffix     = ifelse(parts[2] == "bhm", "map", parts[2])
      func_name  = paste0(prefix, ".logml.", suffix)

      if( !exists(func_name, mode = "function") ){
        stop(paste("Function", func_name, "not found"))
      }

      func      = get(func_name, mode = "function")
      arg_names = names(formals(func))
      func_args = list(post.samples = fit.list[[i]], bridge.args = bridge.args)

      if( "iter_warmup" %in% arg_names ){
        func_args$iter_warmup = iter_warmup
      }
      if( "iter_sampling" %in% arg_names ){
        func_args$iter_sampling = iter_sampling
      }
      if( "chains" %in% arg_names ){
        func_args$chains <- chains
      }
      if("..." %in% arg_names ){
        func_args <- c(func_args, dots)
      }
      ## call the function
      do.call(func, func_args)
    })
    logml = vapply(res.logml, function(r){r$logml}, numeric(1))

    if( is.null(prior.prob) ){
      prior.prob = rep(1/K, K)
    }
    if( length(prior.prob) != K ){
      stop("fit.list and prior.prob must have the same length.")
    }
    if( any(prior.prob < 0) || abs(sum(prior.prob) - 1) > 10^(-6) ){
      stop("prior.prob must be non-negative and sum to 1.")
    }

    ## compute BMA weights
    log_wts_unnorm = log(prior.prob) + logml
    log_denom      = log_sum_exp(log_wts_unnorm)
    wts            = exp(log_wts_unnorm - log_denom)

  }else{
    loo.list = lapply(fit.list, function(f){
      loglik = suppressWarnings({
        as.matrix( f[ , grepl("^log_lik", names(f))] )
      })
      do.call(loo::loo, c(list(loglik), loo.args))
    })

    loo.wts.args = c(list(x = loo.list, method = type), loo.wts.args)
    if( type == "pseudobma+" ){
      loo.wts.args$method = "pseudobma"
      loo.wts.args$BB = TRUE
    }else if( type == "pseudobma" ){
      loo.wts.args$BB = FALSE
    }
    # compute pseudo-BMA or pseudo-BMA+ or Stacking weights
    wts = as.numeric(do.call(loo::loo_model_weights, loo.wts.args))
  }

  unique.models = make.unique(models, sep = "_")
  names(wts)    = unique.models
  res           = list(
    weights     = wts,
    type        = type
  )
  if( type == "bma" ){
    res$res.logml = res.logml
  }else{
    res$loo.list = loo.list
  }
  return(res)
}
