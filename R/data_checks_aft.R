#' check if the input data are in appropriate forms for all methods
#' @param formula           a two-sided formula giving the relationship between the response variable and covariates.
#'                          The response is a survival object as returned by the `survival::Surv(time, event)` function,
#'                          where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type
#'                          of censoring is assumed to be right-censoring.
#' @param data.list         a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                          are the historical data sets. For fitting accelerated failure time (AFT) models, all
#'                          historical data sets will be stacked into one historical data set.
#' @param dist              a character indicating the distribution of survival times. Currently, `dist` can be one
#'                          of the following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".
#' @param strata.list       a list of vectors specifying the stratum ID for each observation in the corresponding data set
#'                          in `data.list`. The first element in the list corresponds to the current data, and the rest
#'                          correspond to the historical data sets. Each vector should have the same length as the number
#'                          of rows in the respective data set in `data.list`, with values representing stratum labels
#'                          as positive integers (e.g., 1, 2, 3, ...). Defaults to NULL.
#' @param is.stratified.pp  whether the method is the stratified power prior. Defaults to FALSE.
#'
#' @noRd
data.checks.aft = function(
    formula, data.list, dist = "weibull", strata.list = NULL, is.stratified.pp = FALSE
) {
  if ( !inherits(formula, "formula") )
    stop('formula must be of type "formula"')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  varnames = all.vars(formula)
  if ( !( is.list(data.list) ) )
    stop("data.list must be a list of data.frames")
  for( i in seq_len( length(data.list) ) ){
    if ( !( is.data.frame(data.list[[i]]) ) )
      stop("element ", i, " in data.list must be a data.frame")
    if ( any( is.na(data.list[[i]]) ) )
      stop("element ", i, " in data.list cannot contain missing values")
    if ( !( all( varnames %in% names(data.list[[i]]) ) ) )
      stop("formula contains terms not in element ", i, " in data.list")
  }
  if( !is.character(dist) ){
    stop('dist must be of type "character"')
  }else if( ! tolower(dist) %in% c("lognormal", "loglogistic", "weibull") ){
    stop('Distribution not supported. dist must be one of the following values: "lognormal", "loglogistic", or "weibull"')
  }

  if( is.stratified.pp ){
    if ( !( is.list(strata.list) ) )
      stop("strata.list must be a list of vectors giving the stratum ID for each observation")
    if ( length(strata.list) != length(data.list) )
      stop("strata.list and data.list must have equal lengths")
    for( i in seq_len( length(strata.list) ) ){
      if ( !( is.vector(strata.list[[i]]) ) )
        stop("element ", i, " in strata.list must be a vector")
      if ( any( is.na(strata.list[[i]]) ) )
        stop("element ", i, " in strata.list cannot contain missing values")
      if ( length(strata.list[[i]]) != nrow(data.list[[i]]) )
        stop("the length of element ", i, " in strata.list must be equal to the number of rows in element ", i, " in data.list")
    }
  }
}


#' turn character value of dist into an integer in 1, 2, 3
#'
#' @param dist        a character indicating the distribution of survival times. Currently, `dist` can be one
#'                    of the following values: "weibull", "lognormal", or "loglogistic". Defaults to "weibull".
#'
#' @noRd
dist.to.integer = function(dist){
  dist_map = c("lognormal" = 1, "loglogistic" = 2, "weibull" = 3)
  res      = dist_map[tolower(dist)]
  return(as.integer(res))
}


#' transfer a scalar/vector/NULL into a vector of given length and (default) values
#' @param param         a scalar or a vector if param is not NULL
#' @param default.value default value for param. Defaults to 0.
#' @param len           length (number of elements) of param. Defaults to 1.
#' @noRd
to.vector = function(
    param, default.value = 0, len = 1
) {
  if ( is.null(param) ){
    param = rep(default.value, len)
  }else if ( length(param) == 1 ){
    param = rep(as.numeric(param), len)
  }else {
    param = as.numeric(param)
  }
  return(param)
}


#' change the variable names of the `draws_df` object obtained from the input `CmdStanMCMC` object and
#' reorder the variables so that the updated variable names appear at the top of the `draws_df` object.
#' @param fit      an object of class `CmdStanMCMC`.
#' @param oldnames a vector of `character` giving the parameter/variable names in fit to be changed.
#' @param newnames a vector of `character` giving the new parameter/variable names.
#' @noRd
rename.params = function(
    fit, oldnames, newnames
) {
  pars = fit$metadata()$model_params
  pars = c(pars[1], oldnames, (pars[!pars %in% oldnames])[-1])
  d    = fit$draws(format = 'draws_df', variables = pars)
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}
