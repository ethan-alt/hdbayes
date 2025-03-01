#' check if the input data are in appropriate forms for all methods
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#'                    The response is a survival object as returned by the `survival::Surv(time, event)` function,
#'                    where event is a binary indicator for event (0 = no event, 1 = event has occurred). The type
#'                    of censoring is assumed to be right-censoring.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets. For fitting piecewise exponential (PWE) models, all
#'                    historical data sets will be stacked into one historical data set.
#' @param breaks      a numeric vector specifying the time points that define the boundaries of the piecewise
#'                    intervals. The values should be in ascending order, with the final value being greater
#'                    than or equal to the maximum observed time.
#'
#' @noRd
data.checks.pwe = function(
    formula, data.list, breaks
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
  if( !is.numeric(breaks) ){
    stop('breaks must be a vector of numeric values')
    if( any( order(breaks) != 1:length(breaks) ) ){
      breaks = breaks[order(breaks)] # breaks should be in ascending order
    }
  }
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
