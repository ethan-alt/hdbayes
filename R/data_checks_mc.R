#' Data validation functions for methods with missing covariates
#'
#' @name data_checks_mc
#' @keywords internal
NULL


#' Check if the input data are in appropriate forms for methods with missing covariates
#'
#' This is a modified version of data.checks that ALLOWS missing values in covariates
#' (but not in the response variable or offsets).
#'
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets.
#' @param offset.list a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                    the length of data.list. The length of each element of offset.list is equal to the number
#'                    of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @noRd
data.checks.mc = function(
    formula, family, data.list, offset.list
) {
  if ( !inherits(formula, "formula") )
    stop('formula must be of type "formula"')
  if ( !inherits(family, 'family') )
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  varnames = all.vars(formula)
  yname = varnames[1]
  covariate_names = varnames[-1]

  if ( !( is.list(data.list) ) )
    stop("data.list must be a list of data.frames")

  total_complete = 0

  for( i in seq_len( length(data.list) ) ){
    if ( !( is.data.frame(data.list[[i]]) ) )
      stop("element ", i, " in data.list must be a data.frame")

    # Check that all variables in formula exist
    if ( !( all( varnames %in% names(data.list[[i]]) ) ) )
      stop("formula contains terms not in element ", i, " in data.list")

    # Response variable cannot have missing values
    if ( any( is.na(data.list[[i]][[yname]]) ) )
      stop("response variable '", yname, "' in element ", i, " in data.list cannot contain missing values")

    # Count complete cases for copula estimation
    cov_data = data.list[[i]][, covariate_names, drop = FALSE]
    n_complete = sum(complete.cases(cov_data))
    total_complete = total_complete + n_complete
  }

  # Need at least some complete cases for copula estimation

  if (total_complete < 3) {
    stop("Need at least 3 complete cases (rows with no missing covariates) across all datasets to estimate copula correlation")
  }

  # Check offsets (offsets cannot have missing values)
  if ( !( is.null(offset.list) ) ){
    if ( !( is.list(offset.list) ) )
      stop("offset.list must be a list of vectors if offset.list is not NULL")
    if ( length(offset.list) != length(data.list) )
      stop("offset.list and data.list must have equal lengths if offset.list is not NULL")
    for( i in seq_len( length(offset.list) ) ){
      if ( !( is.vector(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list must be a vector")
      if ( any( is.na(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list cannot contain missing values")
      if ( length(offset.list[[i]]) != nrow(data.list[[i]]) )
        stop("the length of element ", i, " in offset.list must be equal to the number of rows in element ", i, " in data.list if offset.list is not NULL")
    }
  }
}


#' Check if any missing covariates exist in data.list
#'
#' @param data.list a list of data.frames
#' @param formula a two-sided formula
#' @return logical indicating if there are missing covariate values
#' @noRd
has.missing.covariates = function(data.list, formula) {
  covariate_names = all.vars(formula)[-1]
  for (d in data.list) {
    cov_data = d[, covariate_names, drop = FALSE]
    if (any(is.na(cov_data))) {
      return(TRUE)
    }
  }
  return(FALSE)
}


#' Check if the input data are in appropriate forms for LEAP with missing covariates
#'
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets.
#' @param K           the desired number of classes to identify for LEAP implementation.
#' @param offset.list a list of matrices giving the offset for current data followed by historical data. For each
#'                    matrix, the number of rows corresponds to observations and columns correspond to classes.
#' @noRd
data.checks.mc.leap = function(
    formula, family, data.list, K, offset.list
) {
  data.checks.mc(formula, family, data.list, NULL)

  if ( K < 2 )
    stop("K should be at least 2")

  if ( K != 2 ){
    if( !is.numeric(K) )
      stop("K must be a numeric value")
  }

  if( !( is.null(offset.list) ) ){
    if ( !( is.list(offset.list) ) )
      stop("offset.list must be a list of matrices if offset.list is not NULL")
    if ( length(offset.list) != length(data.list) )
      stop("offset.list and data.list must have equal lengths if offset.list is not NULL")
    for( i in seq_len( length(offset.list) ) ){
      if ( !( is.matrix(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list must be a matrix")
      if ( nrow(offset.list[[i]]) != nrow(data.list[[i]]) )
        stop("element ", i, " in offset.list must have the same number of rows as element ", i, " in data.list if offset.list is not NULL")
      if ( ncol(offset.list[[i]]) != K )
        stop("element ", i, " in offset.list must have the same number of columns as K if offset.list is not NULL")
      if ( any( is.na(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list cannot contain missing values")
    }
  }
}
