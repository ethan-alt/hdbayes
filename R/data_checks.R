#' get distribution and link for a GLM
#' @exportS3Method NULL
#' @param family   an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
get.dist.link = function(family) {
  fams  = c('binomial', 'poisson', 'gaussian', 'Gamma', 'inverse.gaussian')
  links = c(
    'identity', 'log', 'logit', 'inverse', 'probit', 'cauchit', 'cloglog'
    ,'sqrt', '1/mu^2'
  )
  fam.id  = which(family$family == fams)
  link.id = which(family$link == links)
  c(fam.id, link.id)
}

#' check if the input data are in appropriate forms
#' @exportS3Method    NULL
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical datasets.
#' @param offset.list a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                    the length of data.list. The length of each element of offset.list is equal to the number
#'                    of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
data.checks = function(
    formula, family, data.list, offset.list
) {
  if ( !inherits(formula, "formula") )
    stop('formula must be of type "formula"')
  if ( !inherits(family, 'family') )
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  yname = formula.tools::lhs.vars(formula)
  if ( length(yname) != 1 )
    stop('formula must contain exactly 1 lhs variable name')
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

#' reshape the input data list into a list which contains the response vector y (from all data sets), the covariate matrix X (from all data sets),
#' and the starting and ending indices of each data set
#' @exportS3Method    NULL
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical datasets.
stack.data = function(
    formula, data.list
) {
  ## get stacked design matrix and response vector using formula
  X = lapply(data.list, function(s){
    stats::model.matrix(formula, s)
  })
  X = do.call(rbind, X)
  y = lapply(data.list, function(s){
    s[, all.vars(formula)[1]]
  })
  y = unlist(y)

  ## reset indices of X
  rownames(X) = NULL

  # get starting and ending indices for each data
  num.obs = sapply(data.list, function(s){
    nrow(s)
  })
  end.index   = cumsum(num.obs)
  start.index = c(1, end.index[-length(data.list)] + 1)

  return(list(X = X, y = y, start.index = start.index, end.index = end.index))
}


#' transfer a scalar/vector/NULL into a vector of given length and (default) values
#' @exportS3Method NULL
#' @param param         a scalar or a vector if param is not NULL
#' @param default.value default value for param. Defaults to 0.
#' @param len           length (number of elements) of param. Defaults to 1.
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
