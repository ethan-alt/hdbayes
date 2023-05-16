library(formula.tools)

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


data.checks = function(
    formula, family, data.list, offset.list
) {
  if ( !( class(formula) == 'formula' ) )
    stop('formula must be of type "formula"')
  if ( class(family) != 'family' )
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  yname = formula.tools::lhs.vars(formula)
  if ( length(yname) != 1 )
    stop('formula must contain exactly 1 lhs variable name')
  varnames = all.vars(formula)
  if ( !( is.list(data.list) ) )
    stop("data.list must be a list of data.frames")
  message("the first element in data.list is regarded as the current data")
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


stack.data = function(
    formula, data.list, include.intercept = TRUE
) {
  ## get stacked design matrix and response vector using formula
  X = sapply(data.list, function(s){
    model.matrix(formula, s)
  })
  X = do.call(rbind, X)
  if ( !include.intercept ) {
    X = X[, -1, drop = F]
  }
  y = sapply(data.list, function(s){
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
