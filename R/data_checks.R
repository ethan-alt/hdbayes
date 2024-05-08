#' get distribution and link for a GLM
#' @param family   an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @noRd
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

#' check if the input data are in appropriate forms for all methods except for LEAP
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets.
#' @param offset.list a list of vectors giving the offsets for each data. The length of offset.list is equal to
#'                    the length of data.list. The length of each element of offset.list is equal to the number
#'                    of rows in the corresponding element of data.list. Defaults to a list of vectors of 0s.
#' @noRd
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

#' check if the input data are in appropriate forms for LEAP
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param family      an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets.
#' @param K           the desired number of classes to identify for LEAP implementation.
#' @param offset.list a list of matrices giving the offset for current data followed by historical data. For each
#'                    matrix, the number of rows corresponds to observations and columns correspond to classes.
#' @noRd
data.checks.leap = function(
    formula, family, data.list, K, offset.list
) {
  data.checks(formula, family, data.list, NULL)

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

#' reshape the input data list into a list which contains the response vector y (from all data sets), the covariate matrix X (from all data sets),
#' and the starting and ending indices of each data set
#' @param formula     a two-sided formula giving the relationship between the response variable and covariates.
#' @param data.list   a list of `data.frame`s. The first element in the list is the current data, and the rest
#'                    are the historical data sets.
#' @noRd
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


#' check if the input `post.samples` are in appropriate forms for computing log marginal likelihood under different priors.
#' @param post.samples      an object of class `draws_df`, `draws_matrix`, `matrix`, or `data.frame` giving posterior
#'                          samples of a GLM under different priors. Each row corresponds to the posterior samples obtained
#'                          from one iteration of MCMC. The column names of `post.samples` should include the names of
#'                          covariates for regression coefficients, such as "(Intercept)", and "dispersion" for the
#'                          dispersion parameter, if applicable.
#' @param covariate.names   a vector of `character` giving the names of the covariates.
#' @param family            an object of class `family`. See \code{\link[stats:family]{?stats::family}}.
#' @noRd
post.samples.checks = function(
    post.samples, covariate.names, family
) {
  if ( !any( class(post.samples) %in% c("draws_df", "draws_matrix", "matrix", "data.frame") ) )
    stop("post.samples must be in one of the following formats: draws_df, draws_matrix, matrix, or data.frame")
  if ( !all( covariate.names %in% colnames(post.samples) ) )
    stop("Column names of post.samples must include the names of covariates for regression coefficients, such as \"(Intercept)\"")
  if ( !family$family %in% c('binomial', 'poisson') ) {
    if ( !("dispersion" %in% colnames(post.samples)) )
      stop("Column names of post.samples must include \"dispersion\" for dispersion parameter")
  }
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
