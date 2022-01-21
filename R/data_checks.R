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
  formula, family, data, histdata, offset, offset0, check.hist = TRUE
) {
  if ( !(is.data.frame(data) ) )
    stop('data must be a data.frame')
  if ( any( is.na(data) ) )
    stop("data cannot contain missing values")
  if ( ! ( class(formula) == 'formula' ) )
    stop('formula must be of type "formula"')
  if ( class(family) != 'family')
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  yname = formula.tools::lhs.vars(formula)
  if (length(yname) != 1)
    stop('formula must contain exactly 1 lhs variable name')
  varnames = all.vars(formula)
  if ( !( all( varnames %in% names(data) ) ) )
    stop("formula contains terms not in data")
  if (!is.null(offset))
    if(length(offset) != nrow(data))
      stop('offset and data must have equal lengths if offset is not NULL')
  if(check.hist) {
    if ( !(is.data.frame(histdata) ) )
      stop('histdata must be a data.frame')
    if ( any(is.na(histdata)) )
      stop("histdata cannot contain missing values")
    if ( !( all( varnames %in% names(histdata) ) ) )
      stop("formula contains terms not in histdata")
    if (!is.null(offset0))
      if(length(offset0) != nrow(histdata))
        stop('offset0 and histdata must have equal lengths if offset0 is not NULL')
  }
}
