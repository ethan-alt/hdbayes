#' International Breast Cancer Study Group (IBCSG) Trial VI Data
#' from patients under 40 (i.e., age < 40)
#' We treat it as the historical data
#'
#' @name IBCSG_hist
#' @docType data
#' @author International Breast Cancer Study Group
#' @references Original data paper (DOI: 10.1200/JCO.1996.14.6.1885)
#' an analysis using the data (DOI: 10.1111/j.1541-0420.2005.00448.x)
#' @keywords data
#' @format A data frame with 103 rows and 8 variables:
#' \describe{
#'   \item{PHYS18}{outcome variable, integer scores between 0 and 100 measuring
#'                the patients' physical well-being on month 18, with a higher
#'                score indicating a better physical well-being}
#'   \item{PHYS1}{physical well-being scores assessed at the start of the study}
#'   \item{N_INIT_CYCLES}{number of initial cycles of oral cyclophosphamide,
#'                        methotrexate, and fluorouracil (CMF), equal to 3 or 6}
#'   \item{REINTRODUCTION}{indicator of reintroduction of chemotherapy, 0 = No
#'                         reintroduction, 1 = having reintroduction}
#'   \item{AGE}{age, in years}
#'   \item{COUNTRY}{country, ANZ = New Zealand/Australia, CH = Switzerland,
#'                  SWED = Sweden}
#'   \item{NODEGRP}{indicator of number of positive nodes being greater than or
#'                  equal to 4, 0 = less than 4, 1 = 4+}
#'   \item{ER}{ER status indicator, 0 = negative, 1 = positive}
#' }
NULL
