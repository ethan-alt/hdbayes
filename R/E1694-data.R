#' ECOG E1694 Data
#' from patients receiving the vaccine (GMK) only
#'
#' @name E1694
#' @docType data
#' @author Joseph Ibrahim \email{ibrahim@@bios.unc.edu}
#' @references DOI: 10.1200/JCO.2001.19.9.2370
#' @keywords data
#' @format A data frame with 167 rows and 4 variables:
#' \describe{
#'   \item{log_igm28}{continuous outcome variable, log of antibody immunoglobulin (IgM)
#'                    levels on day 28, more precisely, log(IgM + 1) on day 28}
#'   \item{age}{age, in years}
#'   \item{sex}{gender indicator, 0 = male, 1 = female}
#'   \item{perform}{ECOG performance status indicator, 0 = fully active patient, able to
#'                  carry on all pre-disease performance without restriction,
#'                  1 = restricted in physically strenuous activity, but are ambulatory
#'                  and able to carry out work of a light or sedentary nature}
#' }
NULL
