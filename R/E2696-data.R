#' ECOG E2696 Data
#' from patients without nodal metastases
#'
#' @name E2696
#' @docType data
#' @author Joseph Ibrahim \email{ibrahim@@bios.unc.edu}
#' @references DOI: 10.1200/JCO.2001.19.5.1430
#' @keywords data
#' @format A data frame with 105 rows and 6 variables:
#' \describe{
#'   \item{failtime}{relapse-free surviva (RFS) times (in months)}
#'   \item{failind}{event indicator, 0 = censored, 1 = failed}
#'   \item{trt}{treatment indicator, 0 = GMK, 1 = IFN and GMK}
#'   \item{sex}{gender indicator, 0 = male, 1 = female}
#'   \item{age}{age, in years}
#'   \item{perform}{ECOG performance status indicator, 0 = fully active patient, able to
#'                  carry on all pre-disease performance without restriction,
#'                  1 = restricted in physically strenuous activity, but are ambulatory
#'                  and able to carry out work of a light or sedentary nature}
#' }
NULL
