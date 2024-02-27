#' ECOG E2696 Trial
#'
#' A data set from the ECOG E2696 trial comparing the combination of the GM2-KLH/QS-21 (GMK) vaccine
#' and high-dose interferon alfa-2b (IFN) therapy with the GMK vaccine alone in resected high-risk
#' melanoma patients. The study results were described in Kirkwood et al. (2001) <doi:10.1200/JCO.2001.19.5.1430>.
#' This data set only includes patients without nodal metastasis.
#'
#' @name E2696
#' @docType data
#' @usage data(E2696)
#' @keywords data
#' @format A data frame with 105 rows and 6 variables:
#' \describe{
#'   \item{failtime}{relapse-free survival (RFS) times (in months)}
#'   \item{failind}{relapse indicator, 0 = right censored, 1 = relapse}
#'   \item{treatment}{treatment indicator, 0 = GMK, 1 = GMK and IFN}
#'   \item{sex}{gender indicator, 0 = male, 1 = female}
#'   \item{age}{patient age in years}
#'   \item{perform}{ECOG performance status indicator, 0 = fully active patient, able to
#'                  carry on all pre-disease performance without restriction,
#'                  1 = restricted in physically strenuous activity, but are ambulatory
#'                  and able to carry out work of a light or sedentary nature}
#' }
#' @references
#'  Kirkwood, J. M., Ibrahim, J., Lawson, D. H., Atkins, M. B., Agarwala, S. S., Collins, K., Mascari, R., Morrissey, D. M., and Chapman, P. B. (2001). High-dose interferon alfa-2b does not diminish antibody response to GM2 vaccination in patients with resected melanoma: Results of the multicenter eastern cooperative oncology group phase II trial E2696. Journal of Clinical Oncology, 19(5), 1430–1436.
NULL
