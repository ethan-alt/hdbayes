#' ECOG E1690 Trial
#'
#' A data set from the ECOG E1690 trial evaluating the effectiveness of the interferon alfa-2b
#' (IFN) therapy compared to the observation in high-risk melanoma patients. There were three
#' arms in the trial: high-dose IFN, low-dose IFN, and observation. The study results were
#' described in Kirkwood et al. (2000) <doi:10.1200/JCO.2000.18.12.2444>. Here, we only consider
#' the high-dose IFN arm and the observation arm so that this data set has the same variables as
#' the E1684 data set. We can use the E1684 data as the historical data and the E1690 data as the
#' current data.
#'
#' @name E1690
#' @docType data
#' @usage E1690
#' @keywords data
#' @format A data frame with 426 rows and 8 variables:
#' \describe{
#'   \item{failtime}{time to relapse in years}
#'   \item{failcens}{censoring indicator for time to relapse, 0 = did not relapse, 1 = relapsed}
#'   \item{survtime}{time to death in years}
#'   \item{survcens}{censoring indicator for time to death, 0 = alive, 1 = dead}
#'   \item{treatment}{treatment indicator, 0 = observation, 1 = high-dose IFN}
#'   \item{sex}{gender indicator, 0 = male, 1 = female}
#'   \item{age}{patient age in years}
#'   \item{node_bin}{indicator for having more than one cancerous lymph node,
#'               0 = with one or no cancerous lymph nodes,
#'               1 = with more than one cancerous lymph node}
#' }
#' @references
#'  Kirkwood, J. M., Ibrahim, J. G., Sondak, V. K., Richards, J., Flaherty, L. E., Ernstoff, M. S., Smith, T. J., Rao, U., Steele, M., and Blum, R. H. (2000). High- and low-dose interferon alfa-2b in high-risk melanoma: First analysis of intergroup trial E1690/S9111/C9190. Journal of Clinical Oncology, 18(12), 2444–2458.
NULL
