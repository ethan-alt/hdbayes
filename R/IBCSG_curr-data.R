#' International Breast Cancer Study Group (IBCSG) Trial VI Data
#'
#' A data set from the IBCSG Trial VI investigating both the duration of adjuvant chemotherapy
#' (3 versus 6 initial cycles of oral cyclophosphamide, methotrexate, and fluorouracil (CMF)) and
#' the reintroduction of single courses of delayed chemotherapy in node-positive premenopausal
#' breast cancer patients. The study results were described by IBCSG (1996) <doi:10.1200/JCO.1996.14.6.1885>
#' and Hürny et al. (1992) <doi:10.1016/0959-8049(92)90399-m>. This data set only includes patients
#' above the age of 40 (i.e., age \eqn{\geq} 40) and treats the measurements of patients' physical well-being
#' on month 18 as the outcome. The IBCSG_hist data set includes patients from the same study but
#' with age < 40. We can use the IBCSG_hist data as the historical data and the IBCSG_curr data as
#' the current data.
#'
#' @name IBCSG_curr
#' @docType data
#' @usage IBCSG_curr
#' @keywords data
#' @format A data frame with 488 rows and 8 variables:
#' \describe{
#'   \item{phys18}{outcome variable, integer scores between 0 and 100 measuring
#'                the patients' physical well-being on month 18, with a higher
#'                score indicating a better physical well-being}
#'   \item{phys1}{physical well-being scores assessed at the start of the study}
#'   \item{n_init_cycles}{number of initial cycles of CMF, equal to 3 or 6}
#'   \item{reintroduction}{indicator of reintroduction of chemotherapy, 0 = no
#'                         reintroduction, 1 = having reintroduction}
#'   \item{age}{patient age in years}
#'   \item{country}{country, ANZ = New Zealand/Australia, CH = Switzerland,
#'                  SWED = Sweden}
#'   \item{nodegp}{indicator of number of positive nodes being greater than or
#'                  equal to 4, 0 = less than 4, 1 = 4+}
#'   \item{ER}{estrogen receptor (ER) status indicator, 0 = negative, 1 = positive}
#' }
#' @references
#'  International Breast Cancer Study Group. (1996). Duration and reintroduction of adjuvant chemotherapy for node-positive premenopausal breast cancer patients. Journal of Clinical Oncology, 14(6), 1885–1894.
#'
#'  Hürny, C., Bernhard, J., Gelber, R. D., Coates, A., Castiglione, M., Isley, M., Dreher, D., Peterson, H., Goldhirsch, A., and Senn, H.-J. (1992). Quality of life measures for patients receiving adjuvant therapy for breast cancer: An international trial. European Journal of Cancer, 28(1), 118–124.
#'
#'  Chi, Y.-Y. and Ibrahim, J. G. (2005). Joint models for multivariate longitudinal and Multivariate Survival Data. Biometrics, 62(2), 432–445.
NULL
