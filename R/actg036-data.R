#' AIDS Clinical Trial ACTG036
#'
#' A data set from the AIDS clinical trial ACTG036 (\url{https://clinicaltrials.gov/study/NCT00001104}) comparing
#' zidovudine (AZT) with a placebo in patients with hereditary coagulation disorders and HIV infection. The study
#' results were described in Merigan et al. (1991) <doi:10.1182/blood.V78.4.900.900>. This data set has the same
#' variables as the actg019 data set. We can use the actg019 data as the historical data and the actg036 data as
#' the current data.
#'
#' @name actg036
#' @docType data
#' @usage actg036
#' @keywords data
#' @format A data frame with 183 rows and 5 variables:
#' \describe{
#'   \item{outcome}{outcome variable with 1 indicating death, development of AIDS or AIDS-related complex (ARC) and 0 otherwise}
#'   \item{age}{patient age in years}
#'   \item{treatment}{treatment indicator, 0 = placebo, 1 = AZT}
#'   \item{race}{race indicator, 0 = non-white, 1 = white}
#'   \item{cd4}{CD4 cell count}
#' }
#' @references
#'  Merigan, T., Amato, D., Balsley, J., Power, M., Price, W., Benoit, S., Perez-Michael, A., Brownstein, A., Kramer, A., and Brettler, D. (1991). Placebo-controlled trial to evaluate zidovudine in treatment of human immunodeficiency virus infection in asymptomatic patients with hemophilia. NHF-ACTG 036 Study Group. Blood, 78(4), 900–906.
#'
#'  Chen, M.-H., Ibrahim, J. G., and Yiannoutsos, C. (1999). Prior elicitation, Variable Selection and Bayesian computation for Logistic Regression Models. Journal of the Royal Statistical Society Series B: Statistical Methodology, 61(1), 223–242.
NULL
