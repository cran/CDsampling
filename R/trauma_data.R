#' Trauma data with multinomial response
#'
#' The data frame saves data from the trauma trial data from Chuang-Stein and Agresti (1997).
#'
#' @format A data frame with 802 rows and 5 variables:
#'  \describe{
#'        \item{Severity}{severity of the trauma symptoms, mild or moderate/severe}
#'        \item{Dose}{dose levels applied to the patients, 4 levels, placebo, low, medium and high}
#'        \item{Label}{stratification group in terms of severity and dose}
#'        \item{Outcome}{treatment outcome, 5 levels, death, vegetative state, major disability, minor disability and good recovery}
#'        \item{ID}{patient ID, 1-802}
#' }
#' @source {Chuang-Stein and Agresti (1997)}
#' @examples
#' data(trauma_data) #lazy loading
"trauma_data"
