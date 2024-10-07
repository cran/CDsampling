#' Generated clinical trial data with binary response
#'
#' Generated with logistic regression model:
#'          \deqn{logit(P(Y_{ij} = 1 | gender_i, age_{i1}, age_{i2})) = 3*gender_i +3*age_{i1} +3*age_{i2}}
#' The data frame can be used to run GLM clinical trial example in Huang, Tong, Yang (2023)
#'
#' @format A data frame with 500 rows and 6 variables:
#'  \describe{
#'        \item{gender}{gender of the patients}
#'        \item{age_1}{1 or 0, whether or not the patient belongs to 18-25 age group}
#'        \item{age_2}{1 or 0, whether or not the patient belongs to 26-64 age group}
#'        \item{label}{stratification group in terms of gender and age, 1 to 6}
#'        \item{Y}{treatment effective or not, Y=1 means treatment is effective to the patient}
#'        \item{ID}{patient ID, 1-500}
#' }
#' @source {Generated pseudo clinical trial data to serve as an example.}
#' @examples
#' data(trial_data) #lazy loading
#'
"trial_data"
