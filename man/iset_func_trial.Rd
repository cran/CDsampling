% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iset_func_trial.R
\name{iset_func_trial}
\alias{iset_func_trial}
\title{trial_data example (see Huang, Tong, Yang (2023)) specific function for finding index set that if allocation of that index add "1", the new allocation still falls within the constraint
Used in approxtoexact_constrained_func()}
\usage{
iset_func_trial(allocation)
}
\arguments{
\item{allocation}{the exact allocation}
}
\value{
list of TRUE and FALSE, if TRUE, it means the allocation of this index will fall out of the constraint with more subject; if TURE, it means the allocation of this index can add more subjects
}
\description{
trial_data example (see Huang, Tong, Yang (2023)) specific function for finding index set that if allocation of that index add "1", the new allocation still falls within the constraint
Used in approxtoexact_constrained_func()
}
\examples{
iset_func_trial(allocation=c(50,30,10,100,100,40))

}
