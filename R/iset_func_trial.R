#' trial_data example (see Huang, Tong, Yang (2023)) specific function for finding index set that if allocation of that index add "1", the new allocation still falls within the constraint
#' Used in approxtoexact_constrained_func()
#'
#' @param allocation the exact allocation
#'
#' @return
#' list of TRUE and FALSE, if TRUE, it means the allocation of this index will fall out of the constraint with more subject; if TURE, it means the allocation of this index can add more subjects
#' @export
#'
#' @examples
#' iset_func_trial(allocation=c(50,30,10,100,100,40))
#'


iset_func_trial <- function(allocation){
  Ni = c(50, 40, 10, 200, 150, 50)
  return(allocation < Ni)
}
