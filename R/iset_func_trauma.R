#' trauma_data example (see Huang, Tong, Yang (2023)) specific function for finding index set that if allocation of that index add "1", the new allocation still falls within the constraint
#' Used in approxtoexact_constrained_func()
#'
#' @param allocation the exact allocation
#'
#' @return
#' list of TRUE and FALSE, if TRUE, it means the allocation of this index will fall out of the constraint with more subject; if TURE, it means the allocation of this index can add more subjects
#' @export
#'
#' @examples
#' iset_func_trauma(allocation=c(50,30,10,10,100,100,200,10))
#'
#'
#'

iset_func_trauma <- function(allocation){
  iset = rep(TRUE,8)
  if(sum(allocation[1:4])>=392){iset[1:4]=FALSE}
  if(sum(allocation[5:8])>=410){iset[5:8]=FALSE}
  return(iset)
}
