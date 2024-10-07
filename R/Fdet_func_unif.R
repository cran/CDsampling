#' Determinant function to be used for finding constrained uniform samplings
#'
#' @param w allocation (can be exact or approximate)
#' @param beta use NULL (default to be NULL)
#' @param X use NULL (default to be NULL)
#' @param link use NULL (default to be NULL)
#'
#' @return product of all allocation
#' @export
#'
#' @examples
#' Fdet_func_unif(w=c(0.2,0.2,0.2,0.2,0.2))
#'
#'

Fdet_func_unif <- function(w, beta=NULL, X=NULL, link=NULL){
  return(prod(w))
}
