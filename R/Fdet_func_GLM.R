#' Determinant of Fisher information matrix for GLM
#'
#' @param w  allocation (can be exact or approximate)
#' @param beta GLM model covariate coefficient
#' @param X model matrix
#' @param link link function, default"logit", choose from "logit", "cloglog", "loglog", "probit", and "identity"(for regular linear regression)
#'
#' @return
#' the determinant of Fisher information matrix given X and model parameter beta
#' @export
#'
#' @examples
#' w = c(1/3,1/3, 1/3)
#' beta = c(0.5, 0.5, 0.5)
#' X = matrix(data=c(1,-1,-1,1,-1,1,1,1,-1), byrow=TRUE, nrow=3)
#' Fdet_func_GLM(w=w, beta=beta, X=X, link='logit')
#'


Fdet_func_GLM = function(w, beta, X, link="logit"){
  W_matrix=W_func_GLM(X=X, beta=beta, link=link)
  Fisher_matrix = t(X * (w*W_matrix)) %*% X
  det_Fisher = det(Fisher_matrix)
  return(det_Fisher)
}
