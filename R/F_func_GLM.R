#' Fisher information matrix of generalized linear model (GLM)
#'
#' @param w  allocation (can be exact or approximate)
#' @param beta GLM model covariate coefficient
#' @param X model matrix
#' @param link link function, default"logit", choose from "logit", "cloglog", "loglog", "probit", and "identity"(for regular linear regression)
#'
#' @return
#' Fisher information matrix given X and model parameter beta
#' @export
#'
#' @examples
#' w = c(1/3,1/3, 1/3)
#' beta = c(0.5, 0.5, 0.5)
#' X = matrix(data=c(1,-1,-1,1,-1,1,1,1,-1), byrow=TRUE, nrow=3)
#' Fdet_func_GLM(w=w, beta=beta, X=X, link='logit')
#'
#'
#'


F_func_GLM = function(w, beta, X, link="logit"){
  W_matrix=W_func_GLM(X=X, beta=beta, link=link)
  Fisher_matrix = t(X * (w*W_matrix)) %*% X
  return(Fisher_matrix)
}
