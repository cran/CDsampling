#' Determinant of Fisher information matrix of multinomial logistic model (MLM)
#'
#' @param w allocation (can be exact or approximate)
#' @param beta MLM model covariate coefficient
#' @param X MLM model matrix
#' @param link link function of Multinomial logistic regression model, options are "baseline", "cumulative", "adjacent", or "continuation"
#'
#' @return
#' Determinant of the Fisher information matrix of MLM model
#' @export
#'
#' @examples
#' w = rep(1/8, 8)
#' Xi=rep(0,5*12*8) #response levels * num of parameters * num of design points
#' dim(Xi)=c(5,12,8)
#' #design matrix
#' Xi[,,1] = rbind(c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,2] = rbind(c( 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,3] = rbind(c( 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,4] = rbind(c( 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,5] = rbind(c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
#'                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,6] = rbind(c( 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,7] = rbind(c( 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 1),
#'                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#'
#' Xi[,,8] = rbind(c( 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#' bvec_temp = c(-4.047, -2.225, -0.302, 1.386, 4.214, 3.519,
#' 2.420, 1.284, -0.131, -0.376, -0.237, -0.120)
#' link_temp = "cumulative"
#'
#' Fdet_func_MLM(w=w, beta=bvec_temp, X=Xi, link=link_temp)
#'
#'
#'

Fdet_func_MLM <- function(w, beta, X, link){
  J = dim(X)[1]
  d = dim(X)[2]
  m = dim(X)[3]
  Fi <- rep(0, d*d*m);  dim(Fi)=c(d,d,m)    # F_i matrix = Fi[,,i]
  nFi <- rep(0, d*d*m);  dim(nFi)=c(d,d,m)
  for(i in 1:m) {
    Fi[,,i]=Fi_func_MLM(X_x=X[,,i], beta=beta, link=link)$F_x
    nFi[,,i]=Fi[,,i]*w[i]
  }

  Fisher_matrix=apply(nFi,c(1,2),sum)        # F = sum_i n_i*F_i
  Fdet=det(Fisher_matrix) # |F| at (p_1, p_2,...,p_m)=(n_1,...,n_m)/n
  return(Fdet)
}


