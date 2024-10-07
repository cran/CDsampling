#' Calculate the diagonal elements nu of Fisher information matrix
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @param X Model matrix
#' @param beta Parameters of GLM model
#' @param link GLM link function, default is "logit", options are "logit", "probit", "cloglog", "loglog", "identity", identity is the same as ordinary linear regression
#'
#' @return the diagonal element nu of GLM Fisher information matrix, can be used as w in liftone_constrained_GLM
#' @export
#'
#' @examples
#' beta = c(0, 3, 3, 3) #main effect model beta_0, beta_1, beta_21, beta_22
#' #gives the 6 categories (0,0,0), (0,1,0),(0,0,1),(1,0,0),(1,1,0),(1,0,1)
#' X.liftone=matrix(data=c(1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1), ncol=4, byrow=TRUE)
#' #calculate diagonal elements of W based on beta's under logit link
#' W=W_func_GLM(X= X.liftone, beta=beta, link="logit")
#'

W_func_GLM = function(X,beta,link="logit"){
  if(dim(X)[2] != length(beta)) {
    message("\nThe dimensions do not match!\n");
    return(0);
  };
  if(link=="identity"){
    return(1)
  }
  if(link=="logit"){
    bx = as.vector(X %*% beta);  # b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4
    return(exp(bx)/(1+exp(bx))^2);
  }
  if(link=="probit"){
    bx = as.vector(X %*% beta);  # b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4
    return(((dnorm(bx))^2) / (pnorm(bx)*(1-pnorm(bx))));
  }
  if(link=="cloglog"){
    bx = as.vector(X %*% beta);  # b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4
    return((exp(2*bx - exp(bx)))/(1-exp(-exp(bx))));
  }
  if(link=="loglog"){
    bx = as.vector(X %*% beta);  # b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4
    return((exp(2*bx - exp(bx)))/(1-exp(-exp(bx))));
  }

}
