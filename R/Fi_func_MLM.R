#' Generate Fisher information matrix F_x at a design point x_i for Multinomial logistic regression model
#'
#' @param X_x model matrix given design point x_i (for example, X_x = h.func(x_i), where h.func transforms a design point to a model matrix)
#' @param beta parameter coefficients in the Multinomial logistic regression model, the order of coefficients in bvec and the order of design points in X_x should be consistent
#' @param link link function of Multinomial logistic regression model, options are "baseline", "cumulative", "adjacent", or "continuation"
#'
#' @return F_x is the Fisher information matrix at design point x_i (with model matrix X_x);
#' @return U_x is a middle step matrix for calculation of F_x, details see Corollary 3.1 in Bu, X., Majumdar, D., & Yang, J. (2020). D-optimal designs for multinomial logistic models
#' @export
#' @examples
#' X_x_temp = rbind(c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#' bvec_temp = c(-4.047, -2.225, -0.302, 1.386, 4.214, 3.519,
#' 2.420, 1.284, -0.131, -0.376, -0.237, -0.120)
#' link_temp = "cumulative"
#' Fi_func_MLM(X_x=X_x_temp, beta=bvec_temp, link=link_temp)



Fi_func_MLM = function(X_x, beta, link){
  ## function to generate fisher information matrix at xi for multinomial generalized linear model #update 2022/09/18
  ## input
  ##        X_x-- design point x_i, X_x=h.func(xi)
  ##        bvec-- beta coefficients in the model
  ##        link-- multinomial generalized linear model link function name "baseline", "cumulative", "adjacent", or"continuation"
  ## output
  ##        F_x-- Fisher information matrix at x_i
  ##        U_x-- U matrix at x_i
  p = length(beta)
  eta_xi =  X_x %*% beta #eta
  J = length(eta_xi)
  pi = rep(NA, J) #calculate pi
  product_element = 1/(exp(eta_xi)+1)
  #adjacent denominator
  sum_element = rep(NA, J-1)
  for(i in 1:(J-1)){
    sum_element[i] = exp(sum(eta_xi[i:(J-1)]))
  }
  for(i in 1:(J-1)){
    if(link=="continuation"){pi[i]=exp(eta_xi[i])*prod(product_element[1:i])} #end of continuation
    if(link == "cumulative"){
      if(i==1){pi[i]=exp(eta_xi[i])/(1+exp(eta_xi[i]))}else{
        pi[i]=((exp(eta_xi[i]))/(1 + exp(eta_xi[i])))-((exp(eta_xi[i-1]))/(1+exp(eta_xi[i-1])))
      }
    }#end of cumulative
    if(link=="baseline"){pi[i]=(exp(eta_xi[i]))/(sum(exp(eta_xi)[1:(J-1)])+1)}#end of baseline
    if(link=="adjacent"){pi[i]=(exp(sum(eta_xi[i:(J-1)])))/(sum(sum_element)+1)}
  }#end of for loop
  if(link=="continuation"){pi[J] = prod(product_element[1:(J-1)])}
  if(link=="cumulative"){pi[J] = 1/(1+exp(eta_xi[J-1]))}
  if(link=="baseline"){pi[J] = 1/(sum(exp(eta_xi)[1:(J-1)])+1)}
  if(link=="adjacent"){pi[J] = 1/(sum(sum_element)+1)}
  #calculate U matrix, U is symmetric matrix
  U = matrix(data = NA, nrow= J, ncol = J)
  U[1:(J-1), J] = 0
  U[J, J] = 1
  #diagonal
  if(link=="continuation"){U[1,1]= pi[1]*(1-pi[1])}
  if(link=="cumulative"){U[1,1]=(pi[1]^2)*((1-pi[1])^2)*(1/pi[1] + 1/pi[2])}
  if(link=="baseline"){U[1,1]=pi[1]*(1-pi[1])}
  if(link=="adjacent"){U[1,1]=pi[1]*(1-pi[1])}
  for(s in 2:(J-1)){
    if(link=="continuation"){U[s,s]=pi[s]*((1-sum(pi[1:s])))*((1-sum(pi[1:(s-1)]))^(-1))}
    if(link=="cumulative"){U[s,s]=(sum(pi[1:s])^2)*((1-sum(pi[1:s]))^2)*(1/pi[s] + 1/pi[s+1])}
    if(link=="baseline"){U[s,s]=pi[s]*(1-pi[s])}
    if(link=="adjacent"){U[s,s]=sum(pi[1:s])*(1-sum(pi[1:s]))}
  }#end of for loop

  #upper triangle
  for(s in 1 : (J-2)){
    for(t in (s+1):(J-1)){
      if(link=="continuation"){U[s,t]=0}#end of continuation
      if(link=="cumulative"){
        if(t - s == 1){U[s,t] = -(sum(pi[1:s])*sum(pi[1:t]))*(1-sum(pi[1:s]))*(1-sum(pi[1:t]))*(1/pi[t])}
        if(t - s > 1){U[s,t]=0}
      }#end of cumulative
      if(link=="baseline"){U[s,t]=-pi[s]*pi[t]}#end of baseline
      if(link=="adjacent"){U[s,t]=sum(pi[1:s])*(1-sum(pi[1:t]))}
    }
  }# end of double for loop

  #symmetric U codes, lower triangle = transpose lower triangle (update 11/17/2024)
  ind_temp=lower.tri(U)
  U[ind_temp]=t(U)[ind_temp]

  F_x = t(X_x) %*% U %*% X_x

  #define S3 class
  output <-list(F_x = F_x, U_x = U)
  class(output)<-"matrix_list"
  return(output)
} #end of Fi.func
