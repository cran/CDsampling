#' Convert the approximate allocation (proportion) to exact allocation (integer) with bounded constraint (ni <= Ni)
#'
#' @param n Sample size, must be a positive integer
#' @param w Approximate allocation/proportion, must be a real-valued vector, can get from running liftone_constrained_GLM or liftone_constrained_MLM
#' @param m The number of sampling groups
#' @param beta Model parameter coefficients, default to be NULL for use in constrained uniform sampling
#' @param X Design matrix of the model for GLM or MLM, default to be NULL for use in constrained uniform sampling
#' @param link Link function of GLM or MLM, if used for GLM model (GLM_T is T), options are "identity", "logit", "probit", "cloglog", "loglog". If used for MLM (GLM_T is F), options are "continuation", "cumulative", "adjacent", and "baseline"
#' @param Fdet_func determinant of Fisher information matrix function, Fdet_func can be self-defined, or use "Fdet_func_GLM", "Fdet_func_MLM" in the package, default is Fdet_func_GLM
#' @param iset_func self-defined function for checking which index of sampling group fall within constraint if add 1 more subject (I set, see Algorithm 2 in Huang, Tong, Yang (2023)), two example functions are provided in the package, iset_func_trial and iset_func_trauma
#' @param label A vector of text strings for subgroups' names, default value NULL
#'
#' @return allocation is the exact allocation or integer value of the number of subjects sampled from the group
#' @return allocation.real is the proportion or the approximate allocation of the number of subjects sampled from the group
#' @return det.maximum is the maximum of |F| from the current exact allocation
#' @export
#'
#' @examples
#'
#' beta = c(0, 3, 3, 3) #main effect model beta_0, beta_1, beta_21, beta_22
#' X.liftone=matrix(data=c(1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1), ncol=4, byrow=TRUE)
#' exact_design = approxtoexact_constrained_func(n=200, w=c(0.25, 0.20, 0.05, 0.50, 0.00, 0.00),
#' m=6, beta=beta, link='logit', X=X.liftone, Fdet_func=Fdet_func_GLM, iset_func=iset_func_trial)
#'

approxtoexact_constrained_func <- function(n, w, m, beta=NULL, link=NULL, X=NULL, Fdet_func=Fdet_func_GLM, iset_func=NULL, label=NULL) {
  # n > 0 is the targeted sample size, must be a positive integer
  # w=w[1:m] is a real-valued approximate allocation, w_i >=0, sum(w_i) = 1, w_i <= Ni/n
  # X is the design matrix for GLM it is m*d matrix, for MLM it is J*d*m matrix
  # beta is the parameter
  # link is the link function of GLM or MLM
  # m is the number of sampling groups
  # d is the number of parameter
  # J is the response levels for MLM, for GLM use NULL
  # Ni[1:m] are constraints
  # F_func is the Fisher information matrix function of design, can use F_func_GLM, F_func_MLM, a special case is to realize constrained uniform design, then F_func will be F_func_unif which is the product of all allocation, with initial approximate allocation of a list of 0
  # output: allocation=integer-valued allocation
  #         allocation.real[1:m]=real-valued allocation
  #         det.maximum--maximized value of "det"

  w = w[1:m];
  if(min(w) < 0) {message("\nw_i needs to be nonnegative!\n"); return(0);};
  if(sum(w)==0){message("\n summation of w can not be 0!\n")};
  #w = w/sum(w)
  ftemp=function(nvec){
    Fdet = Fdet_func(w=nvec, beta=beta, X=X, link=link)
    return(Fdet)                   # |F| at (w_1, w_2,...,w_m)=(n_1,...,n_m)/n
  }

  allocation=floor(n*w);
  det.maximum = ftemp(allocation);
  k=n-sum(allocation);
  if(is.null(iset_func)){iset=rep(1,m)}else{iset = iset_func(allocation)};      # I
  while(k>0) {
    dtemp=rep(0,m);                  # d1,...,dm
    for(i in 1:m) if(iset[i]==1) {
      ntemp=allocation;
      ntemp[i]=ntemp[i]+1;
      dtemp[i]=ftemp(ntemp);
    };
    det.maximum=max(dtemp);
    istar=which.max(dtemp);
    allocation[istar]=allocation[istar]+1;
    k=k-1;
    if(is.null(iset_func)){iset=rep(1,m)}else{iset = iset_func(allocation)}; #check new allocation for I set
  };

  #define S3 class
  output<-list(allocation=allocation, allocation.real=w, det.maximum=det.maximum, label=label);
  class(output)<-"list_output"
  return(output)
}

