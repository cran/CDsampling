#' Unconstrained lift-one algorithm to find D-optimal allocations for GLM
#' @importFrom stats rexp
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param W Diagonal of W matrix in Fisher information matrix, can be calculated from W_func_GLM in package
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 500
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of random initial points, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#' @param w00 Specified initial design proportion; default to be NULL, this will generate a random initial design
#' @param label A vector of text strings for subgroups' names, default value NULL
#'
#' @return w is the approximate D-optimal design
#' @return w0 is the initial design used to get optimal design w
#' @return Maximum is the maximized |F| value
#' @return itmax is the number of iterations
#' @return convergence is TRUE or FALSE, if TRUE means the reported design is converged
#' @export
#'
#' @examples
#' beta = c(0.5, 0.5, 0.5)
#' X = matrix(data=c(1,-1,-1,1,-1,1,1,1,-1), byrow=TRUE, nrow=3)
#' W_matrix = W_func_GLM(X=X, beta=beta)
#' w00 = c(1/6, 1/6, 2/3)
#' approximate_design = liftone_GLM(X=X, W=W_matrix, reltol=1e-10, maxit=100,
#' random=FALSE, nram=3, w00=w00)
#'



liftone_GLM <- function(X, W, reltol=1e-5, maxit=500, random=TRUE, nram=3, w00=NULL, label=NULL)  {   ## W=W[1,2,...,m] are strictly positive
  # if random=T, run 5 random initial points and pick up the best; default initial p1=p2=...=1/m
  # output: w=w--optimal design based on "det"
  #         Maximum--maximized value of "det"
  #         convergence -- "T" indicates success
  #         w0 -- initial weight
  #         itmax -- number of iterations
  m = dim(X)[1];
  d = dim(X)[2];
  W = W[1:m];
  if(min(W) <= 0) {
    message("\nW's need to be strictly positive!\n");
    return(0);
  };
  ftemp <- function(p) { det(t(X * (p*W)) %*% X);};
  if(is.null(w00)) w00=rep(1/m,m);
  maximum = ftemp(w00);
  maxvec = rexp(m);
  convergence=FALSE;
  p = w00;
  ind = 0;
  while((ind < maxit) && ((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 > reltol)) {
    io = sample(x=(1:m), size=m);
    for(i in 1:m) {    # run updating in random order of w
      if(p[io[i]]>0) {
        ptemp1 = p/(1-p[io[i]]);
        ptemp1[io[i]] = 0;
        b = ftemp(ptemp1);      # b=fs(0)
        a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
      } else {         # p[io[i]]==0
        b = maximum;
        ptemp1 = p/2;
        ptemp1[io[i]] = 1/2;          # for fs(1/2)
        a = ftemp(ptemp1)*2^d - b;
      }
      if(a > b*d) x=(a-b*d)/((a-b)*d) else x=0;
      ptemp1 = p*(1-x)/(1-p[io[i]]);
      ptemp1[io[i]] = x;
      if(a > b*d) maximum = ((d-1)/(a-b))^(d-1)*(a/d)^d else maximum=b;
      p = ptemp1;
      maxvec[io[i]] = maximum;
    }
    ind = ind+1;
  }
  p.ans=p; maximum.ans=maximum; if((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 <= reltol) convergence=TRUE;itmax=ind;
  if(random) for(j in 1:nram) {
    p0=rexp(m);p0=p0/sum(p0);
    p=p0;
    maxvec = rexp(m);
    maximum = ftemp(p);
    ind = 0;
    while((ind < maxit) && ((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 > reltol)) {
      io = sample(x=(1:m), size=m);
      for(i in 1:m) {    # run updating in random order of w
        if(p[io[i]]>0) {
          ptemp1 = p/(1-p[io[i]]);
          ptemp1[io[i]] = 0;
          b = ftemp(ptemp1);      # b=fs(0)
          a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
        } else {         # p[io[i]]==0
          b = maximum;
          ptemp1 = p/2;
          ptemp1[io[i]] = 1/2;          # for fs(1/2)
          a = ftemp(ptemp1)*2^d - b;
        }
        if(a > b*d) x=(a-b*d)/((a-b)*d) else x=0;
        ptemp1 = p*(1-x)/(1-p[io[i]]);
        ptemp1[io[i]] = x;
        if(a > b*d) maximum = ((d-1)/(a-b))^(d-1)*(a/d)^d else maximum=b;
        p = ptemp1;
        maxvec[io[i]] = maximum;
      }
      ind = ind+1;
    }
    if(maximum > maximum.ans) {
      maximum.ans=maximum;
      p.ans=p;
      convergence=FALSE;
      if((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 <= reltol) convergence=TRUE;
      w00=p0;
      itmax=ind;
    }
  }
  #define S3 class
  if((max(maxvec,na.rm=T)/min(maxvec,na.rm=T))-1 <= reltol) convergence=TRUE;
  output <- list(w=p.ans, w0=w00, Maximum=maximum.ans, itmax=itmax, convergence=convergence, label=label);  # convergence=T indicates success
  class(output)<-"list_output"
  return(output)
}
