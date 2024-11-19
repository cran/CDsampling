#' Unconstrained lift-one algorithm to find D-optimal allocations for MLM
#' @importFrom stats rexp
#' @importFrom stats optim
#' @param m The number of design points; it is usually the number of combinations of all the stratification factors
#' @param p The number of parameters in the MLM model
#' @param Xi Model matrix, a J by p by m 3D array of predictors for separate response category at all design points(input to determine ppo,npo,po)
#' @param J The number of response levels
#' @param beta A p*1 vector, parameter coefficients for MLM, the order of beta should be consistent with Xi
#' @param link Link function of MLM, default to be "cumulative", options from "continuation", "cumulative", "adjacent", and "baseline"
#' @param Fi.func A function for calculating Fisher information at a specific design point, default to be Fi_func_MLM function in the package
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 500
#' @param w00 Specified initial design proportion; default to be NULL, this will generate a random initial design
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of random initial points, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#'
#' @return w is the approximate D-optimal design
#' @return w0 is the initial design used to get optimal design
#' @return Maximum is the maximized |F| value
#' @return itmax is the number of iterations
#' @return convergence is TRUE or FALSE, if TRUE means the reported design is converged
#' @export
#'
#' @examples
#' J = 5    # number of categories,  >= 3
#' p = 12    # number of parameters
#' m = 8    # number of design points
#' Xi=rep(0,J*p*m) #J*p*m=5*12*8
#' dim(Xi)=c(J,p,m)
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
#'
#' thetavec = c(-4.047, -0.131, 4.214, -2.225, -0.376, 3.519,
#' -0.302, -0.237,  2.420, 1.386,  -0.120,  1.284)
#'
#' liftone_MLM(m=m, p=p, Xi=Xi, J=J, beta=thetavec, link = "cumulative",
#' Fi.func=Fi_func_MLM, reltol=1e-5, maxit=5000, w00=NULL, random=TRUE, nram=3)
#'
#'




liftone_MLM <- function(m, p, Xi, J, beta, link = "continuation", Fi.func=Fi_func_MLM, reltol=1e-5, maxit=500, w00=NULL, random=TRUE, nram=3) {
  if(is.null(w00)){w00=rexp(m); w00=w00/sum(w00);}
  Fi <- rep(0, p*p*m);  dim(Fi)=c(p,p,m)
  nFi <- rep(0, p*p*m);  dim(nFi)=c(p,p,m)

  for(i in 1:m) {
    Fi[,,i]=Fi.func(Xi[, ,i], beta=beta, link=link)$F_x
    nFi[,,i]=w00[i]*Fi[,,i]
  }
  F=apply(nFi,c(1,2),sum)
  Fdet=det(F)

  Bn1 <- matrix(1, J-1, J-1)     # B^{-1}
  for(j in 2:(J-1)) Bn1[,j]=(1:(J-1))^(j-1);
  Bn1 = solve(Bn1);
  fdet <- function(p) {   # |F|=|sum_i p_i F_i|, p[1:m], need "Fi"
    atemp=p[1]*Fi[,,1];
    for(i in 2:m) atemp=atemp+p[i]*Fi[,,i];
    det(atemp);
  }
  fiz <- function(z, p, i) {   # f_i(z), need "fdet"
    p1=p*(1-z)/(1-p[i]);
    p1[i]=z;
    fdet(p1);
  }
  if(is.null(w00)) w00=w00;        # default initial point is uniform design
  maximum = fdet(w00);
  maxvec = rexp(m);
  convergence=FALSE;
  p0 = w00;
  ind = 0;

  while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
    io = sample(seq(1:m));
    for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
      avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
      avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
      cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
      for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
      avec[J:2]=Bn1%*%cvec;
      ftemp <- function(z) {   # -f_i(z)
        obj=-(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
        # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
        return(obj)
      }
      ftemp1 <- function(z) {  # -f'_i(z)
        #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
        -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
      }
      temp=optim(par=0.5, fn=ftemp, gr=ftemp1, method="L-BFGS-B", lower=0, upper=1, control=list(maxit=maxit, factr=1e5));
      zstar=temp$par;         # z_*
      fstar=-temp$value;
      if(fstar <= avec[1]) {zstar=0; fstar=avec[1];};
      ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
      ptemp1[io[ia]] = zstar;
      if(fstar > maximum) {maximum = fstar; p0=ptemp1;};
      maxvec[io[ia]] = maximum;
    }
    ind = ind+1;
    #cat("\nmaxit", maxit, "\nmax(maxvec)", max(maxvec), "\nmin(maxvec)", min(maxvec)) #delete
  }# end of "while"
  w00.ans = w00;
  p0.ans=p0;
  maximum.ans=maximum;
  #maximum.adj=maximum*n^p;
  #fdet.adj=Fdet*n^p;
  if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;
  itmax=ind;
  effi=(Fdet/maximum.ans)^(1/p)
  #random initial weights
  if(random){
    for(j in 1:nram){
      w00=rexp(m)
      w00=w00/sum(w00)
      Fi <- rep(0, p*p*m);  dim(Fi)=c(p,p,m)
      nFi <- rep(0, p*p*m);  dim(nFi)=c(p,p,m)

      for(i in 1:m) {
        Fi[,,i]=Fi.func(Xi[, ,i], beta=beta, link=link)$F_x
        nFi[,,i]=w00[i]*Fi[,,i]
      }
      F=apply(nFi,c(1,2),sum)
      Fdet=det(F)

      Bn1 <- matrix(1, J-1, J-1)     # B^{-1}
      for(j in 2:(J-1)) Bn1[,j]=(1:(J-1))^(j-1);
      Bn1 = solve(Bn1);
      fdet <- function(p) {   # |F|=|sum_i p_i F_i|, p[1:m], need "Fi"
        atemp=p[1]*Fi[,,1];
        for(i in 2:m) atemp=atemp+p[i]*Fi[,,i];
        det(atemp);
      }
      fiz <- function(z, p, i) {   # f_i(z), need "fdet"
        p1=p*(1-z)/(1-p[i]);
        p1[i]=z;
        fdet(p1);
      }
      if(is.null(w00)) w00=w00;        # default initial point is uniform design
      maximum = fdet(w00);
      maxvec = rexp(m);
      convergence=FALSE;
      p0 = w00;
      ind = 0;

      while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
        io = sample(seq(1:m));
        for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
          avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
          avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
          cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
          for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
          avec[J:2]=Bn1%*%cvec;
          ftemp <- function(z) {   # -f_i(z)
            obj=-(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
            # cat("\navec", avec, "\nz",z, "\nsum", sum(avec*z^(0:(J-1))*(1-z)^((J-1):0)),"\nobject",obj,"\n") #delete
            return(obj)
          }
          ftemp1 <- function(z) {  # -f'_i(z)
            #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
            -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
          }
          temp=optim(par=0.5, fn=ftemp, gr=ftemp1, method="L-BFGS-B", lower=0, upper=1, control=list(maxit=maxit, factr=1e5));
          zstar=temp$par;         # z_*
          fstar=-temp$value;
          if(fstar <= avec[1]) {zstar=0; fstar=avec[1];};
          ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
          ptemp1[io[ia]] = zstar;
          if(fstar > maximum) {maximum = fstar; p0=ptemp1;};
          maxvec[io[ia]] = maximum;
        }
        ind = ind+1;
        #cat("\nmaxit", maxit, "\nmax(maxvec)", max(maxvec), "\nmin(maxvec)", min(maxvec)) #delete
      }# end of "while"
      if(maximum > maximum.ans){
        w00.ans = w00;
        p0.ans=p0;
        maximum.ans=maximum;
        #maximum.adj=maximum*n^p;
        #fdet.adj=Fdet*n^p;
        if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;
        itmax=ind;
        #effi=(Fdet/maximum.ans)^(1/p)
      }

    }#end of for loop of 1:nram

  }#end of if(random)

  #assign S3 class
  if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;
  output<-list(w=p0.ans, w0=w00.ans, Maximum=maximum.ans, convergence=convergence, itmax=itmax);
  class(output)<-"list_output"
  return(output)
}
