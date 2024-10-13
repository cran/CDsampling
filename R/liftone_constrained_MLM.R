#' Find constrained D-optimal designs for Multinomial Logit Models (MLM)
#' @importFrom stats rexp
#' @importFrom stats uniroot
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom lpSolve lp
#' @importFrom Rglpk Rglpk_solve_LP
#' @param m The number of design points; it is usually the number of combinations of all the stratification factors
#' @param p The number of parameters in the MLM model
#' @param Xi Model matrix, a J by p by m 3D array of predictors for separate response category at all design points(input to determine ppo,npo,po)
#' @param J The number of response levels
#' @param beta A p*1 vector, parameter coefficients for MLM, the order of beta should be consistent with Xi
#' @param lower.bound A function to determine lower bound r_i1 in Step 3 of Constrained lift-one algorithm from Yifei, H., Liping, T., Yang, J. (2023) Constrained D-optimal design for paid research study
#' @param upper.bound A function to determine upper bound r_i2 in Step 3 of Constrained lift-one algorithm from Yifei, H., Liping, T., Yang, J. (2023) Constrained D-optimal design for paid research study
#' @param g.con A matrix of numeric constraint coefficients, one row per constraint, on column per variable (to be used in as const.mat lp() and mat in Rglpk_solve_LP())
#' @param g.dir Vector of character strings giving the direction of the constraint: each value should be one of "<," "<=," "=," "==," ">," or ">=". (In each pair the two values are identical.) to be used as const.dir in lp() and dir in Rglpk_solve_LP()
#' @param g.rhs Vector of numeric values for the right-hand sides of the constraints. to be used as const.rhs in lp() and rhs in Rglpk_solve_LP()
#' @param w00 Specified initial design proportion; default to be NULL, this will generate a random initial design
#' @param link Link function of MLM, default to be "cumulative", options from "continuation", "cumulative", "adjacent", and "baseline"
#' @param Fi.func A function for calculating Fisher information at a specific design point, default to be Fi_func_MLM function in the package
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 500
#' @param delta A very small number, used in alpha_star calculation, default to be 1e-6.
#' @param epsilon A very small number, for comparison of >0, <0, ==0, to reduce errors, default 1e-12
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of random initial points, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#'
#' @return w is the approximate D-optimal design
#' @return w0 is the initial design used to get optimal design w
#' @return Maximum is the maximized |F| value
#' @return itmax is the number of iterations
#' @return convergence is TRUE or FALSE, if TRUE means the reported design is converged
#' @return deriv.ans is the derivative from step 6 of constrained lift-one algorithm
#' @return gmax is the maximum g function in step 8 of constrained lift-one algorithm
#' @return reason is the lift-one loops break reason, either "all derivatives <=0" or "gmax <=0"
#' @export
#'
#' @examples
#' #Example 8 of Trauma data example in Yifei, H., Liping, T., Yang, J. (2025)
#' #Constrained D-optimal design for paid research study
#'
#' J = 5    # number of categories,  >= 3
#' p = 12    # number of parameters
#' m = 8    # number of design points
#' nsample=600 #collect 600 samples finally from the 802 subjects
#' lower.bound <- function(i, w0){
#' n = 600
#' constraint = c(392,410)
#' if(i <= 4){
#'   a.lower <- (sum(w0[5:8])-(constraint[2]/n)*(1-w0[i]))/(sum(w0[5:8]))
#' }
#' else{
#'   a.lower <- (sum(w0[1:4])-(constraint[1]/n)*(1-w0[i]))/(sum(w0[1:4]))
#' }
#' a.lower
#' }

#' upper.bound <- function(i, w0){
#'   n = 600
#'   constraint = c(392,410)
#'   if(i <= 4){
#'     b.upper <- ((constraint[1]/n)*(1-w0[i]) - (sum(w0[1:4])-w0[i]))/(1-sum(w0[1:4]))
#'   }
#'   else{
#'     b.upper <- ((constraint[2]/n)*(1-w0[i]) - (sum(w0[5:8])-w0[i]))/(1-sum(w0[5:8]))
#'   }
#'   b.upper
#' }
#'
#'
#' constraint = c(392,410)
#' g.con = matrix(0,nrow=length(constraint)+1+m, ncol=m)
#' g.con[2:3,] = matrix(data=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1), ncol = m, byrow=TRUE)
#' g.con[1,] = rep(1, m)
#' g.con[4:(length(constraint)+1+m), ] = diag(1, nrow=m)
#' g.dir = c("==", "<=","<=", rep(">=",m))
#' g.rhs = c(1, ifelse((constraint/nsample<1),constraint/nsample,1), rep(0, m))
#' Xi=rep(0,J*p*m)
#' dim(Xi)=c(J,p,m)
#' Xi[,,1] = rbind(c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
#'                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,2] = rbind(c( 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,3] = rbind(c( 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,4] = rbind(c( 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0),
#'               c( 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,5] = rbind(c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,6] = rbind(c( 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,7] = rbind(c( 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' Xi[,,8] = rbind(c( 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 1),
#'                 c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#' thetavec = c(-4.3050, -0.0744,  4.3053, -2.3334, -0.3290, 3.4773,
#' -0.1675, -0.3609, 2.7358, 1.2935, -0.1612, 1.4899)

#' set.seed(123)
#' liftone_constrained_MLM(m=m, p=p, Xi=Xi, J=J, beta=thetavec, lower.bound=lower.bound,
#' upper.bound=upper.bound, g.con=g.con,g.dir=g.dir, g.rhs=g.rhs, w00=NULL,
#' link='cumulative', Fi.func = Fi_func_MLM, reltol=1e-5, maxit=500,
#' delta = 1e-6, epsilon=1e-8, random=TRUE, nram=1)
#'






liftone_constrained_MLM <- function(m, p, Xi, J, beta, lower.bound, upper.bound, g.con, g.dir, g.rhs, w00=NULL, link='cumulative', Fi.func = Fi_func_MLM, reltol=1e-5, maxit=500, delta = 1e-6, epsilon=1e-8, random=TRUE, nram=3) {
  # inputs: m -- number of design points
  #         p --number of parameters in the model
  #         Xi-- J*p*m 3D array of predictors for separate response category at all design points(input to determine ppo,npo,po)
  #         J-- total number of response categories
  #         beta -- p*1 vector of parameters for separate and common response category(input to determine ppo,npo,po)
  #         w00 -- specified initial design proportion
  #         nsample -- sample size
  #         link -- link function of MLM, chose from "continuation", "cumulative", "adjacent", and "baseline"
  #         Fi.func -- function for calculating Fisher information at a specific design point
  #         g.con -- matrix of numeric constraint coefficients, one row per constraint, on column per variable (to be used in lp as const.mat)
  #         g.dir -- Vector of character strings giving the direction of the constraint: each value should be one of "<," "<=," "=," "==," ">," or ">=". (In each pair the two values are identical.) to be used as const.dir in lp()
  #         g.rhs -- 	Vector of numeric values for the right-hand sides of the constraints. to be used as const.rhs in lp()
  #         reltol -- relative convergence tolerance, default value 1e-5
  #         maxit -- maximum number of iterations, default value 100
  #         delta -- a very small number, used in defining alpha_star
  #         epsilon -- a very small number, mainly for comparison like >0, <0, ==0, to reduce errors
  #         random -- generate random initial points
  #         nram -- if random=T, nram is the number of random initial points
  # output: w0 -- original proportion (old p0)
  #         w -- D-optimal design proportion (old p00)
  #         Maximum -- maximized |F| value
  #         convergence -- "T" indicates success, if given "F", consider increasing "maxit"
  #         itmax -- number of iterations used to achieve optimal design
  # NOTE: p is a number, and p0 and w00 are vectors!!!!!

  # n=sum(ni)
  ## step 1: calculate F_i for i=1,...,m, Fisher information F = sum_i^m n_i*F_i, in Theorem 2.1
  # aji <- rep(0, J*m);  dim(aji)=c(J,m)     # (eta_1, eta_2, ..., eta_m)
  # for(i in 1:m) { aji[,i]=Xi[,,i]%*%beta }
  # pji <- rep(0, J*m);  dim(pji)=c(J,m)     # pi_ij=pji[j,i] for logit(pi_ij) = Xi*theta
  # gji <- rep(0, J*m);  dim(gji)=c(J,m)     # gamma_ij=gji[j,i]=pi_i1 + ... + pi_ij
  # for(i in 1:m) {
  #   pji[1,i]=exp(aji[1,i])/(1+exp(aji[1,i]))   # logit^{-1}(eta_ij), j=1,...,J-1
  #   pji[J,i]=1/(1+exp(aji[J-1,i]))
  #   for(j in 2:(J-1)) {
  #     pji[j,i]=exp(aji[j,i])/(1+exp(aji[j,i]))-exp(aji[j-1,i])/(1+exp(aji[j-1,i]))}
  #   gji[,i]= cumsum(pji[,i])
  # }
  # invm <- rep(0, J*J*m);  dim(invm)=c(J,J,m)
  # derm <- rep(0, J*p*m);  dim(derm)=c(J,p,m)
  # tderm <- rep(0, p*J*m);  dim(tderm)=c(p,J,m)
  # dpi <- rep(0, J*J*m);  dim(dpi)=c(J,J,m)
  Fi <- rep(0, p*p*m);  dim(Fi)=c(p,p,m)    # F_i matrix = Fi[,,i]
  #nFi <- rep(0, p*p*m);  dim(nFi)=c(p,p,m)
  # for(i in 1:m) {
  #   for(j in 1:J) {invm[j,J,i]= pji[j,i]}
  #   for(j in 1:(J-1)) {invm[j,j,i]= gji[j,i]*(1-gji[j,i])}
  #   for(j in 1:(J-1)) {invm[j+1,j,i]= -gji[j,i]*(1-gji[j,i])}
  # }
  for(i in 1:m) {
    # derm[, ,i] =invm[, ,i]%*%Xi[, ,i]
    # tderm[, ,i]=t(derm[, ,i])
    # diag(dpi[1:J,1:J,i])=1/pji[1:J,i]
    # Fi[,,i]=tderm[,,i]%*%dpi[,,i]%*%derm[,,i]
    Fi[,,i]=Fi.func(X_x=Xi[,,i], beta=beta, link=link)$F_x
    #nFi[,,i]=Fi[,,i]*ni[i]/n     # F_i*n_i/n = F_i*p_i
  }
  ## end of Step 1

  #F=apply(nFi,c(1,2),sum)        # F = sum_i n_i*F_i
  #Fdet=det(F)                    # |F| at (p_1, p_2,...,p_m)=(n_1,...,n_m)/n

  ## Step 2: Preparation for using Theorem S.9, calculating B_{J-1}
  Bn1 <- matrix(1, J-1, J-1)     # B_{J-1}^{-1}
  for(j in 2:(J-1)) Bn1[,j]=(1:(J-1))^(j-1);
  Bn1 = solve(Bn1);
  ## end of Step 2
  fdet <- function(p) {   # |F|=|sum_i p_i F_i|, p[1:m], need "Fi"
    atemp=p[1]*Fi[,,1];
    for(i in 2:m) atemp=atemp+p[i]*Fi[,,i];
    det(atemp);
  }  # no need to change

  ## Step 3: Calculate f_i(z) in equation (S.9), no need to change
  fiz <- function(z, p, i) {   # f_i(z), need "fdet"
    p1=p*(1-z)/(1-p[i]);
    p1[i]=z;
    fdet(p1);
  }

  ## Step 4: lift-one iteration
  if(is.null(w00)){
    obj_num=sample(seq(nram*2),1)
    max_TF = sample(c(TRUE,FALSE),1)
    rand_matrix = matrix(data=rnorm(nram*m*2),ncol=m)
    w00 = Rglpk_solve_LP(obj=rand_matrix[obj_num,], mat=g.con, dir=g.dir, rhs=g.rhs, max=max_TF)$solution
  } #if no initial design w00, then generate a random initial design w00
  maximum = fdet(w00);
  maxvec = rexp(m);
  convergence = F;
  p0 = w00;
  ind = 0;
  iter = 0;
  #start iteration, two break condition, 1:all derivative <= 0(step7 in algorithm5)  2:max g(w) <=0 (step8 in algorithm5) #updated on 4/09/2022
  repeat{
    iter = iter + 1 #updated 5/09/2022
    #cat(iter,"th, w00=", w00, "; maximum=",maximum, "; p0=", p0, "\n") #updated 5/09/2022
    while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
      io = sample(seq(1:m));
      for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
        avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
        avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
        cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
        for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
        avec[J:2]=Bn1%*%cvec;
        ftemp <- function(z) {   # -f_i(z)
          -(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
        }
        ftemp1 <- function(z) {  # -f'_i(z)
          #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
          -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
        }
        ### determine a.low and b.upper such that z in [a.lower, b.upper] for f_i(z), i=io[ia]
        # a.lower <- lower.bound(io[ia], p0, nsample, constraint);
        # b.upper <- upper.bound(io[ia], p0, nsample, constraint);
        a.lower <- lower.bound(io[ia], p0);
        b.upper <- upper.bound(io[ia], p0);
        initial <- (a.lower+b.upper)/2
        temp=optim(par=initial, fn=ftemp, gr=ftemp1, method="L-BFGS-B", control = list(factr=0, maxit=1e4, pgtol=0, lmm=8, ndeps=1e-1), lower=max(0,a.lower), upper=min(1,b.upper));
        #temp=optim(par=initial, fn=ftemp, gr=ftemp1, method="Brent", control = list(reltol=1e-20, maxit=1e4), lower=max(0,a.lower), upper=min(1,b.upper));
        zstar=temp$par;         # z_*
        fstar=-temp$value;
        ### need to check if 0 in [a.lower, b.upper]
        if((fstar <= avec[1]) & (0 > a.lower) & (0 < b.upper)){zstar=0; fstar=avec[1];}; ### need to check 0
        ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
        ptemp1[io[ia]] = zstar;
        if(fstar > maximum) {maximum = fstar; p0=ptemp1;}; #updated on 4/24/2022, derivative calculation after finding w*
        maxvec[io[ia]] = maximum;
      }
      ind = ind+1;
    }# end of "while"
    p0.ans=p0; maximum.ans=maximum;w00.ans = w00;

    #calculate derivatives #updated 4/24/2022
    derivative = rep(NA, m)
    for(i in 1:m){
      derivative[io[i]] = ftemp1(p0.ans[io[i]])
    }

    # check if all the derivatives smaller than 0 updated 4/02/2022
    po = NA
    alpha_star = NA
    if(sum(derivative>epsilon)<=0){ #updated on 5/09/2022 if fi'(w_i*) <= 0 break
      #cat("\nall derivative <=0, break\n")
      reason = "all derivative <=0"
      deriv.ans = derivative #updated 4/24/2022
      gmax = NA #updated 5/05/2022
      break;
    }
    else{ #updated on 4/09/2022
      g.obj = (1-p0.ans)*derivative
      # if(is.null(g.con)){
      #   g.con = matrix(0,nrow=length(constraint)+1+m, ncol=m)
      #   g.con[2:3,] = matrix(data=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1), ncol = m, byrow=TRUE)
      #   g.con[1,] = rep(1, m)
      #   g.con[4:(length(constraint)+1+m), ] = diag(1, nrow=m)
      # }
      # if(is.null(g.dir)){
      #   g.dir = c("==", "<=","<=", rep(">=",m))
      # }
      # if(is.null(g.rhs)){
      #   g.rhs = c(1, ifelse((constraint/nsample<1),constraint/nsample,1), rep(0, m))
      # }
      solution=lp("max", g.obj, g.con, g.dir, g.rhs)
      po = solution$solution
      gmax = sum(g.obj%*%po) #updated on 5/09/2022
      if(gmax <= epsilon){ #updated on 4/26/2022 change 0 to epsilon
        #cat("\ngmax <=0, break\n")
        reason = "gmax <=0"
        deriv.ans = derivative #updated 4/24/2022
        break;
      }
      else{
        #Bp matrix calculation  #updated 4/02/2022
        B = matrix(NA, nrow=p, ncol=p)
        for(s in 1:p){
          for(t in 1:p){
            B[s,t] = (s/p)^t
          }
        }

        #f(w) function f(w) = det(sum(wi*F_i)) for i from 1 to m  #updated 4/02/2022
        fwfunction = function(wi){
          fw <- rep(0, p*p*m)
          dim(fw)=c(p,p,m)
          for(i in 1:m) {
            fw[,,i]=Fi[,,i]*wi[i]   # F_i*w_i
          }
          Fw=apply(fw,c(1,2),sum)        # F = sum_i n_i*F_i
          Fwdet=det(Fw) # |F| at (p_1, p_2,...,p_m)=(n_1,...,n_m)/n
          return(Fwdet)
        }

        #h(1/p), ..., h((p-1)/p) calculation #updated 4/02/2022
        h_const = rep(0, p)
        for(i in 1:p){
          wi = (1-(i/p))*p0.ans + (i/p)*po
          h_const[i] = fwfunction(wi)
        }
        h_const = h_const - fwfunction(p0.ans)
        #constant in h function in 9th step
        C = solve(B)%*% h_const
        #hprime function #updated 4/02/2022
        hprime = function(alpha){
          param = seq(1:p)
          A = param * alpha^(param-1)
          h_prime = t(C) %*% A
          return(h_prime)
        }
        #h function
        h = function(alpha){ #updated 4/25/2022
          param = seq(1:p)
          A = alpha^(param)
          h_result = t(C) %*% A + fwfunction(p0.ans)
          return(h_result)
        }

        #updated on 5/09/2022
        alpha_star = 1-delta
        #cat("hprime(1):", hprime(1), ", hprime(0):", hprime(0), ", gmax:", gmax)
        if((hprime(1) >=0)&(h(1) > 0)){alpha_star = 1};
        if((h(1) > epsilon)&(hprime(1) < -epsilon)){alpha_star = uniroot(hprime,c(0,1))$root};
        if((hprime(1)< -epsilon)  & (h(1) < epsilon)){alpha_star = uniroot(hprime,c(0,1))$root};

        #set new starting point and repeat lift-one #updated on 4/02/2022
        p0 = (1-alpha_star)*p0.ans + alpha_star * po
        maximum = fdet(p0);
        convergence = F;
      }

    }
  }
  w00.original = w00.ans
  p0.report=p0.ans
  maximum.report=maximum.ans
  derivative.report = derivative
  po.report = po
  gmax.report = gmax
  alpha_star.report = alpha_star
  reason.report = reason

  #if random initial points is needed, run the following {...} #updated 5/12/2022
  if(random==T){
    for(u in 1:nram){
      #generate random starting points that satisfy constraints
      #in trauma example, first 4 <= constraint[1]/nsample, last 4 <= constraint[2]/nsample, sum = 1
      # w00 = rep(0, m)
      # for(ui in 1:m){
      #   if(ui <= 4){
      #     w00[ui] = round(runif(1, min = 0, max = min(constraint[1]-sum(w00[1:4]), nsample-sum(w00[1:8])) ))
      #   }else{
      #     w00[ui] = round(runif(1, min=0, max = min(constraint[2]-sum(w00[5:8]), nsample-sum(w00[1:8])) ))
      #   }
      # }
      # w00 = w00/sum(w00)
      # #updated 5/29/2022 in case sum(w00) != 600
      # if(sum(w00[1:4])>constraint[1]/nsample){w00[5:8] = w00[5:8]+(sum(w00[1:4])-constraint[1]/nsample)/4 ; if(w00[1]-(sum(w00[1:4])-constraint[1]/nsample)>0){w00[1]=w00[1]-(sum(w00[1:4])-constraint[1]/nsample)}else{w00[1:2]=w00[1:2]-(sum(w00[1:4])-constraint[1]/nsample)/2}}
      # if(sum(w00[5:8])>constraint[2]/nsample){w00[1:4] = w00[1:4]+(sum(w00[5:8])-constraint[2]/nsample)/4 ; if(w00[5]-(sum(w00[5:8])-constraint[2]/nsample)>0){w00[5]=w00[5]-(sum(w00[5:8])-constraint[2]/nsample)/2 }else{w00[5:6]=w00[5:6]-(sum(w00[5:8])-constraint[2]/nsample)/2 }}

      obj_num=sample(seq(nram*2),1)
      max_TF = sample(c(TRUE,FALSE),1)
      rand_matrix = matrix(data=rnorm(nram*m*2),ncol=m)
      w00 = Rglpk_solve_LP(obj=rand_matrix[obj_num,], mat=g.con, dir=g.dir, rhs=g.rhs, max=max_TF)$solution
      #if no initial design w00, then generate a random initial design w00
      maximum = fdet(w00);
      maxvec = rexp(m);
      p0 = w00

      repeat{
        iter = iter + 1 #updated 5/09/2022
        #cat(iter,"th, w00=", w00, "; maximum=",maximum, "; p0=", p0, "\n") #updated 5/09/2022
        while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
          io = sample(seq(1:m));
          for(ia in 1:m) {    # run updating in random order of {1,2,...,m}
            avec <- rep(0, J);    # a0, a1, ..., a_{J-1}
            avec[1] = fiz(0, p0, io[ia]);  # a0=f_i(0)
            cvec <- rep(0, J-1);  # c1, c2, ..., c_{J-1}
            for(j in 1:(J-1)) cvec[j]=(j+1)^p*j^(J-1-p)*fiz(1/(j+1), p0, io[ia])-j^(J-1)*avec[1];
            avec[J:2]=Bn1%*%cvec;
            ftemp <- function(z) {   # -f_i(z)
              -(1-z)^(p-J+1)*sum(avec*z^(0:(J-1))*(1-z)^((J-1):0));
            }
            ftemp1 <- function(z) {  # -f'_i(z)
              #   -(1-z)^(p-J)*sum(((1:(J-1))-P*z)*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+p*avec[1]*(1-z)^(p-1);
              -(1-z)^(p-J+1)*sum((1:(J-1))*avec[2:J]*z^(0:(J-2))*(1-z)^((J-2):0))+(1-z)^(p-J)*sum((p:(p-J+1))*avec[1:J]*z^(0:(J-1))*(1-z)^((J-1):0));
            }
            ### determine a.low and b.upper such that z in [a.lower, b.upper] for f_i(z), i=io[ia]
            # a.lower <- lower.bound(io[ia], p0, nsample, constraint);
            # b.upper <- upper.bound(io[ia], p0, nsample, constraint);
            a.lower <- lower.bound(io[ia], p0)
            b.upper <- upper.bound(io[ia], p0)
            initial <- (a.lower+b.upper)/2
            temp=optim(par=initial, fn=ftemp, gr=ftemp1, method="L-BFGS-B", control = list(factr=0, maxit=1e4, pgtol=0, lmm=8, ndeps=1e-1), lower=max(0,a.lower), upper=min(1,b.upper));
            #temp=optim(par=initial, fn=ftemp, gr=ftemp1, method="Brent", control = list(reltol=1e-20, maxit=1e4), lower=max(0,a.lower), upper=min(1,b.upper));
            zstar=temp$par;         # z_*
            fstar=-temp$value;
            ### need to check if 0 in [a.lower, b.upper]
            if((fstar <= avec[1]) & (0 > a.lower) & (0 < b.upper)){zstar=0; fstar=avec[1];}; ### need to check 0
            ptemp1 = p0*(1-zstar)/(1-p0[io[ia]]);
            ptemp1[io[ia]] = zstar;
            if(fstar > maximum) {maximum = fstar; p0=ptemp1;}; #updated on 4/24/2022, derivative calculation after finding w*
            maxvec[io[ia]] = maximum;
          }
          ind = ind+1;
        }# end of "while"
        p0.ans=p0; maximum.ans=maximum;w00.ans = w00;

        #calculate derivatives #updated 4/24/2022
        derivative = rep(NA, m)
        for(i in 1:m){
          derivative[io[i]] = ftemp1(p0.ans[io[i]])
        }

        # check if all the derivatives smaller than 0 updated 4/02/2022
        po = NA
        alpha_star = NA
        ####can lift restriction####
        if(sum(derivative>epsilon)<=0){ #updated on 5/09/2022 if fi'(w_i*) <= 0 break
          #cat("\nall derivative <=0, break\n")
          reason = "all derivative <=0"
          deriv.ans = derivative #updated 4/24/2022
          gmax = NA #updated 5/05/2022
          break;
        }
        else{ #updated on 4/09/2022
          g.obj = (1-p0.ans)*derivative
          # if(is.null(g.con)){
          #   g.con = matrix(0,nrow=length(constraint)+1+m, ncol=m)
          #   g.con[2:3,] = matrix(data=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1), ncol = m, byrow=TRUE)
          #   g.con[1,] = rep(1, m)
          #   g.con[4:(length(constraint)+1+m), ] = diag(1, nrow=m)
          # }
          # if(is.null(g.dir)){
          #   g.dir = c("==", "<=","<=", rep(">=",m))
          # }
          # if(is.null(g.rhs)){
          #   g.rhs = c(1, ifelse((constraint/nsample<1),constraint/nsample,1), rep(0, m))
          # }
          solution=lp("max", g.obj, g.con, g.dir, g.rhs)
          po = solution$solution
          gmax = sum(g.obj%*%po) #updated on 5/09/2022
          if(gmax <= epsilon){ #updated on 4/26/2022 change 0 to epsilon
            #cat("\ngmax <=0, break\n")
            reason = "gmax <=0"
            deriv.ans = derivative #updated 4/24/2022
            break;
          }
          else{
            #Bp matrix calculation  #updated 4/02/2022
            B = matrix(NA, nrow=p, ncol=p)
            for(s in 1:p){
              for(t in 1:p){
                B[s,t] = (s/p)^t
              }
            }

            #f(w) function f(w) = det(sum(wi*F_i)) for i from 1 to m  #updated 4/02/2022
            fwfunction = function(wi){
              fw <- rep(0, p*p*m)
              dim(fw)=c(p,p,m)
              for(i in 1:m) {
                fw[,,i]=Fi[,,i]*wi[i]   # F_i*w_i
              }
              Fw=apply(fw,c(1,2),sum)        # F = sum_i n_i*F_i
              Fwdet=det(Fw) # |F| at (p_1, p_2,...,p_m)=(n_1,...,n_m)/n
              return(Fwdet)
            }

            #h(1/p), ..., h((p-1)/p) calculation #updated 4/02/2022
            h_const = rep(0, p)
            for(i in 1:p){
              wi = (1-(i/p))*p0.ans + (i/p)*po
              h_const[i] = fwfunction(wi)
            }
            h_const = h_const - fwfunction(p0.ans)
            #constant in h function in 9th step
            C = solve(B)%*% h_const
            #hprime function #updated 4/02/2022
            hprime = function(alpha){
              param = seq(1:p)
              A = param * alpha^(param-1)
              h_prime = t(C) %*% A
              return(h_prime)
            }
            #h function
            h = function(alpha){ #updated 4/25/2022
              param = seq(1:p)
              A = alpha^(param)
              h_result = t(C) %*% A + fwfunction(p0.ans)
              return(h_result)
            }

            #updated on 5/09/2022
            alpha_star = 1-delta
            #cat("\nhprime(1):", hprime(1), ", hprime(0):", hprime(0), ", gmax:", gmax, ", po:", po, ", derivative:", derivative, ", g.obj", g.obj, ", p0:", p0.ans, ", w00:", w00.ans)
            if((hprime(1) >=0)&(h(1) > 0)){alpha_star = 1};
            if((h(1) > epsilon)&(hprime(1) < -epsilon)){alpha_star = uniroot(hprime,c(0,1))$root};
            if((hprime(1)< -epsilon)  & (h(1) < epsilon)){alpha_star = uniroot(hprime,c(0,1))$root};


            #set new starting point and repeat lift-one #updated on 4/02/2022
            p0 = (1-alpha_star)*p0.ans + alpha_star * po
            maximum = fdet(p0);
            convergence = F;
          }

        }
      }

      if(maximum.ans > maximum.report) {
        maximum.report=maximum.ans;
        p0.report=p0.ans;
        convergence=F;
        if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=T;
        w00.original=w00.ans;
        derivative.report = derivative
        po.report = po
        gmax.report = gmax
        alpha_star.report = alpha_star
        reason.report = reason
        itmax=ind;
      }
    } #end "for" loop of u in 1:nram
  }# end of random initial

  #cat(iter,"th, p0=", p0, "; maximum=",maximum, "\n") #updated 5/09/2022
  #maximum.adj=maximum.report*n^p;
  #fdet.adj=Fdet*n^p;
  if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=T;
  itmax=ind;
  #effi=(Fdet/maximum.report)^(1/p)
  list(w=p0.report, w0=w00.original, Maximum=maximum.report, iteration=ind, convergence=convergence, deriv.ans = derivative.report, gmax=gmax.report, reason = reason.report); # updated on 4/02/2022
}

