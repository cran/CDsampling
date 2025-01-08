#' Find constrained D-optimal approximate design for generalized linear models (GLM)
#' @importFrom stats rexp
#' @importFrom stats uniroot
#' @importFrom stats runif
#'
#' @param X Model matrix, with nrow = num of design points and ncol = num of parameters
#' @param W Diagonal of W matrix in Fisher information matrix, can be calculated from W_func_GLM in package
#' @param g.con A matrix of numeric constraint coefficients, one row per constraint, on column per variable (to be used in as const.mat lp() and mat in Rglpk_solve_LP())
#' @param g.dir Vector of character strings giving the direction of the constraint: each value should be one of "<," "<=," "=," "==," ">," or ">=". (In each pair the two values are identical.) to be used as const.dir in lp() and dir in Rglpk_solve_LP()
#' @param g.rhs Vector of numeric values for the right-hand sides of the constraints. to be used as const.rhs in lp() and rhs in Rglpk_solve_LP().
#' @param lower.bound A function to determine lower bound r_i1 in Step 3 of Constrained lift-one algorithm from Yifei, H., Liping, T., Yang, J. (2023) Constrained D-optimal design for paid research study
#' @param upper.bound A function to determine upper bound r_i2 in Step 3 of Constrained lift-one algorithm from Yifei, H., Liping, T., Yang, J. (2023) Constrained D-optimal design for paid research study
#' @param reltol The relative convergence tolerance, default value 1e-5
#' @param maxit The maximum number of iterations, default value 500
#' @param random TRUE or FALSE, if TRUE then the function will run with additional "nram" number of random initial points, default to be TRUE
#' @param nram When random == TRUE, the function will generate nram number of initial points, default is 3
#' @param w00 Specified initial design proportion; default to be NULL, this will generate a random initial design
#' @param epsilon A very small number, for comparison of >0, <0, ==0, to reduce errors, default 1e-12
#'
#' @return w is the approximate D-optimal design
#' @return w0 is the initial design used to get optimal design w
#' @return maximum is the maximized |F| value
#' @return itmax is the number of iterations
#' @return convergence is TRUE or FALSE, if TRUE means the reported design is converged
#' @return deriv.ans is the derivative from step 6 of constrained lift-one algorithm
#' @return gmax is the maximum g function in step 8 of constrained lift-one algorithm
#' @return reason is the lift-one loops break reason, either "all derivatives <=0" or "gmax <=0"
#' @export
#'
#' @examples
#' #Example 6 in Section 3.4 of Yifei, H., Liping, T., Yang, J. (2025)
#' #Constrained D-optimal design for paid research study
#'
#' #main effect model beta_0, beta_1, beta_21, beta_22
#' beta = c(0, -0.1, -0.5, -2)
#'
#' #gives the 6 categories (0,0,0), (0,1,0),(0,0,1),(1,0,0),(1,1,0),(1,0,1)
#' X.liftone=matrix(data=c(1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1),
#' ncol=4, byrow=TRUE)
#'
#' #calculate W matrix based on beta's under logit link
#' W_matrix=W_func_GLM(X= X.liftone, b=beta)
#'
#' m=6 #number of categories
#' nsample = 200
#' rc = c(50, 40, 10, 200, 150, 50)/nsample
#' g.con = matrix(0,nrow=(2*m+1), ncol=m)
#' g.con[1,] = rep(1, m)
#' g.con[2:(m+1),] = diag(m)
#' g.con[(m+2):(2*m+1), ] = diag(m)
#' g.dir = c("==", rep("<=", m), rep(">=", m))
#' g.rhs = c(1, rc, rep(0, m))
#'
#' lower.bound=function(i, w){
#'   nsample = 200
#'   rc = c(50, 40, 10, 200, 150, 50)/nsample
#'   m=length(w) #num of categories
#'   temp = rep(0,m)
#'   temp[w>0]=1-pmin(1,rc[w>0])*(1-w[i])/w[w>0];
#'   temp[i]=0;
#'   max(0,temp);
#' }

#' upper.bound=function(i, w){
#'   nsample = 200
#'   rc = c(50, 40, 10, 200, 150, 50)/nsample
#'   m=length(w) #num of categories
#'   rc[i];
#'   min(1,rc[i]);
#' }
#'
#' approximate_design = liftone_constrained_GLM(X=X.liftone, W=W_matrix,
#' g.con=g.con, g.dir=g.dir, g.rhs=g.rhs, lower.bound=lower.bound,
#' upper.bound=upper.bound, reltol=1e-10, maxit=100, random=TRUE, nram=4,
#' w00=NULL, epsilon = 1e-8)
#'
#'
#'

liftone_constrained_GLM <- function(X, W, g.con, g.dir, g.rhs, lower.bound, upper.bound, reltol=1e-5, maxit=500, random=TRUE, nram=3, w00=NULL, epsilon = 1e-12){
  # W=W[1,2,...,m] are strictly positive; diagonal of W matrix in Fisher information matrix
  # rc=rc[1,2,...,m] are strictly positive; each group's proportion (n_group / n_total)
  # if random=T, run 5 random initial points and pick up the best; default initial p1=p2=...=1/m
  # epsilon -- a very small number, mainly for comparison in the algorithm e.g. > 0, < 0, ==0, to avoid errors
  # output: w=w--optimal design based on "det"
  #         Maximum--maximized value of "det"
  #         convergence -- "T" indicates success
  #         w0 -- initial w
  #         itmax -- number of iterations
  m = dim(X)[1]; #num of categories
  d = dim(X)[2]; #num of parameters
  W = W[1:m];
  #rc = pmin(1, rc[1:m]);
  ftemp <- function(p) { det(t(X * (p*W)) %*% X);};
  fiprime <- function(a,b,x) {((1-x)^(d-2))*(a-b*d+(b-a)*d*x); }; #updated on 4/20/2022
  #generate random initial approximate allocations #updated on 03/02/2024
  if(is.null(w00)){
    obj_num=sample(seq(m),1)
    max_TF = sample(c(T,F),1)
    w00 = Rglpk::Rglpk_solve_LP(obj=diag(m)[obj_num,], mat=g.con, dir=g.dir, rhs=g.rhs, max=max_TF)$solution
  }

  # if(is.null(w00)) w00=rep(1/m,m);
  # w00=pmin(w00,rc);             # updated on 7/16/2021
  # if(sum(w00)<1) w00 = w00 + (rc-w00)*((1-sum(w00))/sum(rc-w00));
  maximum = ftemp(w00);
  maxvec = rexp(m);     # initial values
  convergence=FALSE;
  p = w00;
  ind = 0;

  #liftone iteration: two break condition, 1:all derivative <= 0(step7 in algorithm5)  2:max g(w) <=0 (step8 in algorithm5) #updated on 4/09/2022
  repeat{
    while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
      io = sample(seq(1:m));
      for(i in 1:m) {    # run updating in random order of p
        if(p[io[i]]>0) { # if p_i > 0
          ptemp1 = p/(1-p[io[i]]);
          ptemp1[io[i]] = 0;
          b = ftemp(ptemp1);      # b=fs(0)
          a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
        } else {         # p_i=0
          b = maximum;
          ptemp1 = p/2;
          ptemp1[io[i]] = 1/2;          # for fs(1/2)
          a = ftemp(ptemp1)*2^d - b;
        }
        # temp=rep(0,m);
        # temp[p>0]=1-pmin(1,rc[p>0])*(1-p[io[i]])/p[p>0];
        # temp[io[i]]=0;         # updated on 7/13/2021
        # z1 <- max(0,temp);
        # z2 <- min(1, rc[io[i]]);
        #update 03/02/2024 define new boundary
        z1 = max(0, lower.bound(io[i], p))
        z2 = min(1, upper.bound(io[i], p))

        if(a > b*d) {
          z0=(a-b*d)/((a-b)*d);
          if(z0<z1) x=z1 else {if(z0>z2) x=z2 else x=z0;};
        } else x=z1;
        ptemp1 = p*(1-x)/(1-p[io[i]]);   # new p
        ptemp1[io[i]] = x;
        temp=ftemp(ptemp1);              # updated on 7/13/2021
        if(temp>maximum) {maximum = temp; p = ptemp1;};  # updated on 4/19/2022
        maxvec[io[i]] = maximum;
      }
      ind = ind+1;
    }
    p.ans=p; maximum.ans=maximum; p0.ans=w00; if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;itmax=ind; # updated on 4/19/2022

    #calculate the derivatives, updated on 4/19/2022
    derivative = rep(NA, m);
    p=p.ans;
    maximum=maximum.ans;
    for(i in 1:m) {    # run updating in random order of weight
      if(p[io[i]]>0) {
        ptemp1 = p/(1-p[io[i]]);
        ptemp1[io[i]] = 0;
        b = ftemp(ptemp1);      # b=fs(0)
        a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
      } else {         # p[io[i]]=0
        b = maximum;
        ptemp1 = p/2;
        ptemp1[io[i]] = 1/2;          # for fs(1/2)
        a = ftemp(ptemp1)*2^d - b;
      }
      derivative[io[i]]=fiprime(a,b,p[io[i]]);
    }
    #check the derivative, if all derivative <=0, break, ow continue;
    if(sum(derivative>epsilon)<=0){ #updated 5/05/2022 change 0 to epsilon, f'_i(w_i*) <= 0, break
      reason.ans = "all derivative <= 0"
      deriv.ans=derivative        # updated on 4/19/2022
      po = NA #updated 5/07/2022
      gmax = NA #updated 5/07/2022
      break;
    }
    else{ #updated on 4/09/2022
      g.obj = (1-p.ans)*derivative
      #g.con = matrix(0,nrow=(2*m+1), ncol=m)
      #g.con[1,] = rep(1, m)
      #g.con[2:(m+1),] = diag(m)
      #g.con[(m+2):(2*m+1), ] = diag(m)
      #g.dir = c("=", rep("<=", m), rep(">=", m))
      #g.rhs = c(1, rc, rep(0, m))
      solution=lpSolve::lp(direction="max", objective.in=g.obj, const.mat=g.con, const.dir=g.dir, const.rhs=g.rhs)
      po = solution$solution
      gmax = solution$objval
      if(gmax <= epsilon){#updated on 4/26/2022 change 0 to epsilon, g(wo) <= 0, break
        reason.ans = "gmax <= 0"
        deriv.ans=derivative        # updated on 4/19/2022
        break;
      }
      else{
        #Bp matrix calculation  #updated 4/02/2022
        B = matrix(NA, nrow=d, ncol=d)
        for(s in 1:d){
          for(t in 1:d){
            B[s,t] = (s/d)^t
          }
        }
        h_const = rep(0, d)
        for(i in 1:d){
          wi = (1-(i/d))*p.ans + (i/d)*po
          h_const[i] = ftemp(wi)
        }
        h_const = h_const - ftemp(p.ans)
        #constant in h function in 9th step
        C = solve(B)%*% h_const
        #hprime function #updated 4/02/2022
        hprime = function(alpha){
          param = seq(1:d)
          A = param * alpha^(param-1)
          h_prime = t(C) %*% A
          return(h_prime)
        }
        h = function(alpha){ #updated 4/25/2022
          param = seq(1:d)
          A = alpha^(param)
          h_result = t(C) %*% A + ftemp(p.ans)
          return(h_result)
        }
        #theorem 4.6 alpha_star determined based on hprime(1) and h(1) #updated 4/26/2022
        if((hprime(1) >= -epsilon) & (h(1) > epsilon)){alpha_star = 1};
        if(((hprime(1) < -epsilon) & (h(1) > epsilon)) | (abs(h(1)) <= epsilon)){alpha_star = uniroot(hprime,c(0,1-epsilon))$root}; #updated 5/05/2022
        if((ftemp(po) < ftemp(p.ans)) & (alpha_star==1)){warning("f(wo) < f(w*) and alpha*=1")} #updated 4/24/2022 check if f(po) < f(p.ans) and alpha*=1 report it
        #set new starting point and repeat lift-one #updated on 4/02/2022
        p = (1-alpha_star)*p.ans + alpha_star * po
        convergence=FALSE;
        maximum = ftemp(p)
      }
    }
  }
  p.report=p.ans; maximum.report=maximum.ans; p0.report = w00; deriv.report = deriv.ans; po.report = po; gmax.report = gmax; reason.report = reason.ans

  if(random==T) for(j in 1:nram) {
    #p0=runif(m, min=0.5, max=500); p0=p0/sum(p0);
    #p0=pmin(p0, rc);             # updated on 7/16/2021
    #if(sum(p0)<1) p0 = p0 + (rc-p0)*((1-sum(p0))/sum(rc-p0));
    #generate random initial allocations #updated on 03/02/2024
    obj_num=sample(seq(m),1)
    max_TF = sample(c(T,F),1)
    w00 = Rglpk::Rglpk_solve_LP(obj=diag(m)[obj_num,], mat=g.con, dir=g.dir, rhs=g.rhs, max=max_TF)$solution

    p=w00;
    maxvec = runif(m, min=0.5, max=500);
    maximum = ftemp(p);
    ind = 0;

    repeat{
      while((ind < maxit) && ((max(maxvec)/min(maxvec))-1 > reltol)) {
        io = sample(seq(1:m));
        for(i in 1:m) {    # run updating in random order of weight
          if(p[io[i]]>0) {
            ptemp1 = p/(1-p[io[i]]);
            ptemp1[io[i]] = 0;
            b = ftemp(ptemp1);      # b=fs(0)
            a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
          } else {         # p[io[i]]=0
            b = maximum;
            ptemp1 = p/2;
            ptemp1[io[i]] = 1/2;          # for fs(1/2)
            a = ftemp(ptemp1)*2^d - b;
          }
          # temp=rep(0,m);
          # temp[p>0]=1-pmin(1,rc[p>0])*(1-p[io[i]])/p[p>0];
          # temp[io[i]]=0;            # updated on 7/13/2021
          # z1 <- max(0,temp);
          # z2 <- min(1, rc[io[i]]);
          #update 03/02/2024 define new boundary
          z1 = max(0, lower.bound(io[i], p))
          z2 = min(1, upper.bound(io[i], p))
          if(a > b*d) {
            z0=(a-b*d)/((a-b)*d);
            if(z0<z1) x=z1 else {if(z0>z2) x=z2 else x=z0;};
          } else x=z1;
          ptemp1 = p*(1-x)/(1-p[io[i]]);   # new p
          ptemp1[io[i]] = x;
          temp=ftemp(ptemp1);       # updated on 7/13/2021
          if(temp>maximum) {maximum = temp; p = ptemp1;}; # updated on 4/19/2022
          maxvec[io[i]] = maximum;
        }
        ind = ind+1;
      }
      p.ans=p; maximum.ans=maximum; p0.ans=w00; if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;itmax=ind; # updated on 5/21/2022

      #calculate the derivatives, updated on 4/19/2022
      derivative = rep(NA, m);
      p=p.ans;
      maximum=maximum.ans;
      for(i in 1:m) {    # run updating in random order of w
        if(p[io[i]]>0) {
          ptemp1 = p/(1-p[io[i]]);
          ptemp1[io[i]] = 0;
          b = ftemp(ptemp1);      # b=fs(0)
          a = (maximum - b*(1-p[io[i]])^d)/(p[io[i]]*(1-p[io[i]])^(d-1));
        } else {         # p[io[i]]=0
          b = maximum;
          ptemp1 = p/2;
          ptemp1[io[i]] = 1/2;          # for fs(1/2)
          a = ftemp(ptemp1)*2^d - b;
        }
        derivative[io[i]]=fiprime(a,b,p[io[i]]);
      }
      #check the derivative, if all derivative <=0, break, ow continue;
      if(sum(derivative>epsilon)<=0){ #updated 5/05/2022 change 0 to epsilon, f'_i(w_i*) <= 0, break
        reason.ans = "all derivative <= 0"
        deriv.ans=derivative        # updated on 4/19/2022
        po = NA #updated 5/07/2022
        gmax = NA #updated 5/07/2022
        break;
      }
      else{ #updated on 4/09/2022
        g.obj = (1-p.ans)*derivative
        # g.con = matrix(0,nrow=(2*m+1), ncol=m)
        # g.con[1,] = rep(1, m)
        # g.con[2:(m+1),] = diag(m)
        # g.con[(m+2):(2*m+1), ] = diag(m)
        # g.dir = c("=", rep("<=", m), rep(">=", m))
        # g.rhs = c(1, rc, rep(0, m))
        solution=lpSolve::lp(direction="max", objective.in=g.obj, const.mat=g.con, const.dir=g.dir, const.rhs=g.rhs)
        po = solution$solution
        gmax = solution$objval
        if(gmax <= epsilon){#updated on 4/26/2022 change 0 to epsilon, g(wo) <= 0, break
          reason.ans = "gmax <= 0"
          deriv.ans=derivative        # updated on 4/19/2022
          break;
        }
        else{
          #Bp matrix calculation  #updated 4/02/2022
          B = matrix(NA, nrow=d, ncol=d)
          for(s in 1:d){
            for(t in 1:d){
              B[s,t] = (s/d)^t
            }
          }
          h_const = rep(0, d)
          for(i in 1:d){
            wi = (1-(i/d))*p.ans + (i/d)*po
            h_const[i] = ftemp(wi)
          }
          h_const = h_const - ftemp(p.ans)
          #constant in h function in 9th step
          C = solve(B)%*% h_const
          #hprime function #updated 4/02/2022
          hprime = function(alpha){
            param = seq(1:d)
            A = param * alpha^(param-1)
            h_prime = t(C) %*% A
            return(h_prime)
          }
          h = function(alpha){ #updated 4/25/2022
            param = seq(1:d)
            A = alpha^(param)
            h_result = t(C) %*% A + ftemp(p.ans)
            return(h_result)
          }
          #theorem 4.6 alpha_star determined based on hprime(1) and h(1) #updated 4/26/2022
          if((hprime(1) >= -epsilon) & (h(1) > epsilon)){alpha_star = 1};
          if(((hprime(1) < -epsilon) & (h(1) > epsilon)) | (abs(h(1)) <= epsilon)){alpha_star = uniroot(hprime,c(0,1-epsilon))$root}; #updated 5/05/2022
          if((ftemp(po) < ftemp(p.ans)) & (alpha_star==1)){warning("f(wo)<f(w*) and alpha*=1")} #updated 4/24/2022 check if f(po) < f(p.ans) and alpha*=1 report it
          #set new starting point and repeat lift-one #updated on 4/02/2022
          p = (1-alpha_star)*p.ans + alpha_star * po
          convergence=FALSE;
          maximum = ftemp(p)
        }
      }
    }#end of repeat

    if(maximum.ans > maximum.report) {
      maximum.report=maximum.ans;
      p.report=p.ans;
      deriv.report = deriv.ans
      po.report = po
      gmax.report = gmax
      convergence=FALSE;
      if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;
      p0.report=p0.ans;
      itmax=ind;
      reason.report = reason.ans
    }
  }

  if((max(maxvec)/min(maxvec))-1 <= reltol) convergence=TRUE;itmax=ind;
  #define S3 class
  output <- list(w=p.report, w0=p0.report, maximum=maximum.report, convergence=convergence, itmax=itmax, deriv.ans=deriv.report, gmax=gmax.report, reason=reason.report);
  class(output)<-"list_output"
  return(output)
}




