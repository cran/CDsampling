## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CDsampling)

## -----------------------------------------------------------------------------
beta = c(0.5, 0.5, 0.5) #coefficients
X = matrix(data=c(1,-1,-1,1,-1,1,1,1,-1), byrow=TRUE, nrow=3) #design matrix
w = c(1/3,1/3,1/3) #approximate allocation
CDsampling::F_func_GLM(w=w, beta=beta, X=X, link='logit')

## -----------------------------------------------------------------------------
J=5; p=12; m=8; #response levels; num of parameters; num of design points
beta = c(-4.047, -0.131, 4.214, -2.225, -0.376, 3.519, -0.302,
    -0.237,  2.420, 1.386,  -0.120,  1.284) 
Xi=rep(0,J*p*m) #design matrix   
dim(Xi)=c(J,p,m)
Xi[,,1] = rbind(c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

Xi[,,2] = rbind(c( 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

Xi[,,3] = rbind(c( 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,4] = rbind(c( 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,5] = rbind(c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,6] = rbind(c( 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,7] = rbind(c( 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,8] = rbind(c( 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
 alloc = rep(1/8,m) #approximate allocation
 CDsampling::F_func_MLM(w=alloc, beta=beta, X=Xi, link='cumulative')

## -----------------------------------------------------------------------------
beta = c(0, 3, 3, 3) #coefficients 
#design matrix X
X=matrix(data=c(1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,1), ncol=4, byrow=TRUE) 
nsample=200 #sample size

## -----------------------------------------------------------------------------
W_matrix=CDsampling::W_func_GLM(X=X, b=beta, link="logit") #W matrix

## -----------------------------------------------------------------------------
rc = c(50, 40, 10, 200, 150, 50)/nsample #constraints for each subgroup
m = 6
g.con = matrix(0,nrow=(2*m+1), ncol=m)
g.con[1,] = rep(1, m)
g.con[2:(m+1),] = diag(m)
g.con[(m+2):(2*m+1), ] = diag(m)
g.dir = c("==", rep("<=", m), rep(">=", m))
g.rhs = c(1, rc, rep(0, m))

lower.bound=function(i, w){
  nsample = 200
  rc = c(50, 40, 10, 200, 150, 50)/nsample
  m=length(w) #num of categories
  temp = rep(0,m)
  temp[w>0]=1-pmin(1,rc[w>0])*(1-w[i])/w[w>0];
  temp[i]=0;
  max(0,temp);
}

upper.bound=function(i, w){
  nsample = 200
  rc = c(50, 40, 10, 200, 150, 50)/nsample
  m=length(w) #num of categories
  rc[i];
  min(1,rc[i])
}

## -----------------------------------------------------------------------------
label = c("F, 18-25", "F, 26-64", "F, >=65", "M, 18-25", "M, 26-64", "M, >=65")

## -----------------------------------------------------------------------------
set.seed(2025)
approximate_design = CDsampling::liftone_constrained_GLM(X=X, W=W_matrix, g.con=g.con, g.dir=g.dir, 
                                               g.rhs=g.rhs, lower.bound=lower.bound, 
                                               upper.bound=upper.bound, label=label, reltol=1e-10, 
                                               maxit=100, random=TRUE, nram=4, w00=NULL, 
                                               epsilon=1e-8)

print(approximate_design)

## -----------------------------------------------------------------------------
exact_design = CDsampling::approxtoexact_constrained_func(n=200, w=approximate_design$w, m=6, 
                                                           beta=beta, link='logit', X=X, 
                                                           Fdet_func=Fdet_func_GLM, 
                                                          iset_func=iset_func_trial, label=label) 

print(exact_design)

## -----------------------------------------------------------------------------
w00 = rep(1/200, 6)
unif_design = CDsampling::approxtoexact_constrained_func(n=200, w=w00, m=6, beta=NULL, 
  link=NULL, X=NULL, Fdet_func=Fdet_func_unif, iset_func=iset_func_trial)

print(unif_design)

## -----------------------------------------------------------------------------
J=5; p=12; m=8; #response levels; num of parameters; num of design points
beta = c(-4.047, -0.131, 4.214, -2.225, -0.376, 3.519, -0.302,
    -0.237,  2.420, 1.386,  -0.120,  1.284) #coefficients
Xi=rep(0,J*p*m) #design matrix
dim(Xi)=c(J,p,m)
Xi[,,1] = rbind(c( 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

Xi[,,2] = rbind(c( 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

Xi[,,3] = rbind(c( 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,4] = rbind(c( 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,5] = rbind(c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,6] = rbind(c( 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,7] = rbind(c( 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

 Xi[,,8] = rbind(c( 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 1),
                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

## -----------------------------------------------------------------------------
nsample=600 #sample size
constraint = c(392, 410)  #mild:severe

#boundaries r_i1 and r_i2 in constrained lift-one
lower.bound <- function(i, w0){
  n = 600
  constraint = c(392,410)
  if(i <= 4){
    a.lower <- (sum(w0[5:8])-(constraint[2]/n)*(1-w0[i]))/(sum(w0[5:8]))
  }
  else{
    a.lower <- (sum(w0[1:4])-(constraint[1]/n)*(1-w0[i]))/(sum(w0[1:4]))
  }
  a.lower
}

upper.bound <- function(i, w0){
  n = 600
  constraint = c(392,410)
  if(i <= 4){
    b.upper <- ((constraint[1]/n)*(1-w0[i]) - (sum(w0[1:4])-w0[i]))/(1-sum(w0[1:4]))
  }
  else{
    b.upper <- ((constraint[2]/n)*(1-w0[i]) - (sum(w0[5:8])-w0[i]))/(1-sum(w0[5:8]))
  }
  b.upper
}

#constraints of the example
g.con = matrix(0,nrow=length(constraint)+1+m, ncol=m)
g.con[2:3,] = matrix(data=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1), ncol = m, byrow=TRUE)
g.con[1,] = rep(1, m)
g.con[4:(length(constraint)+1+m), ] = diag(1, nrow=m)
g.dir = c("==", "<=","<=", rep(">=",m))
g.rhs = c(1, ifelse((constraint/nsample<1),constraint/nsample,1), rep(0, m))


## -----------------------------------------------------------------------------
label = c("Placebo-Mild", "Low-Mild", "Medium-Mild", "High-Mild", "Placebo-Severe", "Low-Severe", "Medium-Severe", "High-Severe")

## -----------------------------------------------------------------------------
set.seed(123)
approx_design = CDsampling::liftone_constrained_MLM(m=m, p=p, Xi=Xi, J=J, beta=beta, 
                                       lower.bound=lower.bound, upper.bound=upper.bound, 
                                        g.con=g.con, g.dir=g.dir, g.rhs=g.rhs, label=label, w00=NULL, 
                                        link='cumulative', Fi.func = Fi_func_MLM, 
                                        reltol=1e-5, maxit=500, delta = 1e-6, 
                                        epsilon=1e-8, random=TRUE, nram=1)

print(approx_design)

## -----------------------------------------------------------------------------
exact_design = CDsampling::approxtoexact_constrained_func(n=600, w=approx_design$w, m=8, 
 beta=beta, link='cumulative', X=Xi, Fdet_func=Fdet_func_MLM, 
 iset_func=iset_func_trauma, label=label)

print(exact_design)

## -----------------------------------------------------------------------------
iset_func_trauma <- function(allocation){
  iset = rep(1,8)
  if(sum(allocation[1:4])>=392){iset[1:4]=0}
  if(sum(allocation[5:8])>=410){iset[5:8]=0}
  return(iset)
 }

unif_design = CDsampling::approxtoexact_constrained_func(n=600, w=rep(1/600,8), m=8, beta=NULL, 
 link=NULL, X=NULL, label=label, Fdet_func=Fdet_func_unif, iset_func=iset_func_trauma)

unif_design

