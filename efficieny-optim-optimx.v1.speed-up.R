#This is even slower than the old one: efficiency-simulation.v4.R
#But need to compare the efficiency with LS, perhaps it is more stable and efficient

# try to use DEoptim first time in the simulation of efficiency of RPD by Y. Zuo on
# 5/2/19

# geaneral idea:
# generate N beta's from the thyperplane determined by p sample points, assume that
# data are Z=(x', y) in R^p, then obtain the normal vector of the hyperplane 
# determined by the p points, U=Null(matrix), if the dim(U)[2]=1, then if U[p]=0;
# beta= (0,0,..0) else beta=

# obtain the lower and upper range of these beta's using min and max
# lower=apply(beta, 1, min); upper=apply(beta, 1, max)
# Then use out= DEoptim(compute_UF, lower, upper, 
# DEoptim.control(NP = 20, itermax = 40, F = 1.2, CR = 0.7)  X=Z, U=U)
# PRD=out$optim[1], UF=out$optim[2]

rm(list=ls())
library(mvtnorm)
library(compiler)
enableJIT(1)
#--------------------------------------------------------------------
UFVhd1<-function(v,bet,X)
{ #v & bet are 1xp vector, X is nxp with the last column is y,i.e. X[i,]=(x' y)
  p=dim(X)[2]
  W=matrix(cbind(1,X[,1:(p-1)]), ncol=p); #w is a n by p matrix with 1st column is 1
  # W'=(1,x')
  Nvect=X[,p]-W%*%matrix(bet, nrow=p) #  n by 1 vector, residuals
  Dvect= W%*%matrix(v, nrow=p)  #  n by 1 vector
  Dvect[Dvect[]==0]=10^{-20}
  quotient=as.vector(Nvect)/as.vector(Dvect)
  #quotient=matrix(Nvect, nrow=1)
  #if (fg==1) 
  ufv=median(quotient)
#  else { #ufv=wm(quotient)
#    ufv=ptm(quotient) 
  abs(ufv) 
}
UFVhd=cmpfun(UFVhd1)
#----------------------------------------------------------------------------
get_N_beta1<-function(X, N) #produce N beta and unit vector from X  
{ #X is n by p matrix
  n=dim(X)[1]; p=dim(X)[2]; m=X
  B=NULL; U=NULL # for all betas  # store the normal vectors 
  
  if (p==2)
  {
   for (i in 1:(n-1))
    { M= matrix((m-matrix(rep(m[i,],n),ncol=2,byrow=T))[(i+1):n,], ncol=2)
    #left-lower matrix of vectors of Xj-Xi, j>i
    #dir=rbind(dir,M) #stack all directions from Xj-Xi
    # assign smallest number to the zeors of x-coordinate of Xj-Xi
    M[M[,1]==0,1]=10^{-20}
    b1=M[,2]/M[,1] #slope vactor (p-i) by 1
    b0=m[i,2]-m[i,1]*b1 #intercept vector (slope and intercept form)
    #from y-y1=(y2-y1)/(x2-x1)(x-x1) to get b1 and b0
    B=rbind(B,cbind(b0,b1)) #stack all slope and intercept from point Xj and Xi
    v=cbind(-M[,2], M[,1])
    u=v/sqrt(sum(v*v))
  #  print(u);
    U=rbind(U,u)  
    } #beta from all linesegments y=b0+b1x
  }
  else
  { 
  #install.packages("pracma")
  #library(pracma)
  for (i in 1:N){
    id=sample(n,p)
    M=m[id,] #p by p matrix
    
    point_wise_diff=(M-matrix(rep(M[p,], p),byrow=T,nrow=p))[1:(p-1),] # p-1 by p point-wise-different
    if (p==2)v=null(t(point_wise_diff)) #point_wise_diff%*%v=0, 
    else v=null(point_wise_diff) #point_wise_diff%*%v=0, 
    # v is suppose a p by 1 vector if rank of point_wise_diff is (p-1)
    #dimv=dim(v)[2]
    #print(dimv)
    #if (dimv==1) {u=v}
    #else {u=v[,1]}
    u=v; 
    U=rbind(U, u)
    
    #if(u[p]==0){ u[p]=10^{-20} } # take case if the pth compnent of v is zero
    b1=u[1:(p-1)]/u[p]
    b0=M[1,p]-sum(b1*M[1,1:(p-1)]) # figure out the slope and intercenpt
    beta=cbind(b0,b1)  # by the fact <v,(x',y)-(x0',y0)>=0
    B=rbind(B, beta)
    }# end for loop
   }# end of else
   list(B=B, U=U)
} #end function
get_N_beta=cmpfun(get_N_beta1)
#---------------------------------------------------------------
compute_UF1<-function(beta, .X, .U)
{# given a deta set X (n by p), first (p-1) columns are x and the last column is y
 # and a vector beta p by 1,U is the matrix of random directions from get_all_beta
  X=.X; U=.U
  
  #library(mvtnorm)
  ufold=0  
  #n=dim(X)[1]; p=dim(X)[2]; 
  # TM=matrix(0,n, ncol=p)
  # W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
  # 
  # for (i in 1:n)
  # {    
  # W[i,]=cbind(1, t(X[i,1:(p-1)])) #1 by p matrix
  # Res[i]=X[i,p]-sum(W[i,]* beta)  # residual, 1 by 1 matrix
  # if( Res[i]==0 ) {Res[i]=10^{-20}} #take care of zero 
  # TM[i,]=W[i,]/(Res[i]*rep(1,p))  #T_i in the Compu_PRD article
  # }
  
  TM=matrix(0,n, ncol=p); DI=matrix(0, nrow=n,ncol=n)
  W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
  
  W=cbind(1, X[,1:(p-1)])
  Res=X[,p]-W%*%matrix(beta, nrow=p)
  Res[Res==0]=10^{-20}
  Res=rep(1,n)/Res
  
  DI=diag(as.vector(Res))
  TM=DI%*%W
  TM
  #specail directions that could lead to maximum UF
#-----------------------------------------------------------------  
  #consider 2n vectors perpendicular to T_i=(yi-wi'beta)/r_i
  for (i in 1:n)
  {
    id=which(TM[i,]!=0) # subscript of non-zero elements
    if (id[1]==1) {uu1=c(-TM[i,2],TM[i,1]); uu=rep(0,p); uu[1:2]=uu1}
    else {uu1=c(-TM[i,id[1]], TM[i,1]); uu=rep(0,p); uu[c(1,id[1])]=uu1}
    # using the idea of two-component vector (x[1], x[2]) and (-x[2], x[1])
    # are perpendicular with additional non-zero restriction
    # uu is non-zero and perpemdicular to TM[i,]
  
    v=uu/sqrt(sum(uu*uu))
    eps=10^{-8}
  
    v1=v+eps*rep(1,p) #to avoid zero denominator in w'u in UFVhd
    v1=v1/sqrt(sum(v1*v1))
    ufnew=UFVhd(v1,beta,X)
    if (ufnew>ufold){ufold=ufnew}
    
  #  v2=v-eps*rep(1,p) #to avoid zero denominator in w'u in UFVhd
  #  v2=v2/sqrt(sum(v2*v2))
  #  ufnew=UFVhd(v2,beta,X)
  #  if (ufnew>ufold){ufold=ufnew}
  }  
# ------------------------------------------------------  
# consider p axis directions
  D=diag(p)
  for (i in 1:p)
  {
    ufnew=UFVhd(D[,i],beta,X)
    if (ufnew> ufold) {ufold=ufnew}
  }  
#-------------------------------------------------------------
  #consider N additional directions, which are the normal vectors to hyperplanes
  #determined by p points from T
  
  #install.packages("pracma")
  #require(pracma)
  #N=length(U[,1])
  for (i in 1:N){
    id=sample(n,p)
    M=TM[id,] #p by p matrix
    point_wise_diff=(M-matrix(rep(M[p,], p),byrow=T,nrow=p))[1:(p-1),] # p-1 by p point-wise-different
    if (p==2)v=null(t(point_wise_diff)) #point_wise_diff%*%v=0, 
    else v=null(point_wise_diff)
    # v is suppose a p by 1 vector if rank of point_wise_diff is (p-1)
    
    u=v
    #dimv=dim(v)[2]
    #if (dimv==1) {u=v}
    #else {u=v[,1]}
    
    v=u/sqrt(sum(u*u))
    ufv=UFVhd(v,beta,X)
    if (ufv> ufold) {ufold=ufv}
  }  
#----------------------------------------------------------------    
  #consider N additional directions, which are the normal vector to hyperplanes
  #determined by p points from X  
  
  for (i in 1:N )
  { 
    v=U[i,]/sqrt(sum(U[i,]*U[i,]))
    ufv=UFVhd(v,beta,X)
    if (ufv> ufold) {ufold=ufv}
  }
#----------------------------------------------  
  return(ufold)
}

compute_UF=cmpfun(compute_UF1) 
#----------------------------------------------------------------------------
get_all_beta1<-function(X) #just good for p=2
{                   #X is n by p matrix
  n=dim(X)[1]; m=X
  #dir=NULL; 
  B=NULL # for all betas
  
  for (i in 1:(n-1))
  { M= matrix((m-matrix(rep(m[i,],n),ncol=2,byrow=T))[(i+1):n,], ncol=2)
  #left-lower matrix of vectors of Xj-Xi, j>i
  #dir=rbind(dir,M) #stack all directions from Xj-Xi
  # assign smallest number to the zeors of x-coordinate of Xj-Xi
  M[M[,1]==0,1]=10^{-20}
  b1=M[,2]/M[,1] #slope vactor
  b0=m[i,2]-m[i,1]*b1 #intercept vector (slope and intercept form)
  #from y-y1=(y2-y1)/(x2-x1)(x-x1) to get b1 and b0
  B=rbind(B,cbind(b0,b1)) #stack all slope and intercept from point Xj and Xi
  } #beta from all linesegments y=b0+b1x
  B
}
get_all_beta=cmpfun(get_all_beta1)
#--------------------------------------------------------------------

#-----------------main part------------------------------------
library(mvtnorm)
require(pracma)
n=80; p=2; R=1000 #sample size and dimension  # replication number
a=(n+300)*p; b=nchoosek(n,p)
N=min(c(a, b)) # number of beta produced and directions used to compute UF
#Z1=rmvnorm(n, mean=rep(0,p), sigma=diag(p))
beta_T_PRD=matrix(0,R,ncol=p) #UF=rep(0, R)
beta_LS=matrix(0,R, ncol=p); beta_T_RD=beta_LS
#beta_T_PRD1=beta_T_PRD; UF1=UF

t1=Sys.time()
for (i in 1:R)
{Z1=rmvnorm(n, mean=rep(0,p), sigma=diag(p))
tprd1=Sys.time()
 out=get_N_beta(Z1,N)
 Beta=out$B; U=out$U
# lower1=apply(Beta, 2, min); upper1=apply(Beta,2,max)
# lower=rep(min(lower1),p); upper=rep(max(upper1),p)
# meanvect= apply(Beta,2, mean)
 meanvect=rep(0,p)
# medianvect=apply(Beta, 2, median)
# t1=Sys.time()
 require(optimx)
 meth0 <- c("Nelder-Mead")
 out=optimx(meanvect,compute_UF, method=meth0, 
            control=list(trace = F), .X=Z1, .U=U)
# out1=optim(medianvect, compute_UF, .X=Z1, .U=U)
# t2=Sys.time()
 beta_T_PRD[i,]=out$p1
tprd2=Sys.time()
print(tprd2-tprd1)
 # UF[i]=out$value
# beta_T_PRD1[i,]=out1$par; UF1[i]=out1$value
 # library(DEoptim)
 # t11=Sys.time()
 # out1=DEoptim(compute_UF,lower, upper, .X=Z1, .U=U, # DEoptim.control(trace = FALSE),)
 #      DEoptim.control(trace = T, itermax=50, F=1.2,CR=0.7) )
 # t22=Sys.time()
 # beta_T_PRD1[i,]=out1$optim[[1]]; UF1[i]=out1$optim[[2]] 
#----------------------------------------------------------------- 
 m1=Z1
 x=m1[,1]; y=m1[,2]; #this is just good for p=2
 #m1=Z_new
 #m1=cbind(x,y); XYF=c("xrn", "yrn")
 ###########################################################
tls1=Sys.time()
  fit<-lm(y~x)
 beta_LS[i,]=c(fit$coefficients[[1]], fit$coefficients[[2]])
tls2=Sys.time()
 print(tls2-tls1)
 ##---------------------------------------------------------------
trd1=Sys.time()
  Z=get_all_beta(m1)
 k=dim(Z)[1] # total length of Z
 library(mrfDepth)
 depth_beta=rdepth(m1,Z)$depthZ  # rdepth of all beta
 id=order(depth_beta)  # id of beta with the maximum depth
 T_RD=Z[id[k],] #beta with maximum depth
 beta_T_RD[i,]=T_RD
trd2=Sys.time() 
print(trd2-trd1)
 print(i)
  }  
t2=Sys.time()
EMSE_LS=sum(beta_LS*beta_LS)
EMSE_RD=sum(beta_T_RD*beta_T_RD)
EMSE_PRD1=sum(beta_T_PRD*beta_T_PRD)

t2-t1
c(EMSE_LS, EMSE_RD,  EMSE_PRD1)/R

print(EMSE_LS/EMSE_RD)
print(EMSE_LS/EMSE_PRD1)

print(c(n, p, N, R))
print("E-O-O.v1.speed-up.R")

summary(out, order=value)