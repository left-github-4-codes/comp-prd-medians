# exampl-3-1 added on 9/25/19 for revision of compu-PRD for csda
# adopted from codes-4-example-3-1.R for compu-PRD on 9/25/19
# scan the data to m1 first manually 

# modified on 7/29/19 for the final version for stpa after acception on 7/29/19
# modified on 7/8/19 for the 2rd revision for Stpa by Y. Zuo
# codes for the example 2 in Section 4 of RDREF 03/02/19 by Y.Zuo
# need to compute the LS and deepest lines from RD of Rousseeuw and Hurbert (1999)
# and from PRD of Zuo (2018).
##############################
#first run bottom functions 
#############################

rm(list=ls())
#library(mvtnorm)
#N=100 #number of replication

#Get data matrix Lung Cancer & Smoking

# > m1=scan()
# 1: 170  455
# 3: 150  510
# 5: 165  380
# 7: 350  1115
# 9: 465  1145
# 11: 58  220 
# 13: 245  460
# 15: 90  250 
# 17: 115  310
# 19: 250  530
# 21: 190  1280
# 23: 
#   Read 22 items
# > m1=matrix(m1, ncol=2, byrow=T)
# > m1
# [,1] [,2]
# [1,]  170  455
# [2,]  150  510
# [3,]  165  380
# [4,]  350 1115
# [5,]  465 1145
# [6,]   58  220
# [7,]  245  460
# [8,]   90  250
# [9,]  115  310
# [10,]  250  530
# [11,]  190 1280
# switch the column of m1 so that x first and then y
# > m1=cbind(m1[,2], m1[,1])
# > m1
# [,1] [,2]
# [1,]  455  170
# [2,]  510  150
# [3,]  380  165
# [4,] 1115  350
# [5,] 1145  465
# [6,]  220   58
# [7,]  460  245
# [8,]  250   90
# [9,]  310  115
# [10,]  530  250
# [11,] 1280  190


# # generate bivariate normal points with mean (0,0)' and covariance (1 -0.8; -0.8 1)
# # contaminated by (replacing) 34% points with mean (10, 10)' and covariance
# #diag(0.1 0.1)
# 
# m1=rmvnorm(n=N, mean=c(8,0),sigma=matrix(c(9, 0.9, 0.9, 1), ncol=2, byrow=T))
# m2=rmvnorm(n=34, mean=c(1,11),sigma=matrix(c(0.1,0,0,0.1), ncol=2, byrow=T))
# contami_M=rbind(m1[1:66,],m2)
# m3=contami_M

par(mfrow=c(1,2))
plot(m1,xlim=c(100,1400),ylim=c(0,600), xlab="cigarettes consumed in 1930",
     ylab="death per million in 1950",main="all eleven countries")
#plot(m2)
#plot(contami_M,xlim=c(-3,12),ylim=c(-3,12),xlab="x",ylab="y",main="contaminated")
############################LS line####################################
# compute the LS line for a given data set and draw the line
fit<-lm(m1[,2]~m1[,1])
abline(fit, col="red", lty=1)

print(c(fit$coefficients[[1]], fit$coefficients[[2]]))

#abline(c(fit$coefficients[[1]], fit$coefficients[[2]]), col="green"
###########################RD line#########################################
# install.packages(mrfDepth)
library(mrfDepth)
#compute the point-wise line (slope and intercept) call it as z
#use rdepth (X,z) to get the rd of all the lines and then find the one with maximum 
#RD 
#use mrfDepth:::rdepth to replace library call

Z=get_all_beta(m1)

k=dim(Z)[1] # total length of Z
depth_beta=rdepth(m1, Z)$depthZ # rdepth of all beta
id=order(depth_beta)  # id of beta with the maximum depth
#maxRD=depth_beta[id[k]]
#all_deepest_beta=Z[depth_beta==maxRD,]
#T_RD=colMeans(all_deepest_beta) #maximum_rdepth_beta
T_RD=Z[id[k],] #beta with maximum depth
abline(T_RD, col="blue", lty=2)
print(T_RD)
#abline(-0.0434626, -0.1468235, col="red")
############################PRD line##########################################

kk=dim(Z)[1]
#print(kk)
ND=100; Nbet=100; 
RN=300; B=matrix(0, RN, ncol=2)
UFbeta=rep(0, RN)

for (j in 1:RN)
{ii=sample(1:kk, 1);beta=Z[ii,]; B[j,]=beta
UFbeta[j]=compute_UF(m1,beta,ND)
}

T_PRD=compute_deepest_PRD(m1,B,UFbeta,ND,Nbet)
abline(T_PRD, col="black",lty=3)

print(T_PRD)

#The following is added on 9/23/19 for the revision of CSDA by Y.Zuo
################## ltsReg line ##########################################
#install.packages(robustbase)
library(robustbase)
fit1<-ltsReg(m1[,1], m1[,2])

abline(c(fit1$coefficients[[1]], fit1$coefficients[[2]]), col="green",lty=4)
print(c(fit1$coefficients[[1]], fit1$coefficients[[2]]))

print(c(fit$coefficients[[1]], fit$coefficients[[2]]))
print(T_RD)
print(T_PRD)

legend(105, 600, legend=c("LS", expression("T*"["RD"]),expression("T*"["PRD"],"ltsReg")),
       col=c("red", "blue", "black", "green"), lty=1:3, cex=0.7,
       title="Line types", text.font=4, bg='lightblue')

#######################################################################
#m2=rmvnorm(n=34, mean=c(1,11),sigma=matrix(c(0.1,0,0,0.1), ncol=2, byrow=T))
#contami_M=rbind(m1[1:66,],m2)
#m3=contami_M
m3=m1[1:10,]
contami_M=m3

plot(contami_M,xlim=c(100,1400),ylim=c(0,600),  xlab="cigarettes consumed in 1930",
     ylab="death per million in 1950",main="data without USA")
fit1<-lm(m3[,2]~m3[,1])
abline(fit1, col="red",lty=1)
print(c(fit1$coefficients[[1]], fit1$coefficients[[2]]))
###########################PD line#######################################
Z=get_all_beta(m3)

k=dim(Z)[1] # total length of Z
depth_beta=rdepth(m3, Z)$depthZ # rdepth of all beta
id=order(depth_beta)  # id of beta with the maximum depth
maxRD=depth_beta[id[k]]
all_deepest_beta=Z[depth_beta==maxRD,]
if (length(id)>1)
{T_RD=colMeans(all_deepest_beta)} #maximum_rdepth_beta
if (length(id)==1)
  T_RD=Z[id[k],] #beta with maximum depth
abline(T_RD, col="blue",lty=2)
print(T_RD)
########################PRD line###########################################
kk=dim(Z)[1]
#print(kk)
ND=100; Nbet=200; 
RN=500; B=matrix(0, RN, ncol=2)
UFbeta=rep(0, RN)

for (j in 1:RN)
{ii=sample(1:kk, 1);beta=Z[ii,]; B[j,]=beta
UFbeta[j]=compute_UF(m3,beta,ND)
}
T_PRD=compute_deepest_PRD(m3,B,UFbeta,ND,Nbet)
abline(T_PRD, col="black",lty=3)
print(T_PRD)

#The following is added on 9/23/19 for the revision of CSDA by Y.Zuo
################## ltsReg line ##########################################
#install.packages(robustbase)
#library(robustbase)
fit2<-ltsReg(m1[,1], m1[,2])

abline(c(fit2$coefficients[[1]], fit2$coefficients[[2]]), col="green",lty=4)
print(c(fit2$coefficients[[1]], fit2$coefficients[[2]]))

print(c(fit1$coefficients[[1]], fit1$coefficients[[2]]))
print(T_RD)
print(T_PRD)

legend(105, 600, legend=c("LS", expression("T*"["RD"]),expression("T*"["PRD"],"ltsReg")),
       col=c("red", "blue", "black", "green"), lty=1:3, cex=0.7,
       title="Line types", text.font=4, bg='lightblue')

#
############################################################################
library(mvtnorm)
library(compiler)
enableJIT(1)
get_all_beta<-function(X)
{                   #X is n by p matrix
  n=dim(X)[1];m=X
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
###################################################################
UFV1<-function(v,bet,X)
{
  #v,and beta are 1x2 vector, and X is nx2
  W=matrix(cbind(1,X[,1]), ncol=2);
  #wi=(1,xi')',w is a n by 2 matrix with 1st column is 1
  #if(is.na(W[1,1])) print("something wrong")
  Nvect=X[,2]-W%*%(bet) # this is suppose to be a n by 1 vector
  Dvect=matrix(0,dim(X)[1],ncol=1)+ W%*%matrix(v,nrow=2,ncol=1)  # this is suppose to be a n by 1 vector
  Dvect[Dvect[]==0]=10^{-20}
  quotient=matrix(Nvect/Dvect,nrow=1)
  #if (is.na(quotient[1])) print("something wrong1")
  ufv=median(quotient);
  #if (is.na(ufv)) print("something wrong2")
  abs(ufv) 
}
UFV=cmpfun(UFV1)
##########################################################################

compute_UF<-function(X,beta, N)
{# given a deta set X (n by 2) and a vector beta 1 by 2, 
  # N is the random directions will used to compute the UFitness of beta
  library(mvtnorm)
  n=dim(X)[1]; p=dim(X)[2];
  ufold=0  
  
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
    ufnew=UFV(v1,beta,X)
    if (ufnew>ufold){ufold=ufnew}
  }  
  # ------------------------------------------------------  
  # consider p axis directions
  D=diag(p)
  for (i in 1:p)
  {
    ufnew=UFV(D[,i],beta,X)
    if (ufnew> ufold) {ufold=ufnew}
  }  
  #-------------------------------------------------------------
  #consider N additional directions, which are the normal vectors to hyperplanes
  #determined by p points from T
  
  #install.packages("pracma")
  #require(pracma)
  #N=length(U[,1])
  # for (i in 1:N){
  #   id=sample(n,p)
  #   M=TM[id,] #p by p matrix
  #   point_wise_diff=(M-matrix(rep(M[p,], p),byrow=T,nrow=p))[1:(p-1),] # p-1 by p point-wise-different
  #   v=c(-point_wise_diff[,2], point_wise_diff[,1])
  #   #if (p==2)v=null(t(point_wise_diff)) #point_wise_diff%*%v=0, 
  #   #else v=null(point_wise_diff)
  #   # v is suppose a p by 1 vector if rank of point_wise_diff is (p-1)
  #   
  #   u=v
  #   v=u/sqrt(sum(u*u))
  #   ufv=UFV(v,beta,X)
  #   if (ufv> ufold) {ufold=ufv}
  # }  
  
  for (i in 1:N )
  { 
    D=matrix(rmvnorm(1,c(0,0)), ncol=2)
    v=D/sqrt(D[,1]^2+D[,2]^2)
    #print(X)
    ufv=UFV(v,beta,X)
    if (ufv> ufold) {ufold=ufv}
  }
  
  return(ufold)
}
# compute_UF1<-function(beta, .X, .U)
# {# given a deta set X (n by p), first (p-1) columns are x and the last column is y
#   # and a vector beta p by 1,U is the matrix of random directions from get_all_beta
#   X=.X; U=.U
#   
#   #library(mvtnorm)
#   ufold=0  
#   #n=dim(X)[1]; p=dim(X)[2]; 
#   # TM=matrix(0,n, ncol=p)
#   # W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
#   # 
#   # for (i in 1:n)
#   # {    
#   # W[i,]=cbind(1, t(X[i,1:(p-1)])) #1 by p matrix
#   # Res[i]=X[i,p]-sum(W[i,]* beta)  # residual, 1 by 1 matrix
#   # if( Res[i]==0 ) {Res[i]=10^{-20}} #take care of zero 
#   # TM[i,]=W[i,]/(Res[i]*rep(1,p))  #T_i in the Compu_PRD article
#   # }
#   
#   TM=matrix(0,n, ncol=p); DI=matrix(0, nrow=n,ncol=n)
#   W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
#   
#   W=cbind(1, X[,1:(p-1)])
#   Res=X[,p]-W%*%matrix(beta, nrow=p)
#   Res[Res==0]=10^{-20}
#   Res=rep(1,n)/Res
#   
#   DI=diag(as.vector(Res))
#   TM=DI%*%W
#   TM
#   #specail directions that could lead to maximum UF
#   #-----------------------------------------------------------------  
#   #consider 2n vectors perpendicular to T_i=(yi-wi'beta)/r_i
#   for (i in 1:n)
#   {
#     id=which(TM[i,]!=0) # subscript of non-zero elements
#     if (id[1]==1) {uu1=c(-TM[i,2],TM[i,1]); uu=rep(0,p); uu[1:2]=uu1}
#     else {uu1=c(-TM[i,id[1]], TM[i,1]); uu=rep(0,p); uu[c(1,id[1])]=uu1}
#     # using the idea of two-component vector (x[1], x[2]) and (-x[2], x[1])
#     # are perpendicular with additional non-zero restriction
#     # uu is non-zero and perpemdicular to TM[i,]
#     
#     v=uu/sqrt(sum(uu*uu))
#     eps=10^{-8}
#     
#     v1=v+eps*rep(1,p) #to avoid zero denominator in w'u in UFVhd
#     v1=v1/sqrt(sum(v1*v1))
#     ufnew=UFVhd(v1,beta,X)
#     if (ufnew>ufold){ufold=ufnew}
#     
#     #  v2=v-eps*rep(1,p) #to avoid zero denominator in w'u in UFVhd
#     #  v2=v2/sqrt(sum(v2*v2))
#     #  ufnew=UFVhd(v2,beta,X)
#     #  if (ufnew>ufold){ufold=ufnew}
#   }  
#   # ------------------------------------------------------  
#   # consider p axis directions
#   D=diag(p)
#   for (i in 1:p)
#   {
#     ufnew=UFVhd(D[,i],beta,X)
#     if (ufnew> ufold) {ufold=ufnew}
#   }  
#   #-------------------------------------------------------------
#   #consider N additional directions, which are the normal vectors to hyperplanes
#   #determined by p points from T
#   
#   #install.packages("pracma")
#   #require(pracma)
#   #N=length(U[,1])
#   for (i in 1:N){
#     id=sample(n,p)
#     M=TM[id,] #p by p matrix
#     point_wise_diff=(M-matrix(rep(M[p,], p),byrow=T,nrow=p))[1:(p-1),] # p-1 by p point-wise-different
#     if (p==2)v=null(t(point_wise_diff)) #point_wise_diff%*%v=0, 
#     else v=null(point_wise_diff)
#     # v is suppose a p by 1 vector if rank of point_wise_diff is (p-1)
#     
#     u=v
#     #dimv=dim(v)[2]
#     #if (dimv==1) {u=v}
#     #else {u=v[,1]}
#     
#     v=u/sqrt(sum(u*u))
#     ufv=UFVhd(v,beta,X)
#     if (ufv> ufold) {ufold=ufv}
#   }  
#   #----------------------------------------------------------------    
#   #consider N additional directions, which are the normal vector to hyperplanes
#   #determined by p points from X  
#   
#   for (i in 1:N )
#   { 
#     v=U[i,]/sqrt(sum(U[i,]*U[i,]))
#     ufv=UFVhd(v,beta,X)
#     if (ufv> ufold) {ufold=ufv}
#   }
#   #----------------------------------------------  
#   return(ufold)
# }
# 
# compute_UF=cmpfun(compute_UF1) 
########################################################################
compute_deepest_PRD<-function(m1,B,UFbeta,ND,Nbet)
{ 
  #B is the given beta vector n by 2 and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull 
  index=order(UFbeta) # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  #print(index)
  J1=index[1] #the index of the largest UF of some beta
  J2=index[2]
  J3=index[3] #
  
  #print(c(J1,J2,J3))
  # final search over the convex hull
  UFbeta_final=matrix(0,Nbet,ncol=1); # N: the # for all possible betas in the triangle
  beta=matrix(0,Nbet,ncol=2);
  for (k in 1:Nbet)
  {
    b=runif(3); a=b/sum(b); 
    
    beta[k,]=a[1]*B[J1,]+a[2]*B[J2,]+a[3]*B[J3,];
    # produce ND random unit directions for calculation of UF of each beta
    D=matrix(rnorm(2*ND), ncol=2, byrow=T)
    #D=matrix(rmvnorm(D,c(0,0)),ncol=2, byrow=T)
    UD=D/sqrt(D[,1]^2+D[,2]^2)
    # calculate the UF along all UD 
    ufold=0
    for (i in 1:ND) 
    { ufv=UFV(UD[i,],beta[k,], m1)
    if (ufv> ufold) {ufold=ufv}
    }
    UFbeta_final[k]=ufold # more work here as the directions
  }
  ORD=order(UFbeta_final)
  #print(UFbeta_final)
  #return(c(beta[ORD[1],],UFbeta_final[ORD[1]]))
  beta[ORD[1],]  
}  
#############################################################################