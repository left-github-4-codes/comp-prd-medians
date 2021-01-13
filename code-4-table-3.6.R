# Use return to replace break in compute_UF in codes-4-table-3.3.R
#
# speed up version of .2 or .1 by adding UF_min of first (p+1) beta's
# and therefore skipping some computation of future beta's, 10-10-19, Y.Zuo

# modified from codes-4-table-3-2.R by added append and rbind to UFbeta and B
# before calling compute_deepest_PRD, 10/5/19, by Y.Zuo
#
# Codes for Table 3 by Y.Zuo on 10/1/19 for revision of Compu-prd for CSDA
# compute (i) the EMSE of T^*_PRD and T^*_RD and ltsReg and 
# (ii) the average time per calculation replication R=1000, 
# EMSE =1/R\sum_{i=1}^R \|T_i-\bs{\beta}_0\|, with \bs{\beta}_0=0.
# 
###########################################################################
# run entire codes first, then ran top (no rm and no function part any more)
###########################################################################
rm(list=ls())
library(mvtnorm)
library(robustbase)
library(mrfDepth)
library(matlib)
library(compiler)
enableJIT(1)

n=20; p=6# these are tuning parameters, every time one can set one wanted. 
RepN=100; Nbet1=300
epsilon=0.1; n1=epsilon*n

ND=1500; Nbet=1500;  #ND is the total number of random directions used for UFV
#Nbet is the total number of beta slected in the convex hull
RN=1500;  #RN is the total# of selected beta from Z below
N1=1500; a=(n+300)*(p-1); b=choose(n,p)
N=min(c(a, b, N1))  #N is the total numer of beta in y=x'beta determined by p points
#N could  also be inputed from outside
#N=???             #how many beta's you try to selected from choose(n, p) 
c1=min(c(RN, N))
RN=max( c(c1, (p+1)) ) #RN should be no less than (p+1) no greater than N

beta_ltsReg=beta_RD=beta_PRD=matrix(0,RepN, ncol=p)
#beta_T_PRD2=beta_LS
#par(mfrow=c(1,2))
#plot(m1, type="p", pch=16, cex=1, xlab="x-axis", ylab='y-axis', 
#xlim=c(-5,15), ylim=c(-5,15))
t_ltsReg=t_PRD=t_RD=0;
####### three methods to compute the regression lines/hyperplanes ###########
for (i in 1:RepN)
{
  # generating data  
  m1=rmvnorm(n, mean=(rep(0, p)), sigma=diag(1:p))#,sigma=matrix(c(9, 0.9, 0.9, 1), ncol=2, byrow=T))
  if(n1>=1)
  {
    m2=rmvnorm(n1, mean=rep(10,p),sigma=diag(rep(0.1, p)))
    # contami_M=rbind(m1[1:(n-n1),],m2)
    # m1=contami_M
    m1[sample(1:n,n1),]<-rep(m2)
  }
  
#########################
# first ltsReg line
#######################  
if(n>2*p)
{  
 t1=Sys.time()  
 fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
 m=as.numeric(fit1$coefficients)
    #m=as.matrix(coef(fit1), ncol=1)
    #rownames(m)<-NULL
 beta_ltsReg[i,]=as.numeric(fit1$coefficients) #=t(m)
 t2=Sys.time()-t1
 t_ltsReg=t2+t_ltsReg
}
 #abline(c(fit1$coefficients[[1]], fit1$coefficients[[2]]), col="green",lty=4)
#print(c(fit1$coefficients[[1]], fit1$coefficients[[2]]))

###############################
# second RD induced line $T^*_{RD}$
###############################
#install.packages("mrfDepth")
Z=get_N_beta(m1, N)
t1=Sys.time()
# Z=get_all_beta(m1);
k=dim(Z)[1] # total length of Z
depth_beta=rdepth(m1, Z)$depthZ # rdepth of all beta
id=order(depth_beta)  # id of beta with the maximum depth
# maxRD=depth_beta[id[k]]
# all_deepest_beta=Z[depth_beta==maxRD,]
# beta_RD[i,]=colMeans(all_deepest_beta) #maximum_rdepth_beta
beta_RD[i,]=Z[id[k],] #beta with maximum depth
 # abline(T_RD, col="blue", lty=2)
 # print(T_RD)
#fit2=rdepthmedian(m1)$deepest
## m=matrix(fit2, ncol=p)
#beta_RD[i,]=as.numeric(fit2) #=m
t2=Sys.time()-t1
t_RD=t2+t_RD
#abline(c(T_RD[[1]], T_RD[[2]]), col="red", lty=1)

################################
## T^*_PRD line
##########################
UF_min=10^10
B=matrix(0, (RN+Nbet1), ncol=p); UFbeta=rep(0, (RN+Nbet1)) #store UF and beta's
t1=Sys.time()
#Z=get_N_beta(m1, N) #all beta determined by p points from data m1 via y=x'beta

id=rep(0, RN); UF_min_new=rep(0, (RN-p))
for (j in 1:RN) #RN should not be greater than N
{  id[j]=sample(1:N, 1);beta=Z[id[j],]; B[j,]=beta # select a beta from Z to get its UF 
   UFbeta[j]=compute_UF(m1,beta,ND,UF_min)  # to build a vector UFbeta (RN by 1)
   if (j>=(p+1))
   {
    jj=j-p
    UF_min_new[jj]=min(UFbeta[1:j])
    if (UF_min_new[jj]< UF_min) 
      UF_min=UF_min_new[jj]
   }
}
# above UFbeta and B might be too short in some cases, so in the following
# we consider more candidate beta just sample from mvnormal distribution

 UF_temp=rep(0, Nbet1); temp=matrix(0, nrow=Nbet1,ncol=p)
 for (ii in 1:Nbet1)  
 {
  temp[ii,]=rmvnorm(1, rep(0,p)) # replace these by betai-betaj/\|betai-betaj\| for all
                                # i \neq j
  UF_temp[ii]=compute_UF(m1,temp[ii,],ND, UF_min)
  if (UF_temp[ii]>=UF_min) {UFbeta[ii+RN]=UF_temp[ii];B[ii+RN,]=temp[ii,]; next}
  else
    {UF_min=UF_temp[ii]; UFbeta[ii+RN]=UF_temp[ii]; B[ii+RN,]=temp[ii,]}
 }
  
#T_PRD=compute_deepest_PRD(m1,B,UFbeta,ND,Nbet)
beta_PRD[i,]=compute_deepest_PRD(m1,B,UFbeta,ND,Nbet,UF_min)
#abline(T_PRD, col="black",lty=3)
#print(T_PRD)
t2=Sys.time()-t1
t_PRD=t_PRD+t2
print(i)
} #bigfor loop end

EMSE_ltsREg= sum(beta_ltsReg*beta_ltsReg)/RepN
EMSE_RD= sum(beta_RD*beta_RD)/RepN
EMSE_PRD= sum(beta_PRD*beta_PRD)/RepN

print(list("n-p-RepN-N-ND-Nbet-RN-epsilon-Nbet1", 
           c(n,p,RepN,N,ND,Nbet,RN,epsilon,Nbet1)))
print(list("t1_lstReg, t1_RD, t1_PRD", 
           c(t_ltsReg, t_RD, t_PRD)/RepN))
print(list("EMSE_ltsREg, EMSE_RD, EMSE_PRD", 
           c(EMSE_ltsREg, EMSE_RD, EMSE_PRD)))
############################################################################
#first update the following function
#library(mvtnorm)
#library(compiler)
#enableJIT(1)

##############################################################################
# get_all_beta<-function(X)
# { #X is n by p matrix
#   n=dim(X)[1]; p=dim(X)[2]; m=X
#   #dir=NULL; 
#   B=NULL # for all betas
#   
#   for (i in 1:(n-1))
#   { M= matrix((m-matrix(rep(m[i,],n),ncol=p,byrow=T))[(i+1):n,], ncol=p)
#   #left-lower matrix of vectors of Xj-Xi, j>i
#   #dir=rbind(dir,M) #stack all directions from Xj-Xi
#   # assign smallest number to the zeors of x-coordinate of Xj-Xi
#   M[M[,1]==0,1]=10^{-20}
#   b1=M[,2]/M[,1] #slope vactor
#   b0=m[i,2]-m[i,1]*b1 #intercept vector (slope and intercept form)
#   #from y-y1=(y2-y1)/(x2-x1)(x-x1) to get b1 and b0
#   B=rbind(B,cbind(b0,b1)) #stack all slope and intercept from point Xj and Xi
#   } #beta from all linesegments y=b0+b1x
#   B
# }
##############################updated get_N_beta on 10/3/19  ##################
get_N_beta1 = function (X, N) #by Shao wei under my general idea
{  
  # input: 
  #    X: is given data matrix with n rows and p columns, 
  #       the last column represents the y coordinates. 
  #    N: is the total number of beta one wants to obtained, 
  #      largest possible is (n choose p).
  # output:
  #    get_N_beta: all the beta result (N times p matrix)   
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  beta = matrix(0, nrow=N, ncol=p); X_temp=array(0, c(p,p, N))
#  X_num_combn = t( combn(1:n, p) ) # which is chosen number of N*p matrix
  for (ii in 1:N)
  {  
    beta = matrix(0, nrow=N, ncol=p); X_temp=array(0, c(p,p, N))
    #  X_num_combn = t( combn(1:n, p) ) # which is chosen number of N*p matrix
    for (ii in 1:N)
    {  
      X_temp[,,ii] = X[sample(1:n,p),]
      beta[ii,] = c( inv( cbind(1, X_temp[,1:(p-1), ii]) ) %*% X_temp[,p, ii]  )
    }
    return(beta)
  }  
}
get_N_beta=cmpfun(get_N_beta1)
###############################################################
# second update the following
###################################################################
# UFV1<-function(v,bet,X)
# {
#   #v,and beta are 1x2 vector, and X is nx2
#   W=matrix(cbind(1,X[,1]), ncol=2);
#   #wi=(1,xi')',w is a n by 2 matrix with 1st column is 1
#   #if(is.na(W[1,1])) print("something wrong")
#   Nvect=X[,2]-W%*%(bet) # this is suppose to be a n by 1 vector
#   Dvect=matrix(0,dim(X)[1],ncol=1)+ W%*%matrix(v,nrow=2,ncol=1)  # this is suppose to be a n by 1 vector
#   Dvect[Dvect[]==0]=10^{-20}
#   quotient=matrix(Nvect/Dvect,nrow=1)
#   #if (is.na(quotient[1])) print("something wrong1")
#   ufv=median(quotient);
#   #if (is.na(ufv)) print("something wrong2")
#   abs(ufv) 
# }
# UFV=cmpfun(UFV1)
##################updated UFV ###################################
UFVhd1<-function(v,bet,X)
{ #v & bet are 1xp vector, X is nxp with the last column is y,i.e. X[i,]=(x_i' y_i)
  p=dim(X)[2]
  W=matrix(cbind(1,X[,1:(p-1)]), ncol=p); #n by p with 1st column is 1: W'=(1,x')
  Nvect=X[,p]-W%*%matrix(bet, nrow=p) #  n by 1 vector, residuals: y-w'beta
  Dvect= W%*%matrix(v, nrow=p)  #  n by 1 vector: w'v
  Dvect[Dvect[]==0]=10^{-20}    # treat zero denominators in advance
  quotient=as.vector(Nvect)/as.vector(Dvect)
  ufv=median(quotient)
  abs(ufv) 
}
UFV=cmpfun(UFVhd1)

#third update the following
##########################################################################
# compute_UF<-function(X,beta, N)
# {# given a deta set X (n by 2) and a vector beta 1 by 2, 
#   # N is the random directions will used to compute the UFitness of beta
#   library(mvtnorm)
#   n=dim(X)[1]; p=dim(X)[2];
#   ufold=0  
#   
#   # compute the t_i ot T matrix it is n by p matrix. t_i=w_i/r_i
#   #--------------------------------------------------------------
#   TM=matrix(0,n, ncol=p); DI=matrix(0, nrow=n,ncol=n)
#   W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
#   
#   W=cbind(1, X[,1:(p-1)])  # n by p matrix
#   Res=X[,p]-W%*%matrix(beta, nrow=p) # n by 1 vector of residuals
#   Res[Res==0]=10^{-20} # deal with zero residuals
#   Res=rep(1,n)/Res  # reciprocal of residuals
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
#     ufnew=UFV(v1,beta,X)
#     if (ufnew>ufold){ufold=ufnew}
#   }  
#   # ------------------------------------------------------  
#   # consider p axis directions
#   D=diag(p)
#   for (i in 1:p)
#   {
#     ufnew=UFV(D[,i],beta,X)
#     if (ufnew> ufold) {ufold=ufnew}
#   }  
#   
#   for (i in 1:N )
#   { 
#     D=matrix(rmvnorm(1,c(0,0)), ncol=2)
#     v=D/sqrt(D[,1]^2+D[,2]^2)
#     #print(X)
#     ufv=UFV(v,beta,X)
#     if (ufv> ufold) {ufold=ufv}
#   }
#   
#   return(ufold)
# }
##############updated compute_UF on 10/3/19  ##############################
compute_UF1<-function(X,beta, N, UF_min)
{ # given a deta set X (n by p) and a vector beta 1 by p, 
  # N is the random directions will used to compute the UFitness of beta
  n=dim(X)[1]; p=dim(X)[2];
  ufold=0  
  
  # compute the t_i ot T matrix it is n by p matrix. t_i=w_i/r_i
  #--------------------------------------------------------------
  TM=matrix(0,n, ncol=p); DI=matrix(0, nrow=n,ncol=n)
  W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
  
  W=cbind(1, X[,1:(p-1)])  # n by p matrix
  Res=X[,p]-W%*%matrix(beta, nrow=p) # n by 1 vector of residuals
  Res[Res==0]=10^{-20} # deal with zero residuals
  Res=rep(1,n)/Res  # reciprocal of residuals
  
  DI=diag(as.vector(Res))
  TM=DI%*%W
  TM
  #############################################################
  #specail directions that could lead to maximum UF
  # (i) ones perpendicular to T_i 
  # (ii) all p axis
  #  (iii) ND random vector, or better ones perpendicular to hyperplane determined 
  #        by p poins from T_i
  #############################################################
  #-----------------------------------------------------------------  
  #consider n vectors perpendicular to T_i=(yi-wi'beta)/r_i
  for (i in 1:n)
  {
    id=which(TM[i,]!=0) # subscript of non-zero elements
    if (length(id)==1)
      {if (id[1]==1) {uu1=c(TM[i,2],-TM[i,1]); uu=rep(0,p); uu[1:2]=uu1}
      else {uu1=c(-TM[i,id[1]], TM[i,1]); uu=rep(0,p); uu[c(1,id[1])]=uu1}
      }  
    if(length(id)>1){
    uu=rep(0, p)
    uu[c(id[1],id[2])]<-c(-TM[i,id[2]], TM[i, id[1]]) #uui perpendicular to TM[i,]
    }
    #if (id[1]==1) {uu1=c(TM[i,2:p],-TM[i,1]); uu=rep(0,p); uu[1:2]=uu1}
    #else {uu1=c(-TM[i,id[1]], TM[i,1]); uu=rep(0,p); uu[c(1,id[1])]=uu1}
    # using the idea of two-component vector (x[1], x[2]) and (-x[2], x[1])
    # are perpendicular with additional non-zero restriction
    # uu is non-zero and perpemdicular to TM[i,]
    v=uu/sqrt(sum(uu*uu))
    # eps=10^{-8}
    # v1=v+eps*rep(1,p) #to avoid zero denominator in w'u in UFV
    # v1=v1/sqrt(sum(v1*v1))
    ufnew=UFV(v,beta,X)
    if (ufnew>=UF_min){ufold=10^10;return(ufold)}
    else
    if (ufnew>ufold){ufold=ufnew}
  }  
  # ------------------------------------------------------  
  # consider p axis directions (special directions)
  D=diag(p)
  for (i in 1:p)
  {
    ufnew=UFV(D[,i],beta,X)
    if (ufnew>=UF_min){ufold=10^10;return(ufold)}
    else
    if (ufnew> ufold) {ufold=ufnew}
  }  
  #--------------------------------------------------------------
  # conside N random directions
  for (i in 1:N )
  { 
    # D=matrix(rmvnorm(1, rep(0,p)), ncol=2)
    # v=D/sqrt(D[,1]^2+D[,2]^2)
    # use perpendicular directions
    X_temp=X[sample(1:n, p),]
    vv=inv(cbind(1,X_temp[,1:(p-1)]))%*%X_temp[,p] #this is (beta_0, beta_1)' 
    # since y=(1,x')(beta_0, beta_1')'
    
    beta_1=vv[2:p,]
    uu=rbind(-t(beta_1), 1)#nornal vector to the hyperplane determined by p points
    # since (-t(beta_1),1)%*%(x,y)'=beta_0   
    v=uu/sqrt(sum(uu*uu)) #the unit vector perpendicular the the hyperplne 
    
    ufnew=UFV(v,beta,X)
    if (ufnew>=UF_min){ufold=10^10;return(ufold)}
    else
    if (ufnew> ufold) {ufold=ufnew}
  }
  
  return(ufold)
}

compute_UF=cmpfun(compute_UF1)

#fourth update the following updated on 10/3/19
########################################################################
compute_deepest_PRD1<-function(m1,B,UFbeta,ND,Nbet, UF_min)
{ 
  #B is the given beta vector dim(B)[1] by p and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull 
  index=order(UFbeta) # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  J=index[1:(p+1)] # (p+1) least UF's indices
  # J1=index[1] #the index of the largest UF of some beta
  # J2=index[2]
  # J3=index[3] #
  
  #print(c(J1,J2,J3))
  # final search over the convex hull
  UFbeta_final=matrix(0,Nbet,ncol=1); # for all possible betas in the triangle
  beta=matrix(0,Nbet,ncol=p);
  for (k in 1:Nbet)
  {
    b=runif(p+1); a=b/sum(b); 
    
    #beta[k,]=a[1]*B[J1,]+a[2]*B[J2,]+a[3]*B[J3,];
    beta[k,]=a%*%B[J,]
    # produce ND random unit directions for calculation of UF of each beta
    #D=matrix(rnorm(p*ND), ncol=p, byrow=T)  #ND by p matrix
  #D=matrix(rmvnorm(ND,rep(0,p)))
  #temp=diag(1/sqrt(diag(D%*%t(D)))) # diag matrix with diag-entry are
    # 1/(norm of ith row)
  #UD=temp%*%D  # normalized D with unit vector in each row
    # calculate the UF along all UD[i,]
  # ufold=0
  # for (i in 1:ND) 
  # { ufv=UFV(UD[i,],beta[k,], m1)
  # if (ufv> ufold) {ufold=ufv}
  # }
    UFbeta_final[k]=compute_UF(m1,beta[k,],ND, UF_min) # more work here as the directions
  }
  ORD=order(UFbeta_final)
  #print(UFbeta_final)
  #return(c(beta[ORD[1],],UFbeta_final[ORD[1]]))
  beta[ORD[1],]  
}  
compute_deepest_PRD=cmpfun(compute_deepest_PRD1)
##########################################################################