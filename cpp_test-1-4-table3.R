#
# modified on 02/15/20 by Y. Zuo
#rm(list=ls())
library(mvtnorm)
library(robustbase)
library(mrfDepth)
library(matlib)
library(matrixcalc)
library(compiler)
enableJIT(1)

options( warn = -1 )

K=3 #; alpha=0.95
n=65; p=4# these parameters can be set as one wanted every time. 
RepN=100; # number of replication in the simultaion
#Nbet1=300 #adjustable parameters
#epsilon=0.05; # contamination percentage
#n1=floor(epsilon*n) #adjustable parameters

ND=600; Nbet=800;  #ND is the total number of random directions used for UFV, 
#Nbet is the total number of beta slected in the convex hull, adjustable parameters
RN=800;  #RN is the total# of selected beta from Z below
N1=900; a=(n+300)*(p-1); b=choose(n,p)
N=min(a, b, N1)  #N is the total numer of beta in y=x'beta determined by p points
#N could  also be inputed from outside
#N=???             #how many beta's you try to selected from choose(n, p) 
c1=min(RN, N)
RN=max(c1, (p+1)) #RN should be no less than (p+1) no greater than N

beta_ltsReg=beta_RD=beta_PRD1=beta_PRD2=beta_PRD3=matrix(0,RepN, ncol=p)
t_ltsReg=t_PRD1=t_PRD2=t_RD=t_PRD3=0;

####### three methods to compute the regression lines/hyperplanes ###########
#i=1
for (i in 1:RepN)
{
  # generating data  
  #beta_0=-2; beta_1=0.1;beta_2=1; beta_3=5; beta_4=10
  beta_0=-2; beta_1=0.1; beta_2=1; # beta_3=5
  #beta_0=50; beta_1=0.1;beta_2=-2; beta_3=15; beta_4=100
  
  #beta_zero=c(beta_0, beta_1, beta_2, beta_3, beta_4)
  beta_zero=c(beta_0, beta_1, beta_2)#, beta_3)
  
  #x1=matrix(rcauchy(n), nrow=n, ncol=1)
  x1=matrix(rnorm(n), nrow=n, ncol=1)
  #x1=matrix(runif(n), nrow=n, ncol=1)
  x2=matrix(rcauchy(n), nrow=n, ncol=1)
  #x2=matrix(rnorm(n), nrow=n, ncol=1)
  #x2=matrix(runif(n), nrow=n, ncol=1)
  x3=matrix(rcauchy(n), nrow=n, ncol=1)
  #x3=matrix(rnorm(n), nrow=n, ncol=1)
  #x3=matrix(runif(n), nrow=n, ncol=1)
  #x4=matrix(rcauchy(n), nrow=n, ncol=1)
  #x4=matrix(rnorm(n), nrow=n, ncol=1)
  #e=matrix(rnorm(n), nrow=n, ncol=1)
  e=matrix(rcauchy(n), nrow=n, ncol=1)
  #y=matrix(beta_0, nrow=n, ncol=1)+beta_1*x1+beta_2*x2+beta_3*x3+ beta_4*x4 +e
  y=e+matrix(beta_0, nrow=n, ncol=1)+beta_1*x1+beta_2*x2 +beta_3*x3 #+e #+ beta_4*x4 +e
  #m1=matrix(cbind(x1, x2, x3, x4, y), nrow=n, ncol=p)
  m1=matrix(cbind(x1, x2, x3, y), nrow=n, ncol=p)
  #m1=matrix(cbind(x1, x2, y), nrow=n, ncol=p)
  # m1=rmvnorm(n, mean=(rep(0, p)), sigma=diag(1:p))#,sigma=matrix(c(9, 0.9, 0.9, 1), ncol=2, byrow=T))
  # if(n1>=1)
  # {
  #   m2=rmvnorm(n1, mean=rep(10,p),sigma=diag(rep(0.1, p)))
  #   m1[sample(1:n,n1),]<-m2   #m1 now becomes 5% contanimated normal data
  # }

#-----------------------------------------------------------------------
  #########################
  # first ltsReg line
  #######################  
  if(n>2*p) #ltsReg requires n > 2*p
  {  
    t1=Sys.time()  
    fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
    beta_ltsReg[i,]=as.numeric(fit1$coefficients) 
    t2=Sys.time()-t1
    t_ltsReg=t2+t_ltsReg
  }
  
  ###############################
  # second RD induced line $T^*_{RD}$
  ###############################
  #install.packages("mrfDepth")
  Z=get_N_beta(m1, N) #the larger N being the better
  t1=Sys.time()
  # Z=get_all_beta(m1);
  k=dim(Z)[1] # total length of Z
  depth_beta=rdepth(m1, Z)$depthZ # rdepth of all beta from Z
  id=order(depth_beta)  # id of ordered beta with ascending depth
  max_RD=depth_beta[id[k]]
  multiplicity=sum(depth_beta==max_RD)
  all_deepest_beta=Z[depth_beta==max_RD,]
  # if(length(all_deepest_beta[,1])>1)
  #       beta_RD[i,]=colMeans(all_deepest_beta) #maximum_rdepth_beta
  # else 
  beta_RD[i,]=Z[id[k],] #beta with maximum depth
  # fit2=rdepthmedian(m1)$deepest
  # beta_RD[i,]=as.numeric(fit2) 
  t2=Sys.time()-t1
  t_RD=t2+t_RD
  ################################
  ## third line T^*_PRD line (Z19)
  ##########################
  UF_min=10^10
  B1=matrix(0, (RN+0), ncol=p); UFbeta1=rep(0, (RN+0)) #store UF and beta's
  t1=Sys.time()
  #Z=get_N_beta(m1, N) #all beta determined by p points from data m1 via y=x'beta
  
  id=rep(0, RN); UF_min_new=rep(0, (RN-p))
  for (j in 1:RN) #RN should not be greater than N
  {  id[j]=sample(1:N, 1);beta=Z[id[j],]; B1[j,]=beta  #select a beta from Z and build B 
  UFbeta1[j]=compute_UF(m1,beta,ND,UF_min)  # to build a vector UFbeta (RN by 1)
  if (j>=(p+1))
  {
    jj=j-p
    if (jj==1)
      UF_min_new[jj]=min(UFbeta1[1:j]) 
    else                            
      UF_min_new[jj]=min(UF_min, UFbeta1[j])
    if (UF_min_new[jj]< UF_min) 
      UF_min=UF_min_new[jj]
  }
  }
  # above UFbeta and B might be too short in some cases, so in the following
  # we consider more candidate beta just sample from mvnormal distribution
  
  # UF_temp=rep(0, Nbet1); temp=matrix(0, nrow=Nbet1,ncol=p)
  # for (ii in 1:Nbet1)  
  # {
  #  temp[ii,]=rmvnorm(1, rep(0,p)) 
  #  
  #  UF_temp[ii]=compute_UF(m1,temp[ii,],ND, UF_min)
  #  jj=ii+RN; UFbeta[jj]=UF_temp[ii];B[jj,]=temp[ii,]
  #  
  #  if (UF_temp[ii]<UF_min) {UF_min=UF_temp[ii]}
  # }
  
  beta_PRD1[i,]=compute_deepest_PRD1(m1,B1,UFbeta1,ND,Nbet,UF_min)$beta
  t2=Sys.time()-t1
  t_PRD1=t_PRD1+t2
  #print(i)
  ############################################
  #fourth line  T^*_PRD ZZ20 from codes-4-table-3.7-1.R
  #######################################
  N2=1+multiplicity # numbers of (ltsReg + T*_RD)
  dpbeta=rbind(beta_ltsReg[i,], all_deepest_beta)
  B2=rbind(Z[1:RN,], dpbeta) # all beta
  UFbeta2=rep(0, (RN+N2)) #store UF of beta's
  t1=Sys.time()
  
  temp=get_Beta_UFbeta_cpp(m1, Z[1:RN,], dpbeta, ND, N2, RN)
  UFbeta2=temp$"UFbeta2"
  UF_min=temp$"UF_min"
  
  t11=Sys.time()

  beta_PRD2[i,]=c(compute_deepest_PRD2_cpp(m1,B2,UFbeta2,ND,Nbet,UF_min, RN))
  
  t12=Sys.time()
  
  temp1=get_dprd_wprd(B2,UFbeta2, K)
#  beta_PRD2[i,]=temp1$"dprd"
#  beta_PRD2[i,]=temp1$"wprd1"
  beta_PRD3[i,]=temp1$"wprd2"
  
  t13=Sys.time()
  
  t2=t12-t1
  t3=t13-t12+ t11-t1
  t_PRD2=t_PRD2+t2
  t_PRD3=t_PRD3+t3
  
  print(i)
} #bigfor loop end
  
#EMSE_ltsREg= sum(beta_ltsReg*beta_ltsReg)/RepN
#EMSE_RD= sum(beta_RD*beta_RD)/RepN
#EMSE_PRD1= sum(beta_PRD1*beta_PRD1)/RepN
#EMSE_PRD2= sum(beta_PRD2*beta_PRD2)/RepN
#EMSE_PRD3= sum(beta_PRD3*beta_PRD3)/RepN
print(beta_zero)
beta_matrix=matrix(beta_zero, byrow=T, nrow=RepN, ncol=p)
beta=beta_matrix
EMSE_ltsREg= sum((beta_ltsReg-beta)*(beta_ltsReg-beta))/RepN
EMSE_RD= sum((beta_RD-beta)*(beta_RD-beta))/RepN
EMSE_PRD1= sum((beta_PRD1-beta)*(beta_PRD1-beta))/RepN
EMSE_PRD2= sum((beta_PRD2-beta)*(beta_PRD2-beta))/RepN
EMSE_PRD3= sum((beta_PRD3-beta)*(beta_PRD3-beta))/RepN

print(list("n-p-RepN-N-ND-Nbet-RN-epsilon-K", 
           c(n,p,RepN,N,ND,Nbet,RN,epsilon, K)))
print(list("t1_lstReg, t1_RD, t1_PRD1, t1_PRD2, t_PRD3", 
           c(t_ltsReg, t_RD, t_PRD1,  t_PRD2, t_PRD3)/RepN))
print(list("EMSE_ltsREg, EMSE_RD, EMSE_PRD1,  EMSE_PRD2, EMSE_PRD3", 
           c(EMSE_ltsREg, EMSE_RD, EMSE_PRD1, EMSE_PRD2, EMSE_PRD3)))

##########################################################################
# function part need to be run first then run above part
##########################################################################  


library(mvtnorm)
library(robustbase)
library(mrfDepth)
library(matlib)
library(matrixcalc)
library(compiler)
enableJIT(1)
### added function ------------------
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
  y_temp=array(0, c(p,p,N))
  for (ii in 1:N)
  {  
    X_temp[,,ii] = X[sample(1:n,p),]
    y_temp[,,ii]=cbind(1, X_temp[,1:(p-1), ii])
    # while (is.singular.matrix(y_temp[,,ii]))
    # { 
    #   X_temp[,,ii] = X[sample(1:n,p),]
    #   y_temp[,,ii]=cbind(1, X_temp[,1:(p-1), ii])
    #   if (!is.singular.matrix(y_temp[,,ii]))
    #     break
    # }
    beta[ii,] = c( solve( y_temp[,,ii] ) %*% X_temp[,p, ii]  )
  }
  return(beta)
}
get_N_beta=cmpfun(get_N_beta1)
################################################################################
#first function
UFVhd1_R<-function(v,bet,X)
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
UFV=cmpfun(UFVhd1_R)

#---------------------------------------------------------------
#second function

compute_UF1_R<-function(X, beta, N, UF_min)
{ # given a deta set X (n by p) and a vector beta 1 by p, 
  # N is the random directions will used to compute the UFitness of beta
  n=dim(X)[1]; p=dim(X)[2];
  ufold=0  
  
  # compute the t_i or T matrix it is n by p matrix. t_i=w_i/r_i
  #--------------------------------------------------------------
  TM=matrix(0,n, ncol=p); DI=matrix(0, nrow=n,ncol=n)
  W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
  
  W=cbind(1, X[,1:(p-1)])  # n by p matrix
  Res=X[,p]-W%*%matrix(beta, nrow=p) # n by 1 vector of residuals
  Res[Res==0]=10^{-20} # deal with zero residuals
  Res=rep(1,n)/Res  # reciprocal of residuals
  
  DI=diag(as.vector(Res))
  TM=DI%*%W
#  print(TM[1:10,])
  #############################################################
  #specail directions that could lead to maximum UF
  # (i) ones perpendicular to T_i 
  # (ii) all p axis
  #  (iii) ND random vector, or better ones perpendicular to hyperplane determined 
  #        by p poins from T_i
  #     (iv) could also consider the directions of betai-betaj/\|betai-betaj\| for all
  #          i \neq j
  #############################################################
  #-----------------------------------------------------------------  
  #consider n vectors perpendicular to T_i=W_i/(yi-wi'beta)
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
    #if (i==1){print(uu); print(TM[i,])}
    ufnew=UFV(v,beta,X)
   
    #print(c('ufold--ufnew', ufold, ufnew)) 
    if (ufnew>=UF_min){ufold=10^10; return(ufold)}
    else {ufold=max(ufold, ufnew);}# print(c(ufold, ufnew))}
    #if (ufnew> ufold) {ufold=ufnew}
  } 
  #print(ufold)
  # ------------------------------------------------------  
  # consider p axis directions (special directions)
  D=diag(p)
  for (i in 1:p)
  {
    ufnew=UFV(D[,i],beta,X)
    if (ufnew>=UF_min){ufold=10^10;return(ufold)}
    else {ufold=max(ufold, ufnew)}
    #if (ufnew> ufold) {ufold=ufnew}
  }  
  #print(ufold)  
  #--------------------------------------------------------------
  # conside N random directions
  for (i in 1:N )
  { 
    # D=matrix(rmvnorm(1, rep(0,p)), ncol=2)
    # v=D/sqrt(D[,1]^2+D[,2]^2)
    # use perpendicular directions
    X_temp=X[sample(1:n, p),]
    vv=solve(cbind(1,X_temp[,1:(p-1)]))%*%X_temp[,p] #this is (beta_0, beta_1)' 
    # since y=(1,x')(beta_0, beta_1')'
    
    beta_1=vv[2:p,]
    uu=c(-t(beta_1), 1)
    #uu=rbind(-t(beta_1), 1)#nornal vector to the hyperplane determined by p points
    # since (-t(beta_1),1)%*%(x,y)'=beta_0   
    v=uu/sqrt(sum(uu*uu)) #the unit vector perpendicular the the hyperplne 
    
    ufnew=UFV(v,beta,X)
    if (ufnew>=UF_min){ufold=10^10;return(ufold)}
    else {ufold=max(ufold, ufnew)}
    #if (ufnew> ufold) {ufold=ufnew}
  }
  
  return(ufold)
}

compute_UF=cmpfun(compute_UF1_R)
#the third function
#--------------------------------------------------------------------
compute_deepest_PRD1<-function(m1,B,UFbeta,ND,Nbet, UF_min)
{ #m1 is the given data, UF_min is for speeding up and skipping some computation
  #B is the given beta vector dim(B)[1] by p and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull 
  index=order(UFbeta) # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  J=index[1:(p+1)] # (p+1) least UF's indices
  
  # final search over the convex hull
  UFbeta_final=matrix(0,Nbet,ncol=1); # UF of Nbet beta's over the convex hull
  beta=matrix(0,Nbet,ncol=p);
  for (k in 1:Nbet)
  {
    b=runif(p+1); a=b/sum(b); 
    
    #beta[k,]=a[1]*B[J1,]+a[2]*B[J2,]+a[3]*B[J3,];
    beta[k,]=a%*%B[J,]
    
    UFbeta_final[k]=compute_UF(m1,beta[k,],ND, UF_min) # more work here on the directions
  }
  ORD=order(UFbeta_final)
  if (UFbeta_final[ORD[1]]< UFbeta[J[1]])
    return(list(beta=beta[ORD[1],], UF=UFbeta_final[ORD[1]]))
  else 
    return(list(beta=B[J[1],], UF=UFbeta[J[1]]))
}  
compute_deepest_PRD1=cmpfun(compute_deepest_PRD1)
#the fourth function
#-------------------------------------------------------------------------
compute_deepest_PRD2_R<-function(m1,B,UFbeta,ND,Nbet, UF_min, RN)
{ #m1 is the given data, UF_min is for speeding up and skipping some computation
  #B is the given beta vector dim(B)[1] by p and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull formed by beta in B
  p=dim(B)[2]; n=dim(B)[1]; N2=n-RN
  index=order(UFbeta) # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  #  min_UF=min(UFbeta); id_min=which(UFbeta==min_UF)[1]
  id_min=index[1]
  J=index[1:(1+p)] # (p+1) least UF's indices
  
  # final search over the convex hull
  # 
  UFbeta_final=matrix(0,Nbet+N2,ncol=1) # UF of Nbet beta's over the convex hull
  beta=matrix(0,Nbet+N2,ncol=p)
  for (k in 1:Nbet)
  {
    #if (p==2)
    b=runif(p+1) 
    #else     
    #b=runif(dim(B)[1])  
    a=b/sum(b); 
    
    #beta[k,]=a[1]*B[J1,]+a[2]*B[J2,]+a[3]*B[J3,];
    #if (p>2) 
    #    J=sample(1:dim(B)[1], (p+1))
    beta[k,]=a%*%B[J,]
    #else 
    #beta[k,]=a%*%B
    UFbeta_final[k]=compute_UF(m1,beta[k,],ND, UF_min)# more work here on the directions
    UF_min=min(UF_min, UFbeta_final[k])
  }
  for (i in 1:N2)
  {
    b=runif(N2) 
    a=b/sum(b); 
    
    k=i+Nbet; J=(RN+1):n
    beta[k,]=a%*%B[J,]
    UFbeta_final[k]=compute_UF(m1,beta[k,],ND, UF_min) 
    UF_min=min(UF_min, UFbeta_final[k])
  }
  #ORD=order(UFbeta_final)
  min_UF_final=min(UFbeta_final)
  #min_UF_final=UF_min
  id_fin=which(UFbeta_final==min_UF_final)[1]
  #if (UFbeta_final[ORD[1]]< UFbeta[index[1]])
  if (UFbeta_final[id_fin]< UFbeta[id_min])
    return(list(beta=beta[id_fin,], UF=UFbeta_final[id_fin]))
  else 
    return(list(beta=B[id_min,], UF=UFbeta[id_min]))
  #return(list(beta=B[id[1],], UF=min_UF))
}  
compute_deepest_PRD2=cmpfun(compute_deepest_PRD2_R)


