# v4 drop the pwm on 4/27/19
# using med and wm simulataneously 4/23/19
#codes for example 3.1 of Compu-PRD adopted from RDREF efficiency-simulation.v.3.R
# by Y.Zuo on 4/21/19, using the algorithm in section 3, especially the (V) 
# of remarks 3.1 in Compu-PRD
#
#codes form relative efficiency in section 4.3 on 3/6/19 by Y.Zuo

# big loop for (i in 1:RN), RN is the total number of replications, start at 100 
# big matrices beta_LS_{RN X 2} for store b0 and b1 from LS, T_{PD} and T_{PRD}
# generate x_i from normal(0,1) and independtly y_i from N(0,1), i=1,...n
#
# compute LS bo and b1, fit$coefficients[[1]], fit$coefficents[[2]] beta_LS[i,]
# do the same for T_PD_beta[i,], and T_PRD_beta[i,]
#
# use GMCM:::colSds(x) to get the sds for LS_beta, ... T_PRD_beta
# slower but equivalent apply(x, 2, sd)
#
# finally, sd_LS/sd_T_RD and sd_LS/sd_PRD sign to effi_M[i,1], effi_M[i,2]
# end of the i loop
#
# final result colmeans(effi_M)
#############################################################
#pwm function
K=3
r0=3.5
wm<-function(z) #given a vector $z$ return a pwm(z)
{
  med<-median(z)
  ma<-median(abs(z-med))
  
  if(ma==0){r<-(z-med!=0)+0.0} #ris the PD
  else {r<-1/(1+abs(z-med)/ma)}
  w<-(r>=r0)+(r<r0)*(exp(-K*(1-r^2/r0^2)^2)-exp(-K))/(1-exp(-K))
  mu<-sum(w*z)/sum(w)
  return(mu)
}
###########################################################
#ptm function
cutoff=9 #5.5 #equivalent to PD=0.10
ptm<-function(z)
{
 mu=median(z)
 sigma=median(abs(z-mu))
 ptm=mean(z[abs(z-mu)/sigma<=cutoff])
 return(ptm)
}  
#--------------------------------------------------------------------
UFV<-function(v,bet,X,fg)
{ #fg is a flag variable, if fg=1, use med if fg=2 use pwm
  #v,and beta are 1x2 vector, and X is nx2
  W=matrix(cbind(1,X[,1]), ncol=2); #w is a n by 2 matrix with 1st column is 1
  # W'=(1,x')
  Nvect=X[,2]-W%*%(bet) #  n by 1 vector
  Dvect= W%*%matrix(v, nrow=2)  #  n by 1 vector
  Dvect[Dvect[]==0]=10^{-20}
  quotient=matrix(Nvect/Dvect,nrow=1)
  #ufv=median(quotient)
  if (fg==1) {ufv=median(quotient)}
  else { ufv=wm(quotient)}
  ##ufv=ptm(quotient)
  abs(ufv) 
}
##########################################################################
compute_UF<-function(X, beta, N,fg)
{# given a deta set X (n by 2) and a vector beta 1 by 2, 
  #N is the random directions will used to compute the UFitness of beta
  library(mvtnorm)
  ufold=0  
  n=dim(X)[1]; TM=matrix(0,n, ncol=2)
  
  W=matrix(cbind(1, X[,1]), ncol=2)
  Res=X[,2]-W%*%beta
  Res[Res==0]=10^{-20}
  #for (i in 1:n){TM[i,]=W[i,]/Res[i,]}   #t_i in the exact computation of UF 
  TM=W/cbind(Res,Res)
  for (i in 1:n)
  {
  uu=c(-TM[i,2], TM[i,1])
  v=uu/sqrt(uu[1]^2+uu[2]^2)
  eps=10^{-7}
  v1=v+eps*c(1,1)
  v1=v1/sqrt(v1[1]^2+v1[2]^2)
  ufnew=UFV(v1,beta,X, fg)
  if (ufnew>ufold){ufold=ufnew}
  v2=v-eps*c(1,1)
  v2=v2/sqrt(v2[1]^2+v2[2]^2)
  ufnew=UFV(v2,beta,X, fg)
  if (ufnew>ufold){ufold=ufnew}
  }  
  #consider the vectors perpendicular to t_i=(yi-wi'beta)/r_i 
  #consider N=100 perhaps is enough
  for (i in 1:N )
  { 
    D=matrix(rmvnorm(1,c(0,0),diag(2)), ncol=2)
    v=D/sqrt(D[,1]^2+D[,2]^2)
    ufv=UFV(v,beta,X,fg)
    if (ufv> ufold) {ufold=ufv}
  }
  return(ufold)
}
########################################################################
compute_deepest_PRD<-function(m1,B,UFbeta,ND,Nbet,fg)
{ #m1 is the original data matrix 
  #B is the given beta vector n by 2 and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull 
  index=order(UFbeta); # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  J1=index[1]; #the index of the largest UF of some beta
  J2=index[2];J3=index[3] #
  
  # final search over the convex hull
  UFbeta_final=matrix(0,Nbet,ncol=1) # N: the # for all possible betas in the triangle
  beta=matrix(0,Nbet,ncol=2); current_min_UF=UFbeta[J1]; #min UF so far
  for (k in 1:Nbet)
  {
    b=runif(3); a=b/sum(b); 
    
    beta[k,]=a[1]*B[J1,]+a[2]*B[J2,]+a[3]*B[J3,];
    # produce ND random unit directions for calculation of UF of each beta
    D=matrix(rnorm(2*ND), ncol=2, byrow=T);
    #D=matrix(rmvnorm(ND,c(0,0)),ncol=2, byrow=T)
    UD=D/sqrt(D[,1]^2+D[,2]^2);
    # calculate the UF along all UD 
    ufold=0;
    for (i in 1:ND) 
    { ufv=UFV(UD[i,],beta[k,], m1, fg)
    if (ufv> ufold) {
    ufold=ufv 
    cyrrent_min_UF=min(c(current_min_UF, ufold))
    }
    if(ufold>=current_min_UF) {break} #skip unnecessary directions
    }
    UFbeta_final[k]=ufold # more work here as the directions
  }
  ORD=order(UFbeta_final)
  #beta[ORD[1],]
  
  # modified by y.Zuo on 3/21/19 to take care of multiple deepest points
  min_UF=UFbeta_final[ORD[1]]
  #multiplicity=sum(UFbeta==min_UF)
  #print(multiplicity)
  idc=c(UFbeta_final==min_UF);
  nn=sum(idc)
  if (nn==1) 
    T_PRD=beta[ORD[1],]
  else
    T_PRD=colSums(beta[idc,])/sum(idc)
  #beta[ORD[1],]
  T_PRD  
}  
#----------------------------------------------------------------------------
get_all_beta<-function(X)
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


ep=0.0; fg=1
#ep=0.1; fg=1
library(mvtnorm)
N=80
RepN=1000; beta_LS=beta_T_RD=beta_T_PRD1=matrix(0,RepN, ncol=2)
beta_T_PRD2=beta_LS
#ND=100 #RN is the total# of replication, ND is the # of unit directions
#n=100 # n is the sample size
contami_M=NULL;
if (ep>0){contami_M=matrix(4, ncol=2, nrow=ep*N)}

t1=Sys.time()
for (i in 1:RepN)
{
#x=rnorm(n) #
#x=rt(n,2);y=rnorm(n);
DM=rmvnorm(n=N, mean=c(0,0), sigma=diag(2))
kk=(1-ep)*N;
m1=rbind(DM[1:kk,], contami_M)
x=m1[,1]; y=m1[,2];
#m1=Z_new
#m1=cbind(x,y); XYF=c("xrn", "yrn")
###########################################################
fit<-lm(y~x)
beta_LS[i,]=c(fit$coefficients[[1]], fit$coefficients[[2]])
############################################################
Z=get_all_beta(m1)
k=dim(Z)[1] # total length of Z
library(mrfDepth)
depth_beta=rdepth(m1,Z)$depthZ  # rdepth of all beta
id=order(depth_beta)  # id of beta with the maximum depth
T_RD=Z[id[k],] #beta with maximum depth
beta_T_RD[i,]=T_RD
#############################################################
 kk=dim(Z)[1]
 ND=100; Nbet=300;  #ND is the total # of directions, Nbet is the beta used
                    #in convex hull
 RN=800; B=matrix(0, RN, ncol=2) #RN is the total# of selected beta from Z
 UFbeta1=rep(1e+8, RN);  # fg=1 means med will be used
 UFbeta2=rep(1e+8,RN)   # fg=2 pwm be used
 #BR=matrix(0, ncol=2, nrow=2)
 #BR[1,]=1*apply(Z,2,min); BR[2,]=1*apply(Z,2,max)
 for (j in 1:RN)
 { ii=sample(1:kk, 1);beta=Z[ii,]; B[j,]=beta
   #b0=runif(1, BR[1,1], BR[2,1]); b1=runif(1, BR[1,2], BR[2,2])
   #beta=c(b0, b1); B[j,]=beta
   UFbeta1[j]=compute_UF(m1,beta,ND,1)
   #UFbeta2[j]=compute_UF(m1,beta,ND,2)
 }
 T_PRD1=compute_deepest_PRD(m1,B,UFbeta1,ND,Nbet,1)
 #T_PRD2=compute_deepest_PRD(m1,B,UFbeta1,ND,Nbet,2)
 beta_T_PRD1[i,]=T_PRD1
 #beta_T_PRD2[i,]=T_PRD2
}  
t2=Sys.time()
#t2-t1
EMSE_LS=sum(beta_LS*beta_LS)
EMSE_RD=sum(beta_T_RD*beta_T_RD)
EMSE_PRD1=sum(beta_T_PRD1*beta_T_PRD1)
#EMSE_PRD2=sum(beta_T_PRD2*beta_T_PRD2)

t2-t1
c(EMSE_LS, EMSE_RD, EMSE_PRD1)/RepN
#c(EMSE_LS, EMSE_RD)/RepN
print(EMSE_LS/EMSE_RD)

print(EMSE_LS/EMSE_PRD1)
#print(EMSE_LS/EMSE_PRD2)

print(list("RepN-N-RN-r0-K-Nbet-ep", c(RepN,N,RN,r0,K,Nbet,ep)))
