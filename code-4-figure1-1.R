## 12/27/2019 by Y. Zuo
##loading the MASS package, the data set is simply constructed by 

#rm(list=ls())
library(mvtnorm)
library(robustbase)
library(mrfDepth)
library(matlib)
#library(matrixcalc)
library(compiler)
enableJIT(1)


## install.packages(MASS)
library(MASS)
Animals2 <- local({D <- rbind(Animals, mammals); unique(D[order(D$body,D$brain),])})

## some of the following codes from internet
data(Animals2)
## Sensible Plot needs doubly logarithmic scale
plot(Animals2, log = "xy")

## Regression example plot:
#plotbb <- function(bbdat) {
#  d.name <- deparse(substitute(bbdat))
#  plot(log(brain) ~ log(body), data = bbdat, main = d.name)
#  abline(       lm(log(brain) ~ log(body), data = bbdat))
#  abline(MASS::rlm(log(brain) ~ log(body), data = bbdat), col = 2)
#  legend("bottomright", leg = c("lm", "rlm"), col=1:2, lwd=1, inset = 1/20)
#}
#plotbb(bbdat = Animals2)

## The `same' plot for Rousseeuw's subset:
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
#plotbb(bbdat = brain)

brain <- Animals[c(1:24, 26:25, 27:28),] #Rousseeuw and Leroy (1987) subset
#plot(brain, log = "xy")

#--------------------------------------------------------------

#--------------functions needed for RD and PRD-------------------------------------------
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
#-----------------------------------------------------------------------------
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

#---------------------------------------------------------------------
get_all_beta<-function(X) #just good for p=2
{                   #X is n by p matrix
  n=dim(X)[1]; m=X
  B=NULL # for all betas
  
  for (i in 1:(n-1))
  { mm=m-matrix(rep(m[i,],n),ncol=2,byrow=T)
    M= matrix(mm[(i+1):n,], ncol=2)
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
#----------------------------------------------------------------------
#-------------------------------------------------------------
##############updated compute_UF on 10/3/19  ##############################
compute_UF1<-function(X,beta, N, UF_min)
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
  #print(dim(Res))
  DI=diag(as.vector(Res))
  #print(c(dim(Res), dim(DI), dim(W)))
  TM=DI%*%W
  TM
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
    vv=solve(cbind(1,X_temp[,1:(p-1)]))%*%X_temp[,p] #this is (beta_0, beta_1)' 
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
#-------------------------------------------------------------
#fourth update the following updated on 10/3/19
########################################################################
compute_deepest_PRD1<-function(m1,B,UFbeta,ND,Nbet, UF_min)
{ #m1 is the given data, UF_min is for speeding up and skipping some computation
  #B is the given beta vector dim(B)[1] by p and UFbeta is the UF of the beta
  #ND is the total# of unit vectors used; Nbet is the # of beta searched over
  #the convex hull
  n=dim(m1)[1]; p=dim(m1)[2];
 
  index=order(UFbeta) # provide the position of the member in ascending order
  #i.e, indexes of the smallest, and of the 2nd smallest, ..., the largest
  J=index[1:(p+1)] # (p+1) least UF's indices
  
  # final search over the convex hull
  UFbeta_final=matrix(0,Nbet,ncol=1); # UF of Nbet beta's over the convex hull
  beta=matrix(0,Nbet,ncol=p);
  for (k in 1:Nbet)
  {
    b=runif(p+1); a=b/sum(b); 
    #print(c(length(a),a, (p+1)))
    #beta[k,]=a[1]*B[J[1],]+a[2]*B[J[2],]+a[3]*B[J[3],];
    beta[k,]=a%*%B[J,]
    #print(beta[k,])
    #print(c(k, beta[k,]))
    UFbeta_final[k]=compute_UF(m1,beta[k,],ND, UF_min) # more work here on the directions
  }
  ORD=order(UFbeta_final)
  if (UFbeta_final[ORD[1]]< UFbeta[J[1]])
    return(list(beta=beta[ORD[1],], UF=UFbeta_final[ORD[1]]))
  else 
    return(list(beta=B[J[1],], UF=UFbeta[J[1]]))
}  
compute_deepest_PRD=cmpfun(compute_deepest_PRD1)

#------------------------------------------------------------------------------------------------
#
par(mfrow=c(1,2))
plot(log(brain[,1]), log(brain[,2]), xlab="log(body)", ylab="log(brain)", main="Scatter plot of Animals")

plot(log(brain[,1]), log(brain[,2]), xlab="log(body)", ylab="log(brain)", main="Four line fits of Animals")

m1=cbind(log(brain[,1]), log(brain[,2]))

#abline(1.10957, 0.49601, lty=1, col=1, lwd=1) #LS line
#abline(0.86914, 0.75092, lty=1, col=1, lwd=1) #RLS line

#---------get line from LS----------------------------------------
fit1<-lm(m1[,2]~m1[,1])
abline(fit1, col=1, lty=1, lwd=1)
print(c(fit1$coefficients[[1]], fit1$coefficients[[2]]))

#-------------get ltsReg line------------------------------------
#install.packages("robustbase")
#library(robustbase)
#if(n>2*p) #ltsReg requires n > 2*p
#{  
# t1=Sys.time()  
fit2<-ltsReg(m1[,1], m1[,2])
abline(fit2$coefficients[[1]], fit2$coefficients[[2]], col=2, lty=2, lwd=1)
print(c(fit2$coefficients[[1]], fit2$coefficients[[2]]))
# beta_ltsReg[i,]=as.numeric(fit1$coefficients)  
# t2=Sys.time()-t1
# t_ltsReg=t2+t_ltsReg
#}


##------ get line from RD-------------------------------

Z=get_all_beta(m1) #all beta
k=dim(Z)[1] # total row(length) of Z
library(mrfDepth)
depth_beta=rdepth(m1,Z)$depthZ  # rdepth of all beta
id=order(depth_beta)  # id of beta with the maximum depth
#T_RD=Z[id[k],] #beta with maximum depth
max_RD=depth_beta[id[k]]
#multiplicity=sum(depth_beta==max_RD)
#print(multiplicity)
T_RD=colMeans(Z[depth_beta==max_RD,])
print(max_RD)
print(rdepth(m1,t(T_RD))$depthZ)
T_RD1=Z[depth_beta==max_RD,][1,]
#T_RD2=Z[depth_beta==max_RD,][2,]
#T_RD3=Z[depth_beta==max_RD,][3,]
#abline(Z[id[k],1], Z[id[k],2], lty=3, col=3, lwd=1) # RD line
#abline(T_RD, lty=3, col=3, lwd=1) #average deepest RD line
abline(T_RD1, lty=3, col=3, lwd=1) #average deepest RD line
print(T_RD1)
################################
## T^*_PRD line
##########################
UF_min=10^10
ND=500; Nbet=600;  #ND is the total # of directions, Nbet is the beta used in convex hull
RN=dim(Z)[1] #; B=matrix(0, RN, ncol=2) #RN is the total# of selected beta from Z
UFbeta1=rep(1e+8, RN);  # fg=1 means med will be used
#n=dim(m1)[[1]]; p=dim(m1)[[2]]
#a=(n+300)*p; b=choose(n,p)
#N=min(c(a, b)) 
N=RN
#Z1=m1
#t1=Sys.time()
#out=get_N_beta(Z1,RN); Beta=out$B; U=out$U               #modification needed here
for (j in 1:RN) {UFbeta1[j]=compute_UF(m1, Z[j,],ND, UF_min)} # modification neeeded here
out2=compute_deepest_PRD(m1,Z,UFbeta1,ND,Nbet, UF_min)
#t2=Sys.time()
#t2-t1
T_PRD1=out2$beta; UF=out2$UF
print(T_PRD1)
#print(UF)
abline(T_PRD1[1], T_PRD1[2], lty=4, col=4, lwd=1)
legend(-4, 8.8, legend=c("LS", "ltsReg", expression("T*"["RD"]),expression("T*"["PRD"])),    
          col=1:4, lty=1:4, cex=0.5,
          title="Line types", text.font=4, bg='lightblue')

