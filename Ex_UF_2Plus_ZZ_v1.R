#rm(list=ls())
library(MASS)
Animals2 <- local({D <- rbind(Animals, mammals); unique(D[order(D$body,D$brain),])})
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
m1=cbind(log(brain[,1]), log(brain[,2]))
Z=m1

beta1=matrix(c(2.5548981, 0.4959947), 2, 1) # LS
beta2=matrix(c(2.0013467, 0.7508721), 2, 1) #ltsReg
#beta3=matrix(c(2.3899543, 0.6735954), 2, 1) #T_RD
#change to one of the other three lines, the average has the best UF
# results from code-4-table1.1.R
# > print(multiplicity)
# [1] 3
# >print(Z[depth_beta==max_RD,])
# b0        b1
# [1,] 2.258175 0.7028644
# [2,] 2.445328 0.6677692
# [3,] 2.466361 0.6501526
# > print(max_RD)
# [1] 0.4285714
# > print(max_RD)*28
# [1] 12
# > colMeans(Z[depth_beta==max_RD,])
# b0        b1 
# 2.3899543 0.6735954 
beta3=matrix(c(2.258175, 0.7028644), 2, 1) #T_RD
beta4=matrix(c(2.4509785, 0.6491982), 2, 1) # T_PRD

Ex_UF_2Plus_ZZ_v1<-function(Z,beta, Npoint)
{ # Z is data set n by p matrix, its first (p-1) columns are x_i
  # and the pth column is y_i, n is the sample size, beta is 1 x p vector
  # in parameter space (or candidate regression parameter)
  
  n=dim(Z)[1]; p=dim(Z)[2]
  gmatrinit=rep(0, 7);
  
  # compute the t_i or T matrix it is n by p matrix. t_i=w_i/r_i
  #--------------------------------------------------------------
  TM=matrix(0, n, ncol=p); DI=matrix(0, nrow=n, ncol=n)
  W=matrix(0, nrow=n, ncol=p); Res=rep(0,n)
  X=Z
  W=cbind(1, X[,1:(p-1)])  # n by p matrix
  Res=X[,p]-W%*%matrix(beta, nrow=p) # n by 1 vector of residuals
  Res[Res==0]=10^{-20} # deal with zero residuals
  Res=rep(1,n)/Res  # reciprocal of residuals
  DI=diag(as.vector(Res))
  TM=DI%*%W
  #-------------------------------------------------------------------
  #--------------------------------------------------------
  if (p==2)
  {
    # compute the pairwise difference as the vectors/dirs connecting two points
    # and record the corresponding id of the pair
    M=t(TM) #p by n matrix now
    N=n*(n-1)/2 #right upper triangle total
    pairwisedir=matrix(0,p,N) 
    #pairwise_id=matrix(0,N,2) #preassign to speed up
    for (i in 1:(n-1))
    {
      M_diff=M-matrix(M[,i],p, n) #difference to ith point, could use repmat
      s=(i-1)*n-i*(i-1)/2+1; t=i*n-i*(i+1)/2#positions of beginning and end
      pairwisedir[,s:t] =M_diff[,(i+1):n]  #(n-i) columns filled
      #pairwise_id[s:t,]=cbind(matrix(i, (n-i),1), matrix((i+1):n, (n-i),1))
    } #end #for loop entire loop corresponds to EA_UF2 step (2)  
    #--------------------------------------------------------------------------------
    perpend_slop=-pairwisedir[1,]/pairwisedir[2,]  #slope of N perpendicular directions
    alpha=atan(perpend_slop)+(perpend_slop<0)*pi # all alpha is in-between 0 and pi
    alpha_sort=sort(alpha) # could use unique to delete the repeated angles
                           # i.e. in the case of parallet line segments

    # initial direction 
    u0=matrix(c(1, 0), 2, 1) #equivalent to set alpha0=0;
    a0=0 #alpha[0]=0
    
    mid_index=matrix(0, N,2); u=matrix(0, 2, N) #temp=matrix(0, N, 2)
    temp0=get_mid_index(u0, TM)
    current_index=temp0$mid_ind

    u[,1]=matrix(c(cos(alpha_sort[1]), sin(alpha_sort[1])), 2, 1) # alpha1
    a1=alpha_sort[1] #alpha[1]
    #gmatrnew=update_gmatr(gmatrinit, u0, u[,1], TM) #compute the UF along the u0 first.
    gmatrnew=update_gmatr(gmatrinit, a0, a1, TM, Npoint)
    gmatrold=gmatrnew; u0=u[,1]
    a0=a1
#--------------------------------------------------    
     for(k in 1:n)
     { #TM is p by n matrix now   
       uu=c(-TM[k,2], TM[k,1])
       v=uu/sqrt(sum(uu^2))
       epsil=1e-7
       v1=v+epsil
       v1=v1/sqrt(sum(v1^2))
       if (v[1]>v1[1]) {v0=v1} else {v0=v} # v0 is the one with smaller polar angle 
       gmatrnew=update_gmatr(gmatrold, v0, v1, TM, Npoint)
       gmatrold=gmatrnew   
    #   print("v1")
    #   print(max(gmatrold))
    #    # v2=v-epsil
    #    # v2=v2/sqrt(sum(v2^2))
    #    # if (v[1]>v2[1]) {v0=v2} else {v0=v} # v0 is the one with smaller polar angle 
    #    # print(v0)
    #    # gmatrnew=update_gmatr(gmatrold, v0, v2, TM, Npoint)
    #    # gmatrold=gmatrnew
    #    # print("v2")
    #    # print(max(gmatrold))
     } 
#------------------------------------------------------  
    for (i in 1:(N-1)) # this part corresponds to the while part of Ex_UF_2Plus_Z19
    { 
       # u[,i]=matrix(c(cos(alpha_sort[i+1]), sin(alpha_sort[i+1])), 2, 1)
        a1=alpha_sort[i+1]
       # mid_index[i,]=get_mid_index(u0, TM)$mid_ind
       # 
       # if (sum(current_index %in% mid_index[i,])==2) 
       #  { current_index=mid_index[i,];  u0=u[,i]
       #    a0=a1; #print(i)
       #    next #skip this direction  
       #  }                                                                   
       #if (sum(current_index %in% mid_index[i,])<2)  
       #{ 
         #gmatrnew=update_gmatr(gmatrold, u0, u[,i], TM)
         gmatrnew=update_gmatr(gmatrold, a0, a1, TM, Npoint)
         gmatrold=gmatrnew;  
       #   u0=u[,i];
         a0=a1
       #    print(i) 
       #  current_index=mid_index[i,]
       #}   
    }# end of for loop  
  } #end of if (p=2)

  UF=max(gmatrold) 
  UF
 } #end of the main function
#------------------------------------------------------------------------------
get_mid_index=function(u, M)
{ # M is n by p data matrix, u is a p by 1 vector, the function returns the middle
  # two indeice of the projected values and the ordered projected values
  n=dim(M)[1]; n1=floor((n+1)/2); n2=floor((n+2)/2)
  
  proj_val=M%*%u
  temp=sort(proj_val, index.return=T)
  permu=temp$ix; sorted_val=temp$x
  mid_index=c(permu[n1], permu[n2])
  
  list(mid_ind=mid_index, sorted=sorted_val)
}  
#----------------------------------------------------------------------
update_gmatr<-function(gmatrold, a0, a1, TM, Npoint)
{# TM is n by p matrix, a0 and a1 are two angles p by 1 vector
  
  n=dim(TM)[1]; p=dim(TM)[2]; kk=floor((n+1)/2);
  u=matrix(c(cos(a0), sin(a0)), 2, 1)
  proj_val=TM%*%u; #%n by 1 vector
  N_u=sum(proj_val[]<0); #%total number of negative projected values
  #print(N_u)
  proj_val[proj_val==0]=10^{-20} #to avoid 0 in the denominator
  temp=sort(1/proj_val, index.return=T); kv=temp$x; uv=sort(proj_val)#sorted k^v
  permu=temp$ix; n1=kk; n2=floor((n+2)/2);
  
  gmatrnew=gmatrold;

  #construct B matrix, linear constraints Bu\leq 0
  # if (u[1]>0 && u1[1]>0 )
  # {B=matrix(c(u[2], -u[1], -u1[2], u1[1]), 2, 2,byrow=T)}
  # if (u[1]>0 && u1[1]<=0)
  # {B=matrix(c(u[2], -u[1], u1[2], -u1[1]), 2,2, byrow=T)}
  # if (u[1]<0 && u1[1]<0)
  # {B=matrix(c(-u[2], u[1], u1[2], -u1[1]), 2, 2,byrow=T)}    
  #   
  if (N_u==0|| N_u==n)
  { print("T0") 
    if (N_u==0)
  { print("T1")
    b=(TM[permu[n1],]+TM[permu[n2],])/2 #1 by p vector
    A=TM[permu[n1],]%*%t(TM[permu[n2],]) #p by p matrix
    # B=matrix(0, (n-1), p)
    # for (j in 1:(n-1))
    # {B[j,]=-(TM[permu[j],]-TM[permu[j+1],]) }   
    temp=maxi_over_pk(a0, a1, TM, Npoint) # a function returning the maximum valuse of |g(v)| 
    #over the current P_k, corresponding to the current permu
    gmatrnew[1]=max(abs(temp), gmatrold[1])
  }
    
    if ( N_u==n)
    { print("T2")
      b=(TM[permu[n1],]+TM[permu[n2],])/2 #1 by p vector
      A=TM[permu[n1],]%*%t(TM[permu[n2],]) #p by p matrix
      # B=matrix(0, (n-1), p)
      # for (j in 1:(n-1))
      # {B[j,]=(TM[permu[j],]-TM[permu[j+1],]) }   
      temp=maxi_over_pk(a0,a1,-TM, Npoint) # a function returning the maximum valuse of |g(v)| 
      #over the current P_k, corresponding to the current permu
      gmatrnew[2]=max(temp, gmatrold[2])
    }          
  }   
  
  if ((0<N_u)&&(N_u<n))            
  { u1=matrix(c(cos(a1),sin(a1)), 2,1)
    proj_val_1=TM%*%u1
    uv_1=sort(proj_val_1)
    N_u_1=sum(proj_val_1[]<0)
    if (n%%2==1 && kv[kk]<0)
    {gmatrnew[3]=max(gmatrold[3], -1/max(uv[N_u-kk+1], uv_1[N_u_1-kk+1]));   
     print("T[3]")
    }  #this and the one below cases are linear constraints, boundary
    #calculations are enough
  
    if (n%%2==1 && kv[kk]>0)
    {gmatrnew[4]=max(gmatrold[4], 1/min(uv[N_u+kk], uv_1[N_u_1+kk]));   
     print("T[4]") 
    }
  
    if (n%%2==0 && kv[kk]<0 && kv[kk+1]>0)
    { #print("T[5]")
      b=(TM[permu[n1],]+TM[permu[n2],])/2 #1 by p vector
      A=TM[permu[n1],]%*%t(TM[permu[n2],]) #p by p matrix
      # B=matrix(0, (n-1), p)
      # for (j in (1:N_u))
      # {B[j,]=(TM[permu[j],]-TM[permu[j+1],]) } 
      # 
      # for (j in (N_u+1):(n-1))
      # {B[j,]=-(TM[permu[j],]-TM[permu[j+1],]) }
      #print("T33")
      temp=maxi_over_pk(a0, a1, TM, Npoint) # a function returning the maximum valuse of |g(v)| 
      #over the current P_k, corresponding to the current permu
      gmatrnew[5]=max(abs(temp), gmatrold[5])
    } 
    if (n%%2==0 && kv[kk]>0)
    { #print("T[6]")
      b=(TM[permu[n1],]+TM[permu[n2],])/2 #1 by p vector
      A=TM[permu[n1],]%*%t(TM[permu[n2],]) #p by p matrix
      # B=matrix(0, (n-1), p)
      # for (j in (1:N_u))
      # {B[j,]=(TM[permu[j],]-TM[permu[j+1],]) } 
      # 
      # for (j in (N_u+1):(n-1))
      # {B[j,]=-(TM[permu[j],]+TM[permu[j+1],]) }
      temp=maxi_over_pk(a0, a1, TM, Npoint) # a function returning the maximum valuse of |g(v)| 
      #over the current P_k, corresponding to the current permu
      gmatrnew[6]=max(temp, gmatrold[6])
    }    
    if (n%%2==0 && kv[kk]<0 && kv[kk+1]<0)
    { #print("T[7]")
      b=t(TM[permu[n1],]+TM[permu[n2],])/2 #p by 1 vector
      A=TM[permu[n1],]%*%t(TM[permu[n2],]) #p by p matrix
      # B=matrix(0, (n-1), p)
      # for (j in (1:N_u))
      # {B[j,]=(TM[permu[j],]-TM[permu[j+1],]) } 
      # 
      # for (j in (N_u+1):(n-1))
      # {B[j,]=-(TM[permu[j],]+TM[permu[j+1],]) }   
      temp=maxi_over_pk(a0, a1, -TM, Npoint) # a function returning the maximum valuse of |g(v)| 
      #over the current P_k, corresponding to the current permu
      #print(temp)
      gmatrnew[7]=max(temp, gmatrold[7])
    }   
  } 
  gmatrnew       
}   #end the big function
#--------------------------------------------------------------------------------  
maxim_over_pk=function(b,A1,B,flag)
{ #flag is either T or F, i.e. maximization or minimization, b px1, A pxp, B (n-1)xp
  #install.packages("optiSolve")
  library(optiSolve)
  n=dim(B)[1]+1; p=dim(B)[2]  
  
  mycop=cop(f=ratiofun(Q1=matrix(0,p,p), a1=t(b), d1=0,Q2=A1,a2=rep(0,p),
                       d2=0, id=1:p, name="ratio.fun"),
            max=flag,
            lb=lbcon(-1,id=1:p),
            ub=ubcon(1,id=1:p),
            lc=lincon(A=B,dir=rep("<=",(n-1)), val=rep(0,(n-1)), name=row.names(brain)[1:27]),
            #qc=quadcon(Q=diag(rep(1,p)), a=rep(0,p), d=-1, val=0,
            #           dir="<=", id=1:p),
            qc=quadcon(Q=diag(rep(-1,p)),a=rep(0,p),d=1,val=0,
                       dir="<=", id=1:p)
  )
  
  res=solvecop(mycop,quiet=FALSE)
  
  Evaluation <- validate(mycop, res, quiet=TRUE)
  
  as.numeric(Evaluation$obj.fun)
}  
#-----------------------------------------------------------------------------  
#using auglag function do minimization
#install.packages("alabama")
#install.packages("nloptr")
library(alabama); library(nloptr)

mini_over_pk=function(x,b,A,B)
{
  x0=rep(1,length(x))
  fn=function(x){ (b%*%x)/(t(x)%*%A%*%x)}
  hin=function(x){-B%*%x}
  heq=function(x){t(x)%*%x-1}
  gr=function(x){nl.grad(x,fn)}
  hinjac=function(x){nl.jacobian(x,hin)}
  heqjac=function(x){nl.jacobian(x,heq)}
  
  res=auglag(x0,fn, gr=NULL, hin=hin, heq=heq,  localsolver="lbfgs") #slsqp, mma,lbfgs
  #print(res)
  #print(res$value)
  #print(res$iter)
  res$value
}  
#-----------------------------------------------------------
maxi_over_pk=function(a0, a1, TM, Npoint) #most foundamental optimization procedure
{
  a=runif(Npoint, min=a0, max=a1)
  temp=rep(NA, Npoint)
  for (i in 1:Npoint)
  {
  v=matrix(c(cos(a[i]),sin(a[i])), 2, 1)
  #temp[i]=(b%*%v)/(t(v)%*%A%*%v)
  proj_val=TM%*%v
  proj_val[proj_val==0]=10^{-20} #to avoid 0 in the denominator
  temp[i]=abs(median(1/proj_val))
  }
  
  v0=matrix(c(cos(a0),sin(a0)), 2, 1)
  v1=matrix(c(cos(a1),sin(a1)), 2, 1)
  
  proj_val=TM%*%v0
  proj_val[proj_val==0]=10^{-20} #to avoid 0 in the denominator
  temp0=abs(median(1/proj_val))
  proj_val=TM%*%v1
  proj_val[proj_val==0]=10^{-20} #to avoid 0 in the denominator
  temp1=abs(median(1/proj_val))

  #temp0=(b%*%v0)/(t(v0)%*%A%*%v0); temp1=(b%*%v1)/(t(v1)%*%A%*%v1)
  res=max(temp, temp0, temp1)

  res
}