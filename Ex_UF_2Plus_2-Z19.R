Ex_UF_2plus_Z19<-function(Z, beta) # convert from matlab code
{
# given a data set Z which is n by p matrix, its first (p-1) columns are x_i
# and the pth column is y_i, n is the sample size, beta is 1 x p vector
# in parameter space (or candidate regression parameter)
#
  n=dim(Z)[1]; p=dim(Z)[2]
  gmatrinit=c(1e+10,-1e+10,1e+10,-1e+10,-1e+10,1e+10,1e+10,-1e+10,-1e+10);

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
  T=TM
#--------------------------------------------------------
if (p==2)
 {
 # compute the pairwise difference as the vectors/dirs connecting two points
 # and record the corresponding id of the pair
  TM=t(TM) #p by n matrix  
  N=n*(n-1)/2 #right upper triangle total
  pairwisedir=matrix(0,p,N) 
  pairwise_id=matrix(0,N,2) #preassign to speed up
  #dir=[]; pairwise_id=[]; could use this and horzcat and vertcat but slow
  for (i in 1:(n-1))
  {
  M_diff=TM-matrix(TM[,i],p, n) #difference to ith point, could use repmat
  s=(i-1)*n-i*(i-1)/2+1; t=i*n-i*(i+1)/2#positions of beginning and end
  pairwisedir[,s:t] =M_diff[,(i+1):n]  #(n-i) columns filled
  pairwise_id[s:t,]=cbind(matrix(i, (n-i),1), matrix((i+1):n, (n-i),1))
  } #end #for loop entire loop corresponds to EA_UF2 step (2)  
#--------------------------------------------------------------------------------
  perpend_slop=-pairwisedir[1,]/pairwisedir[2,]  #slope of N perpendicular directions
  alpha=atan(perpend_slop)+(perpend_slop<0)*pi # all alpha is in-between 0 and pi
  alpha_sort=sort(alpha) # could use unique to delete the repeated
  # angles i.e. in the case of parallet line segments
  # initial direction 
  u0=matrix(c(1, 0), 2, 1) #m=zeros(1,2);  %equivalent to set alpha0=0;
  #TM=t(TM) #n by p matrix    
  m=update_m(u0,T) #select the two middle subscripts of projected values
      #update UF according to Corollary 2.1, call a function update_UF
      ll=0 #counter of function update_UF called time
      gmatrnew=update_UF(gmatrinit, u0, T) #compute the UF along the u0 first.
      gmatrold=gmatrnew

#----along the directions of perpendicular to T_i (or ti)      
  for(k in 1:n)
  { #TM is p by n matrix now   
  uu=c(-TM[2, k], TM[1, k])
  v=uu/norm(uu)
  epsil=1e-7
  v1=v+epsil
  v1=v1./norm(v1)
  gmatrnew=update_UF(gmatrold, v1, T)
  gmatrold=gmatrnew   
  v2=v-epsil
  v2=v2./norm(v2)
  gmatrnew=update_UF(gmatrold, v2, T)
  gmatrold=gmatrnew 
  }  # end of for loop
  #above 2*n unit directions are used
#-------------------------------------------------------------      
  ll=ll+1 # this is the counter of unit vetors used so far (u0 is counted)
  counter=0 
  while(counter<N)
   { counter=counter+1   
    uk=matrix(c(cos(alpha_sort(counter)),sin(alpha_sort(counter))),2,1)   
      # skipping some directions that are not median sequence 
        if (sum(c(m) %in% pairwise_id[counter,])==0) #there is some issue with
           counter=counter+1                         # this skipping rule
        else 
          { m=pairwise_id[counter,]
            gmatrnew=update_UF(gmatrold, uk, T)
            gmatrold=gmatrnew
            ll=ll+1
          }   
   } #while loop
 #-----------------------------------------------------------------------
    diff=gmatrold-gmatrinit;
    fm=gmatrold(diff!=0);
    g=abs(fm); m=length(g); 
    UF=0;
    for(jj in 1:m)
      {
       UF=max(UF, 1/g[jj]);
      }
 } #%if(p==2)
#---------------end of if p=2-----------------------------------------

  UF  
} #end function Ex_UF_2plus_2  
#--------------------end of the main function-----------------------------
###################################################################
update_m<-function(u,T)
{
  m=matrix(0,1,2)
  proj_val=T%*%u #n by 1 vector
  kk=floor((dim(T)[1]+1)/2)
    # get sorted values and the permutation
    proj_val_sort=sort(proj_val); permu=order(proj_val)
    if (proj_val_sort[kk]==proj_val_sort[kk-1]) 
     {  m[,1]=permu[kk]; m[,2]=permu[kk-1]} 
    else
    if (proj_val_sort[kk]==proj_val_sort[kk+1]) 
       { m[,1]=permu[kk]; m[,2]=permu[kk+1]}  
       else
       { m[,1]=permu[kk]; m[,2]=permu[kk]} #unique case,i.e. u is not a perpendicular
                                            # line segment direction
  m                                           
}#end function update_m
#########################################################################
update_UF<-function(gmatriold,u, T)
{
   n=dim(T)[1]; kk=floor((n+1)/2);
   proj_val=T%*%u; #%n by 1 vector
           N_u=sum(proj_val[]<0); #%total number of negative projected values
           uv=sort(proj_val); kv=sort(rep(1,n)/proj_val); #%u^v and k^v

           gmatrnew=gmatrold;

         if (N_u==0|| N_u==n)
          {  if (mod(n,2)==1 && N_u==0)
              {gmatrnew[1]=min(gmatrold[1], uv[kk]);
              }    
             if (mod(n,2)==1 && N_u==n)
              {gmatrnew[2]=max(gmatrold[2], uv[kk]);
              }
             if (mod(n,2)==0 && N_u==0)
              {gmatrnew[3]=min(gmatrold[3], uv[kk+1]);
              }    
             if (mod(n,2)==0 && N_u==n)
              {gmatrnew[4]=max(gmatrold[4], uv[kk+1]);
              }           
          }   

         if ((0<N_u)&&(N_u<n))            
          {   if (mod(n,2)==1 && kv[kk]<0)
              {gmatrnew[5]=max(gmatrold[5], uv[N_u-kk+1]);   
              }  
  
             if (mod(n,2)==1 && kv[kk]>0)
              {gmatrnew[6]=min(gmatrold[6], uv[N_u+kk]);   
              } 
             if (mod(n,2)==0 && kv[kk]>0)
              {gmatrnew[7]=min(gmatrold[7], uv[N_u+kk]);   
              }    
             if (mod(n,2)==0 && kv[kk]<0 && kv[kk+1]>0)
              {gmatrnew[8]=max(gmatrold[8], uv[1]);   
              }  
             if (mod(n,2)==0 && kv[kk]<0 && kv[kk+1]<0)
              {gmatrnew[9]=max(gmatrold[9], uv[N_u-kk]);   
              }   
         } 
   gmatrinew       
}#end function of function update_UF
  
  