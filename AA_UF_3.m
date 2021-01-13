% AA_UF_3.m by Y.Zuo on 4/28/19 added to compu_PRD using the directions
% normal to hyperplnes determined by p points
%%
function [UF]=AA_UF_3(Z, beta, UN)
% Z is a p by n matrix with Z(:,i)=(x_i', y_i)'(i=1,...n), 
% x_i is a p-1 vector. that is the last row of Z is y, beta is a p by 1
% vector. UN is the total number of unit diretions used in the AA 
[p, n]=size(Z);
w=zeros(p,n);
r=zeros(1,n);
T=zeros(p,n);
UFold=0;

% compute T={t_i}, i=1,... n,
for i=1:n  
   w(:,i)=[1,Z(1:(p-1),i)']';  %p by 1 vector
   r(i)=Z(p,i)-dot(w(:,i),beta); %ri=yi-wi'beta
   if (r(i)==0)  %take care of zero ri
       r(i)=1e-20; 
   end 
   T(:,i)=w(:,i)./(r(i).*ones(p,1));   
end %for loop

% a big loop to update UF
for j=1:UN
    
     M_point=zeros(p,p);
  
     sample_id_1=randperm(n,p);
     M_point(:,:)=T(:,sample_id_1); %p by p matrix
     
     %find the normal vector to the hyperplane formed by p points
     pairwise_diff=M_point(:,1:(p-1))-M_point(:,p);  
     
     if (p==2)
         v1=[-pairwise_diff(2,1);pairwise_diff(1,1)]; v1=v1./norm(v1);
     else  
         v1=null(pairwise_diff');% v1=v1(:,1); 
     end
     
     v=v1./norm(v1);
     dim=size(v,2);
% along v update UF     
if (dim==1) 
      ufvoutput=ufvhd(v,beta,Z);
   if (ufvoutput>= UFold)
      UFold=ufvoutput;
   end
else
   for jj=1:dim
      ufvoutput=ufvhd(v(:,jj),beta,Z);
         if (ufvoutput>= UFold)
           UFold=ufvoutput;
         end         
   end    
end
 
end %for UN loop

UF=UFold;
end %function end

%%
function [UFV]=ufvhd(v, beta, X) %X is p by n, beta and v are p by 1 vector.
[p,n]=size(X);
w=horzcat(ones(n,1),X(1:(p-1),:)');  % w is a n by p matrix
Nvect=X(p,:)'-w*beta; %numinator of UF, n by 1 vector, residual
%size(w)
%size(v)
Dvect=w*v; % denominator of UF , n by 1 vector
Dvect(Dvect==0)=1e-20; %take care of denominaor zero
quotient=Nvect./Dvect;
UFV=abs(median(quotient));
end
%%
