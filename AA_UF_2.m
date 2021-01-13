% function of AA-UF-2 for paper Comp-RPD example 2 AA versus AA, by Y.Zuo
% 4/15/19
%%
function [UF]=AA_UF_2(Z, beta, UN)
% Z is a p by n matrix with Z(:,i)=(x_i', y_i)'(i=1,...n), 
% x_i is a p-1 vector. that is the last row of Z is y, beta is a p by 1
% vector. UN is the total number of unit diretions used in the AA 
[p, n]=size(Z);
w=zeros(p,n);
r=zeros(1,n);
T=zeros(p,n);
UFold=0;
%t=zeros((n-1),1);
%ub=ones(p,1); lb=-ones(p,1);
%gmatrinit=[1e+10,-1e+10,1e+10,-1e+10,-1e+10,1e+10,1e+10,-1e+10,-1e+10];UFold=0;

% compute T={t_i}, i=1,... n,
for i=1:n  
   w(:,i)=[1,Z(1:(p-1),i)']';  %p by 1 vector
   r(i)=Z(p,i)-dot(w(:,i),beta); %ri=yi-wi'beta
   if (r(i)==0)  %take care of zero ri
       r(i)=1e-20; 
   end 
   T(:,i)=w(:,i)./(r(i).*ones(p,1));   
end %for loop

%sample p points from {t_i} and another p points from {t_i} to get non-parallel
%normal vectors v1 and v2 of two typerplanes H1 and H2

N=nchoosek(n,p); % only allow to p=13 when n=100
KK=min([N,UN]);

%big for loop to generate the unit direction v by sample points from {t_i}
for j=1:KK 

 M_point=zeros(p,p,2); %M_permu=zeros(N,p);
 %rng(1); % could repeat the same results in rand sampling
 %rng('shuffle') gives different results every time
 sample_id_1=randperm(n,p);
 M_point(:,:,1)=T(:,sample_id_1); %p by p matrix


 sample_id_2=randperm(n,p);
 %to avoid the repeated id or points
  while( isempty( setdiff(sample_id_2, sample_id_1) ) )
      %(isequal(sort(sample_id_1), sort(sample_id_2)))
      % to avoid sort, use isempty(setdiff(sample_id_2, sample_id_1)) 
      sample_id_2= randperm(n,p);
  end %while inner loop   
  
      M_point(:,:,2)=T(:,sample_id_2);      
      % construct hyperplanes/lines Hi through the selected p points
      % and find their normal vectors vi
      pairwise_diff_1=M_point(:,1:(p-1),1)-M_point(:,p,1);    
      pairwise_diff_2=M_point(:,1:(p-1),2)-M_point(:,p,2);
      %pairwide difference of the selected p points to form p-1 vectors
      if (p==2)
         v1=[-pairwise_diff_1(2,1);pairwise_diff_1(1,1)]; v1=v1./norm(v1);
         v2=[-pairwise_diff_2(2,1);pairwise_diff_2(1,1)]; v2=v2./norm(v2);
      else  
          %pairwise_diff_1
          v1=null(pairwise_diff_1');% v1=v1(:,1); 
          %v1
          v2=null(pairwise_diff_2'); %v2=v2(:,1);
      end
      %normal vectors of the hyperplanes (or lines) Hi, i=1,2.   
      %if v2 is parallel to v1, then continue to find a v2 that is not
   while(range(v1./v2)==0) % v1 and v2 are parallel
      sample_id_2= randperm(n,p);
      M_point(:,:,2)=T(:,sample_id_2); 
      pairwise_diff_2=M_point(:,1:(p-1),2)-M_point(:,2:p,2);
      if (p==2)
         v2=[-pairwise_diff_2(2,1);pairwise_diff_2(1,1)]; v2=v2./norm(v2);
      else   
         v2=null(pairwise_diff_2'); %v2=v2(:,1)
      end
   end

  %[v1 v2]    
% the following lines are unnecessary since ui is parallel to vi      
%{   
   if (p==2)
    
   else %(p>2)

   end
  % construct hyperplanes that are perpendicular to Hi and through
  % the origin and p-2 additional points from M_point (not unique)
  % first find the points Pi=(xi',yi) on Hi that is closest to the origin
  % then by origin, Pi and the p-2 points from M_point to find the normal
  % vectors of ui, note that (xi',yi)=ci vi (vectors are parallel) and
  % ci=vi'(xi',yi)=vi'(x0i',y0i) where (x0i' y0i) on Hi, therefore
  % ci=vi'M_point(:,1,i), since vi'(Pi-Pi*)=0 for any Pi* on Hi; 
  % therefore Pi=ci vi and ci is given by the formula
  % all the above can be done in matlab with solve and subs, syms,dot,etc
  % see http://www2.math.umd.edu/~jmr/241/lines_planes.html
      c1=v1.*M_point(:,1,1); P1=c1.*v1; 
      c2=v2.*M_point(:,1,2); P2=c2.*v2;
  % normal vectors of the perpendicular hyperplanes
      u1=null([P1,M_point(:,3:p,1)]);u2=null([P2, M_point(:,3:p,2)]);
  % points listed above are the difference between them and the the origin
  % above line is invalid unless p>=3, should write as 1:p-2  
  %   u1=u1(:,1); u2=u2(:,1);
  
  %check if u1 snd u2 are parallel, if is get u2 again
 while(range(u1./u2)==0)
      sample_id_2= randperm(n,p);
      M_point(:,:,2)=T(:,sample_id_2); 
      pairwise_diff_2=M_point(:,1:(p-1),2)-M_point(:,2:p,2);
      v2=null(pairwise_diff_2); %v2=v2(:,1);
      c2=v2.*M_point(:,1,2); P2=c2.*v2;
      u2=null([P2, M_point(:,3:p,2)]); %u2=u2(:,1);
 end
%} 
  if (p==2)
    u=v1;
  else
    u=null([v1 v2]'); u=u(:,1);
   %valid only if p>2 if p=2 and rank of [v1 v2] is 2
   %then the null sapce has dimension 0.
  end
 %u=u(:,1); % size of u is supposed to be p by p-2
 %since, the rank of [u1 u2] is 2, hence the rank of null space is p-2
   v=u./norm(u);
  % v
% along each v compute the directional ufv and update UFold
 %rng('shuffle');  
 %rng(0);
% u=randn(p,1); v=u./norm(u); v=v(:,1);
 dim=size(v,2);
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
  
end  %for j=1:UN 

if (UN>KK)
     for i=1:(UN-KK)
       u=randn(p,1); v=u./norm(u);
       ufvoutput=ufvhd(v,beta,Z);
       if (ufvoutput>= UFold)
        UFold=ufvoutput;
       end 
     end
end

UF=UFold;
end % function end 
%%
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
