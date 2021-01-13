%by Y.Zuo on 3/6/20 for the revision of comp_PRD
function [UF]=Ex_UF_HD_no_UN(Z, beta, add) 
% given a data set Z which is p by n matrix, its first (p-1) rows are x_i
% and the pth row is y_i, n is the sample size, beta is 1 x p vector
% in parameter space (or candidate regression parameter), add is a number
% used with the combination number {n choose p} to determine the total #
% of the distinct permutations used. UN is used in Ex_UF_HD.m
[p, n]=size(Z);
w=zeros(p,n);
r=zeros(1,n);
T=zeros(p,n);
t=zeros((n-1),1);
ub=ones(p,1); lb=-ones(p,1);
gmatrinit=[1e+10,-1e+10,1e+10,-1e+10,-1e+10,1e+10,1e+10,-1e+10,-1e+10];

%generate the T matrix p by n
for i=1:n  
   w(:,i)=[1,Z(1:(p-1),i)']';  %p by 1 vector
   %disp(Z(p,i)); disp(w(:,i));
%   disp("i="); disp(i);
   %disp(beta*w(:,i));
   r(i)=Z(p,i)-beta*w(:,i); %ri=yi-beta*wi
   if (r(i)==0)
       r(i)=1e-20; 
   end  %take care of zero ri
   T(:,i)=w(:,i)./(r(i).*ones(p,1));   
end %for loop

if (p>2)
% initializing
%{
q= nchoosek(n,p); N=0;
for j=0:(p-1) 
    N=N+nchoosek(q-1,j);
end
%}
q=nchoosek(n,p);
N=q+add;
%disp([q, N]);
%N= nchoosek(n,p);
k=0; %counters
N1= N; %max([N, UN]);
M_point=zeros(p,p,N1); M_permu=zeros(N1,n);
%rng(1); % could repeat the same results in rand sampling
%rng('shuffle') gives different results every time
sample_id_1=randperm(n,p); %to randomly select p points from T
k=k+1;
M_point(:,:,k)=T(:,sample_id_1); %p by p matrix

%%
gmatrold=gmatrinit; 
%n_distinct_permu=0;
m_loops= 0; % n_distinct_permu;
K=8*(p-2); c=0; %countet of distinct permutations
%lpn=floor(UN/(K+2));
while (c<N) % (m_loops<=lpn ) %big loop
       sample_id_2=randperm(n,p);
       v1=ones(1,p); v2=v1;
       u1=v1; u2=v1;  %below ui's can be skipped, just using vi's, 3/7/20
       while (range(u2./u1)==0) %check to see if they are parallel
         while (range(v2./v1)==0)
           while (isequal(sort(sample_id_1), sort(sample_id_2)))
            sample_id_2= randperm(n,p);
           end %while inner loop   
            k=k+1;
            M_point(:,:,k)=T(:,sample_id_2);
         % construct hyperplanes Hi through the selected p points
           pairwise_diff_1=M_point(:,1:(p-1),k-1)-M_point(:,2:p,k-1);    
           pairwise_diff_2=M_point(:,1:(p-1),k)-M_point(:,2:p,k);
         %pairwide difference of the selected p points to form p vectors
           
           v1=null(pairwise_diff_1'); v2=null(pairwise_diff_2');
         %normal vectors of the hyperplanes Hi
         end %while ~=0) loop for vi   
         %take advantage of vi first  
         gmatrnew=update_UF(gmatrold, v1, T); 
         gmatrold=gmatrnew;
         gmatrnew=update_UF(gmatrold, v2, T); 
         gmatrold=gmatrnew;
         
        % construct hyperplanes that are perpendicular to Hi and through
        % the origin (not unique) and p-2 points from M_point
        % first find the points Pi=(xi',yi)'on Hi, closest to the origin
        % then use Pi and the p-2 points from M_point to find the normal
        % vectors of ui, note that (xi',yi)'=civi (vectors are parallel) &
        % vi'(xi',yi)'=vi'(x0i',y0i)' where (x0i', y0i)' on Hi, therefore
        % ci=vi'(M_point[:,1,k+i-2])/norm(vi)^2, Pi=ci vi
        % all the above can be done in matlab with solve and subs,
        % syms,dot,etc
        % see http://www2.math.umd.edu/~jmr/241/lines_planes.html
          c1=dot(v1, M_point(:,1,k-1))/dot(v1, v1); P1=c1*v1;
          c2=dot(v2, M_point(:,1,k))/dot(v2, v2); P2=c2*v2;

         M1=[zeros(p,1),P1,M_point(:,3:p,k-1)]; %points for new hyperplane
         M2=[zeros(p,1), P2, M_point(:,3:p,k)]; %points for new hyperplane
         M1_diff=M1(:,1:(p-1))-M1(:,2:p); %pointwise difference
         M2_diff=M2(:,1:(p-1))-M2(:,2:p); 
         
         u1=null(M1_diff');u2=null(M2_diff');
        % disp([u1, u2]);
       end %while loop for ui   
      
       u=null([u1'; u2']); % size of u is supposed to be p by p-2 since
       %the rank of [u1'; u2'] is 2, hence the rank of null space is p-2
       %disp(u);
       u0=u;
       u=zeros(p,K);
       
      for i=1:(p-2)
       L=8*(i-1); %u0=u(:,i)/norm(u(:,i));
      % disp(L); disp(u0);
       u(:,(L+1))=u0(:,i)+(v1+v2);
       u(:,L+2)=u0(:,i)+(v1-v2);
       u(:,L+3)=u0(:,i)-(v1+v2); u(:,L+4)=u0(:,i)-(v1-v2);
       u(:,L+5)=-u(:,L+1); u(:,L+6)=-u(:,L+2);
       u(:,L+7)=-u(:,L+3); u(:,L+8)=-u(:,L+4);
      end
      
      %disp(u);
       proj_value=zeros(n,K);
       for i=1:K
          proj_value(:,i)=T'*(u(:,i)/norm(u(:,i))); 
       end
       [B, P]=sort(proj_value, 1); 
      % disp(P(:,1:5));
       disp(["before c=:", c]);
       % update c
       for i=1:K
          if (~ismember(P(:,i)', M_permu, 'row'))
              c=c+1; 
              M_permu(c,:)=P(:,i)';
          end
       end 
       m_loops=m_loops+1;
       disp(["after c=:", c]);
     %update UF or g(u) via corollary 2.1
      for i=1:K
         f=repmat(B(floor((n+1)/2),i), p,1); M=T(:,P(:,i)');
         A=M(:,1:(n-1))-M(:,2:n);
        
         options = optimoptions('linprog','Display','none');
         u=linprog(f, A', t, [], [], lb,ub, options); u=u/norm(u);
         gmatrnew=update_UF(gmatrold, u, T); 
         gmatrold=gmatrnew;
     end % for loop
     disp(["permu#, loop#", [c, m_loops]]);
end %while (c<N) loop    
%%
 diff=gmatrold-gmatrinit;
 fm=gmatrold(diff~=0);
 g=abs(fm); m=size(g,2); 
 UF=0;
 for jj=1:m
  UF=max([UF, 1/g(jj)]);
 end 
end %if(p>2)
%%
end %end of function Ex_UF_HD programm