%% Ex_UF_2Plus_2 in Compu-RPD, 1/3/2019 by Y. Zuo
function [UF]=Ex_UF_2plus_2(Z, beta, UN) 
% given a data set Z which is p by n matrix, its first (p-1) rows are x_i
% and the pth row is y_i, n is the sample size, beta is 1 x p vector
% in parameter space (or candidate regression parameter) UN is used to stop
% while loop and for other purpose when to compare with approximate methods
[p, n]=size(Z);
w=zeros(p,n);
r=zeros(1,n);
T=zeros(p,n);
t=zeros((n-1),1);
ub=ones(p,1); lb=-ones(p,1);
gmatrinit=[1e+10,-1e+10,1e+10,-1e+10,-1e+10,1e+10,1e+10,-1e+10,-1e+10];

for i=1:n  
   w(:,i)=[1,Z(1:(p-1),i)']';  %p by 1 vector
   r(i)=Z(p,i)-dot(w(:,i),beta); %ri=yi-wi'beta
   if (r(i)==0)
       r(i)=1e-20; 
   end  %take care of zero ri
   %r(r==0)=1e-20; %take care of zero ri
   T(:,i)=w(:,i)./(r(i).*ones(p,1));   
end %for loop
%disp(T);
%%
if (p==2)
  % compute the pairwise difference as the vectors/dirs connecting two points
  % and record the corresponding id of the pair
  
   N=n*(n-1)/2; %right upper triangle total
   pairwisedir=zeros(p, N); 
   pairwise_id=zeros(N,2); %preassign to speed up
  %dir=[]; pairwise_id=[]; could use this and horzcat and vertcat but slow
   for i=1:(n-1)
      M_diff=T-T(:,i)*ones(1, n); %difference to ith point, could use repmat
      s=(i-1)*n-i*(i-1)/2+1; t=i*n-i*(i+1)/2;%positions of beginning and end
      pairwisedir(:,s:t) =M_diff(:,(i+1):n);  %(n-i) columns filled
      %dir=horzcat(dir, M_diff(:,(i+1):n)); 
      pairwise_id(s:t,:)=[i*ones((n-i),1), ((i+1):n)']; %(n-i) rows added
      %pairwise_id=vertcat(pairwise_id, [i*ones((n-i),1), ((i+1):n)']);
   end %for loop entire loop corresponds to EA_UF2 step (2)  
    
   perpend_slop=-pairwisedir(1,:)./pairwisedir(2,:);
   alpha=atan(perpend_slop)+(perpend_slop<0)*pi;%perpendicular directions
   [alpha_sort, ~]=sort(alpha);% could use unique to delete the repeated
   % angles i.e. in the case of parallet line segments
   % initial direction 
   u0=[1 0]'; %m=zeros(1,2);  %equivalent to set alpha0=0;
   [m]=update_m(u0,T); %select the two middle subscripts of projected values
  % disp(m);disp(T);
   %%
   %update UF according to Corollary 2.1, call a function update_UF
   ll=0; %counter of function update_UF called time
   gmatrnew=update_UF(gmatrinit, u0, T); %compute the UF along the u0 first.
   gmatrold=gmatrnew;
   %disp(size(gmatrold,2));
   
   for k=1:n
      uu=[-T(2,k); T(1,k)];
      v=uu./norm(uu);
      epsil=1e-7;
      v1=v+epsil.*ones(p,1);
      v1=v1./norm(v1);
      gmatrnew=update_UF(gmatrold, v1, T);
      gmatrold=gmatrnew;   
      v2=v-epsil.*ones(p,1);
      v2=v2./norm(v2);
      gmatrnew=update_UF(gmatrold, v2, T);
      gmatrold=gmatrnew;
%      disp(size(gmatrold,2));
   end    
 %} % above 2*n unit directions are used 
   ll=ll+1; % this is the counter of unit vetors used so far (u0 counted)
   counter=0; 
   while(counter<N)
    counter=counter+1;   
    uk=[cos(alpha_sort(counter)),sin(alpha_sort(counter))]';   
   %disp(uk);
   % skipping some directions that are not median sequence 
     if (sum(ismember(m',pairwise_id(counter,:)', 'row'))==0)
         counter=counter+1; 
     else 
         m=pairwise_id(counter,:); 
         gmatrnew=update_UF(gmatrold, uk, T);
         gmatrold=gmatrnew;
         ll=ll+1;
     end   
   end %while loop
 % disp(ll);
%{ 
   dir=differencedir(Z);     
   for j=1:N 
     u=[-dir(2,j);dir(1,j)];
     v=u./norm(u);
     %disp(v);
     gmatrnew=update_UF(gmatrold, v, T);
     gmatrold=gmatrnew;
   end
 %}
% v=[0.8944; 0.4472]+0.00001.* [1;1];
% gmatrnew=update_UF(gmatrold, v, T);
% gmatrold=gmatrnew;
% disp(ones(1,9)./gmatrold);
 
 diff=gmatrold-gmatrinit;
 fm=gmatrold(diff~=0);
 g=abs(fm); m=size(g,2); 
 UF=0;
 for jj=1:m
  UF=max([UF, 1/g(jj)]);
 end 
end %if(p==2)
%%
%pairwise mm(mm-1)/2 unit directions induced from the pairwise differences
 function [dir]=differencedir(X)
      mm=size(X,2);
      NN=mm*(mm-1)/2;
      pairwisediff=zeros(2, NN);  %preassign to speed up
    for ii=1:(mm-1)
       DiffM=X-X(:,ii)*ones(1, mm); %difference to ith point, could use repmat
       aa=(ii-1)*mm-ii*(ii-1)/2+1; bb=ii*mm-ii*(ii+1)/2;%positions of beginning and end
       pairwisediff(:,aa:bb) =DiffM(:,(ii+1):mm);%(n-i)
    end
     dir=pairwisediff; %./norm(pairwisediff);
     %dir=normc(pairwisediff);
 end
%%
%{
function [m]=update_m(u,T)
    proj_val=T'*u;kk=floor((size(T,2)+1)/2);%n by 1 vector
    % get sorted values and the permutation
    [proj_val_sort, permu]=sort(proj_val);
    if (proj_val_sort(kk)==proj_val_sort(kk-1)) 
       m(:,1)=permu(kk); m(:,2)=permu(kk-1); 
    elseif (proj_val_sort(kk)==proj_val_sort(kk+1)) 
       m(:,1)=permu(kk); m(:,2)=permu(kk+1);  
    else
       m(:,1)=permu(kk); m(:,2)=permu(kk); %unique case
    end
end %end of function update_m
%}
%{
%%
function [g_0v, g_nv, g_0nv1, g_0nv2]= update_UF(u, T)
         g_0v=1e+10;
         g_nv=-1e+10;
         g_0nv1=-1e+10 ;
         g_0nv2=1e+10 ;
         [p, n]=size(T); k=floor((n+1)/2);
         proj_val=T'*u; %n by 1 vector
         N_u=sum(proj_val(:)<0); %total number of negative projected values
         uv=sort(proj_val); kv=sort(ones(1,n)./proj_val); %u^v and k^v
         
         if (N_u==0|| N_u==n)
            if N_u==0
              g_0v=abs( 1/min([g_0v,uv(k)]) );  
            end    
            if N_u==n
              g_nv=abs( 1/max([g_nv,uv(k)]) ); 
            end   
         end    
         if ((0<N_u)&&(N_u<n))            
             if kv(k)<0
               g_0nv1=  abs( 1/max([g_0nv1,uv(N_u-k)]) );
             end             
             if (kv(k)>0) 
               g_0nv2=  abs( 1/min([g_0nv2, uv(N_u+k)]));
             end           
         end    
end  %end of function   update_UF
%% 
%}
if (p>2)
% initializing
%{
q= nchoosek(n,p); N=0;
for j=0:(p-1) 
    N=N+nchoosek(q-1,j);
end
%}
%q= nchoosek(n,p);
%N=q+0;
N= nchoosek(n,p);
k=0; %counters
N1=max([N, UN]);
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
lpn=floor(UN/(K+2));
%while(c<q+add)
while (m_loops<=lpn ) %big loop
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
         end %while loop for vi   
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
        
       end %while loop for ui   
      
       u=null([u1'; u2']); % size of u is supposed to be p by p-2 since
       %the rank of [u1'; u2'] is 2, hence the rank of null space is p-2
       u0=u;
       u=zeros(p,K);
       
      for i=1:(p-2)
       L=8*(i-1); %u0=u(:,i)/norm(u(:,i));
       u(:,(L+1))=u0(:,i)+(v1+v2);
       u(:,L+2)=u0(:,i)+(v1-v2);
       u(:,L+3)=u0(:,i)-(v1+v2); u(:,L+4)=u0(:,i)-(v1-v2);
       u(:,L+5)=-u(:,L+1); u(:,L+6)=-u(:,L+2);
       u(:,L+7)=-u(:,L+3); u(:,L+8)=-u(:,L+4);
      end
      
       proj_value=zeros(n,K);
       for i=1:K
          temp= u(:,i)/norm(u(:,i));
          proj_value(:,i)=T'*temp;
          gmatrnew=update_UF(gmatrold, temp, T); 
          gmatrold=gmatrnew;
       end
       [B, P]=sort(proj_value, 1); 
      % disp(["before c=:", c]);
       % update c
       for i=1:K
          if (~ismember(P(:,i)', M_permu, 'row'))
              c=c+1; 
              M_permu(c,:)=P(:,i)';
              %distinct permutation obtained
              sort_proj_val=B(:,i);
              N_u=sum(sort_proj_val(:)<0); kk=floor((n+1)/2);
              uv=B(:,i);kv=ones(1,n)./uv;
              % generate f'x in the linear programming min_x f'x  
              if(N_u==0 &&mod(n,2)==1)
                 f=T(:,P(kk,i));
              end    
              if(N_u==0 && mod(n,2)==0)
                f=T(:,P(kk+1,i));
              end
              if((0<N_u)&& (N_u<n))
                 if(mod(n,2)==1 && kv(kk)>0)
                  f=T(:,P(N_u+kk, i));                 
                 end
                 if(mod(n,2)==0 && kv(kk)>0)
                  f=T(:,P(N_u+kk,i));
                 end 
              end    
              % take care of max cases in coroll.2.1, convert it to min
              if (N_u==n && mod(n,2)==1)
                f=-T(:,P(kk,i));               
              end
              if (N_u==n && mod(n,2)==0)
                f=-T(:,P(kk+1,i));               
              end
              if (0<N_u)&&(N_u<n)
                 if (mod(n,2)==1 && kv(kk)<0)
                   f=-T(:,P(N_u-kk+1,i));
                 end
                 if (mod(n,2)==0 && kv(kk)<0 && kv(kk+1)>0)
                   f=-T(:,P(1,i));                    
                 end
                 if (mod(n,2)==0 && kv(kk+1)<0)
                   f=-T(:,P(N_u-kk,i));
                 end    
              end
              M=T(:,P(:,i)'); A=M(:,1:(n-1))-M(:,2:n);
     
              options = optimoptions('linprog','Display','none');  
              u=linprog(f, A', t, [], [], lb,ub, options); u=u/norm(u);
              gmatrnew=update_UF(gmatrold, u, T); 
              gmatrold=gmatrnew; 
              
          end % distinct permutation if
       end % for i=1:K loop 
       m_loops=m_loops+1;
       %disp(["after c=:", c]);
     %update UF or g(u) via corollary 2.1
%{     
      for i=1:K
         f=repmat(B(floor((n+1)/2),i), p,1); M=T(:,P(:,i)');
         A=M(:,1:(n-1))-M(:,2:n);
        
         options = optimoptions('linprog','Display','none');
         u=linprog(f, A', t, [], [], lb,ub, options); u=u/norm(u);
         gmatrnew=update_UF(gmatrold, u, T); 
         gmatrold=gmatrnew;
     end % for loop
 %}    
     disp(["permu#, loop#", [c, m_loops]]);
end %while (c<N) loop      
%%
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
end %end of function Ex_UF_2plus_2 programm