function [gmatrnew]= update_UF(gmatrold, u, T)
         [~, n]=size(T); kk=floor((n+1)/2);
         proj_val=T'*u; %n by 1 vector
%         if(~all(proj_val))
%           proj_val=T'*(u + 1e-5.*ones(p,1));
%           N_u=sum(proj_val(:)<0); %total number of negative projected values
%           uv=sort(proj_val); kv=sort(ones(1,n)./proj_val); %u^v and k^v
%           disp(size(gmatrold,2));
%           gmatrnew= get_UF(uv,kv, n, kk, N_u, gmatrold);
%           gmatrold=gmatrnew;
%           disp(size(gmatrold,2));
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%           proj_val=T'*(u - 1e-5.*ones(p,1));
%           N_u=sum(proj_val(:)<0); %total number of negative projected values
%           uv=sort(proj_val); kv=sort(ones(1,n)./proj_val); %u^v and k^v
%          % disp([n,kk,N_u]); disp(uv);
%           gmatrnew= get_UF(uv,kv,n, kk, N_u, gmatrold);
%          % gmatrold=gmatrnew;
%         else
           N_u=sum(proj_val(:)<0); %total number of negative projected values
           uv=sort(proj_val); kv=sort(ones(1,n)./proj_val); %u^v and k^v
%           gmatrnew= get_UF(uv,kv,n, kk, N_u, gmatrold);
           gmatrnew=gmatrold;
%         end    
%end  %end of function update_UF

%function [gmatrnew]=get_UF(uv,kv,n,kk,N_u,gmatrold)
          %disp([n,kk,N_u]);
         if (N_u==0|| N_u==n)
            if (mod(n,2)==1 && N_u==0)
              gmatrnew(1)=min([gmatrold(1), uv(kk)]);
            end    
            if (mod(n,2)==1 && N_u==n)
              gmatrnew(2)=max([gmatrold(2),uv(kk)]);
            end
            if (mod(n,2)==0 && N_u==0)
              gmatrnew(3)=min([gmatrold(3), uv(kk+1)]);
            end    
            if (mod(n,2)==0 && N_u==n)
              gmatrnew(4)=max([gmatrold(4),uv(kk+1)]);
            end           
         end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
         if ((0<N_u)&&(N_u<n))            
             if (mod(n,2)==1 && kv(kk)<0)
              gmatrnew(5)=max([gmatrold(5), uv(N_u-kk+1)]);   
             end  
             %disp([kk, N_u]); disp([size(uv,1), size(gmatrold,2)]);
             %disp(uv);
             if (mod(n,2)==1 && kv(kk)>0)
              gmatrnew(6)=min([gmatrold(6), uv(N_u+kk)]);   
             end 
             if (mod(n,2)==0 && kv(kk)>0)
              gmatrnew(7)=min([gmatrold(7), uv(N_u+kk)]);   
             end    
             if (mod(n,2)==0 && kv(kk)<0 && kv(kk+1)>0)
              gmatrnew(8)=max([gmatrold(8), uv(1)]);   
             end  
             if (mod(n,2)==0 && kv(kk)<0 && kv(kk+1)<0)
              gmatrnew(9)=max([gmatrold(9), uv(N_u-kk)]);   
             end   
         end 
%end %end of function  get_UF         
end  %end of function update_UF