%%
function [m]=update_m(u,T)
    proj_val=T'*u; %n by 1 vector
    % get sorted values and the permutation
   % disp(proj_val);
    kk=floor((size(T,2)+1)/2);
    %disp(kk);
    [proj_val_sort, permu]=sort(proj_val);
    %disp(proj_val_sort);
    m=zeros(1,2);
    if (proj_val_sort(kk)==proj_val_sort(kk-1)) 
       m(:,1)=permu(kk); m(:,2)=permu(kk-1); 
    elseif (proj_val_sort(kk)==proj_val_sort(kk+1)) 
       m(:,1)=permu(kk); m(:,2)=permu(kk+1);  
    else
      %  disp(kk); disp(permu);
       m(:,1)=permu(kk); m(:,2)=permu(kk); %unique case
     %  disp(permu(kk));
    end
end %end of function update_m

