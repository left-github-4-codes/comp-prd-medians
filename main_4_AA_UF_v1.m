%modified to add AA_UF_3 on 4/28/19
function [results]=main_4_AA_UF_v1(p,n,UN,R)
%p=2; n=100; UN=1000; R=1000;
beta=zeros(p,1);
tAA1=0;tAA2=0; tAA3=0;
UF1=zeros(1,R); UF2=zeros(1,R); UF3=zeros(1,R);

for i=1:R
    Z=randn(p,n);
    tic;
    UF1(1,i)=AA_UF_1(Z,beta,UN);
    td=toc;
    tAA1=td+tAA1;
    
    tic;
    UF2(1,i)=AA_UF_2(Z,beta,UN);
    td=toc;
    tAA2=td+tAA2; 
    
    tic;
    UF3(1,i)=AA_UF_3(Z,beta,UN);
    td=toc;
    tAA3=td+tAA3;         
end

%compute the mean and standard deviation of UF's
results=[mean(UF1,2) std(UF1,0,2) tAA1; mean(UF2,2) std(UF2,0,2) tAA2; mean(UF3,2) std(UF3,0,2) tAA3];

end


