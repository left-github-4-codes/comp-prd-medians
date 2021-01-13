function [results]=main_4_AA_VS_Ex(p,n,UN,R)
 %R is the replication number
 % UN is the number of directions used in AA or UN/10 is the while loop # in
 % Ex
 beta=ones(p,1);
 ExUF=zeros(1,R);
 AA1UF=zeros(1,R);
 AA2UF=zeros(1,R);
 AA3UF=zeros(1,R);
 Z=randn(p,n);
 for i=1:R
     ExUF(1,i)= EA_UFHD_final(Z,beta', UN); 
     %ExUF(1,i)=Ex_UF_HD(Z,beta', UN);
     AA1UF(1,i)=AA_UF_1(Z,beta, UN);
     AA2UF(1,i)=AA_UF_2(Z,beta, UN);
     AA3UF(1,i)=AA_UF_3(Z,beta, UN);
     %disp(["i=:", i]);
 end % end for loop   
 temp=[ExUF(1:R)', AA1UF(1:R)', AA2UF(1:R)', AA3UF(1:R)'];
 temp1=temp-ExUF(1:R)';
 res1=mean(temp, 1);
 res2=min(temp, [], 1);
 res3=max(temp,[], 1);
 res4=res3-max(ExUF);
 t1=(temp1(:,1)<0); t2=(temp1(:,2)<0);
 t3=(temp1(:,3)<0); t4=(temp1(:,4)<0);
 res5=sum([t1, t2, t3, t4], 1);
% disp([res1', res2', res3',res3', res4']); 
 results=[res1', res2', res3', res4', res5'];
 disp(["n-p-UN-R=:", [n, p, UN, R]]);
end %end of big function