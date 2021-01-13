% function for AA-UF-1 in paper Compu_PRD for example 2 for AA-UF-1 vs
% AA-UF-2
%%
function [UF]=AA_UF_1(Z, beta, UN)
% Z is a p by n matrix with Z(:,i)=(x_i', y_i)'(i=1,...n), 
% x_i is a p-1 vector. that is the last row of Z is y, beta is a p by 1
% vector. UN is the total number of unit diretions used in the AA 

[p,~]=size(Z);
UFold=0;

% generate UN unit directions v, along each v compute the directional ufv 
for j=1:UN
 %rng('shuffle');  
 %rng(0);
 u=randn(p,1); v=u./norm(u);
 ufvoutput=ufvhd(v,beta,Z);
 if (ufvoutput>= UFold)
    UFold=ufvoutput;
 end
end 
 UF=UFold;
end

function [UFV]=ufvhd(v, beta, X) %X is p by n, beta and v are p by 1 vector.
[p,n]=size(X);
w=horzcat(ones(n,1),X(1:(p-1),:)');  % w is a n by p matrix
Nvect=X(p,:)'-w*beta; %numinator of UF, n by 1 vector, residual
Dvect=w*v; % denominator of UF , n by 1 vector
Dvect(Dvect==0)=1e-20; %take care of denominaor zero
quotient=Nvect./Dvect;
UFV=abs(median(quotient));
end
