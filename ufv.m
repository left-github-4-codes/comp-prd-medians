function [UFV]=ufv(v, beta, X)
n=size(X,2);
w=horzcat(ones(n,1),X(1,:)');
Nvect=X(2,:)'-w*beta; %numinator of UF
Dvect=w*v;   % denominator of UF  
Dvect(Dvect==0)=1e-20; %take care of denominaor zero
quotient=Nvect./Dvect;
UFV=abs(median(quotient));
end