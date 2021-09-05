function X = sigmaPoints(x0,P0,lambda)
%SIGMAPOINTS Summary of this function goes here
%   Detailed explanation goes here
if size(x0,1)==1
    x0=x0';
end
P0=0.5*(P0+P0');
minEig=min(real(eig(P0)));
if minEig < 0
    P0=P0+eye(size(P0,1))*abs(minEig*2);
end
S=chol(P0)';
t=sqrt(lambda+size(P0,1));
X=[x0 x0+t*S x0-t*S];
end

