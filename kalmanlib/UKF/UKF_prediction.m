function [hx,P] = UKF_prediction(f,hx0,P0,Q,UT,varargin)
%Additive-noise UKF
n=size(P0,1);
lambda=UT.alpha^2*(n+UT.kappa)-n;
X=sigmaPoints(hx0,P0,lambda);
%% Propagate the Sigma-points through non-linear map
for i=1:2*n+1
    X(:,i)=f(X(:,i));
end
Wm=[lambda/(lambda+n) ones(1,2*n)*(1/(2*(lambda+n)))]';
Wc=[lambda/(lambda+n)+(1-UT.alpha^2+UT.beta) ones(1,2*n)*(1/(2*(lambda+n)))]';
hx=X*Wm;
epsx=X-hx;
P=epsx*diag(Wc)*epsx'+Q;
P=0.5*(P+P'); %compensate numeric issues
hx=hx';
end
