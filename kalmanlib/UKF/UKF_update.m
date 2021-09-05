function [hx,P,z,v] = UKF_update(h,hx0,P0,R,y,UT)
n=size(P0,1);
m=size(R,1);
lambda=UT.alpha^2*(n+UT.kappa)-n;
X=sigmaPoints(hx0,P0,lambda);
Y=zeros(m,2*n+1);
for i=1:2*n+1
    Y(:,i)=h(X(:,i));
end
Wm=[lambda/(lambda+n) ones(1,2*n)*(1/(2*(lambda+n)))]';
Wc=[lambda/(lambda+n)+(1-UT.alpha^2+UT.beta) ones(1,2*n)*(1/(2*(lambda+n)))]';
hy=Y*Wm;
epsy=Y-hy;
epsx=X-hx0';
S=epsy*diag(Wc)*epsy'+R;
C=epsx*diag(Wc)*epsy';
%%
K=C/S;%Kalman Gain
z=y'-hy; %innovation
v=K*z; % estimation variation
hx=hx0+v'; %mean update
P=P0-K*S*K'; %covariance update
P=real(0.5*(P+P'));
end

