function [Xi,hx,P,z,v] = PF_update(h,Xi,hx0,R,y)
Np=size(Xi,2);
m=size(R,1);
Yi=zeros(m,Np);
W=zeros(Np,1);
iR=eye(m)/R;
for i=1:Np
    Yi(:,i)=h(Xi(:,i));
    z=y-Yi(:,i);
    W(i)=exp(-0.5*z'*iR*z);
end
W=W/sum(W);
hy=mean(Yi,2); %measurement estimate
z=y-hy;
hx=Xi*W; 
v=hx-hx0; % mean estimation variation
errX=Xi-hx;
P=errX*diag(W)*errX'; %covariance update
P=0.5*(P+P'); %compensate numeric issues
[Xi,W]=resampling(Xi,W);
end

