function [hx,P] = unscented_kf(f,h,Q,R,alpha,kappa,beta,x0,P0,y,u)
% function [hx, P, v] = unscented_kalman_filter(A,H,Q,R,x0,P0,y)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal 
% system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions      
%        
% output: hx -> final estimate
%         P -> final Covariance
%         v -> estimation variation
%
%
% Last Update: 01/04/2018
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Site: https://marofe.github.com

%% parameters
n=size(Q,1);
lambda = alpha^2*(n+kappa)-n;
%% Getting the weights
Wm=lambda/(lambda+n);
Wc=lambda/(lambda+n)+(1-alpha^2+beta);
Wm=[Wm 1/(2*(lambda+n))*ones(1,2*n)];
Wc=[Wc 1/(2*(lambda+n))*ones(1,2*n)];
%% Prediction Step
%Getting the sigma points
c=sqrt(n+lambda);
X0=sigmaPoints(x0,P0,c);
hx=zeros(n,1);
    for i=1:2*n+1
        %propagate through the dynamic model
        X(:,i)=f(X0(:,i),u);
    end
    %compute the sampled mean
        hx=X*Wm';
P=Q;
    for i=1:2*n+1
        %compute the sampled covariance
        P=P+Wc(i)*(X(:,i)-hx)*(X(:,i)-hx)';
    end
%% Update Step
p=size(R,1);
hy=zeros(p,1);
    for i=1:2*n+1
        %propagate through the measurement model
        Y(:,i)=h(X(:,i));
    end
    %compute sampled mean output
        hy=Y*Wm';
 S=R;
 C=zeros(n,p);
    for i=1:2*n+1
        %compute output covariance Pyy
        S=S+Wc(i)*(Y(:,i)-hy)*(Y(:,i)-hy)';
        %compute output covariance Pxy
        C=C+Wc(i)*(X(:,i)-hx)*(Y(:,i)-hy)';
    end
% compute the Kalman gain
    K=C/S;
% update the predictions    
    hx=hx+K*(y-hy);
    P=P-K*S*K';
end

function [X]=sigmaPoints(x0,P0,gamma)
M=sqrtm(P0);
n=size(P0,1);
X=x0;
for i=2:n+1
    X(:,i)=x0+gamma*M(:,i-1);
    X(:,i+n)=x0-gamma*M(:,i-1);
end
end
