function [hx,P] = BDU_filter(F,H,M,Ef,Eg,G,Q,R,x0,P0,y)
% function [hx, P] = BDU_filter(A,H,Q,R,x0,P0,y)
%
% One Prediction and Update step of BDU Filter for the nominal 
%system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions     
%         M,Ef,Eg,G -> uncertainties parameters
%        
% output: hx -> final estimate
%         P -> final Covariance
%
% Implementation based on the paper
% Sayed, A. H. (2001). A framework for state space estimation with 
% uncertain models. Ieee Transactions on Automatic Control, 46(7), 1–15. 
% Retrieved from https://ieeexplore.ieee.org/document/935054/
%
% Author: Marcos Rogerio Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 25/10/2018

n=size(F,1);
p=size(y,1);
%% sayed paremeters
% b=y-H*F*x0;
% A=H*[F G];
% Ql=blkdiag(P0\eye(n),Q\eye(n));
% W=R\eye(p);
% Hl=H*M;
% Ea=[Ef Eg];
% Eb=-Ef*x0;
%lopt=BDU_opt(A,Ql,W,b,Ea,Eb,Hl);
%%
l=norm(M'*H'/R*H*M);
lopt=1.5*l; %approx opt gamma
if (lopt==0)
  hQ=Q;
  hG=G;
  hR=R;
  hF=F;
  hP=P0;
else
hQ=Q;
hG=G;
hP=(P0\eye(n)+lopt*Ef'*Ef)\eye(n);
hF=F*(eye(n)-lopt*hP*Ef'*Ef);
hR=R-1/lopt*H*M*M'*H';
end
 %% prediction
    hx=hF*x0;
    P=hF*hP*hF'+hG*hQ*hG';  
    %% update
    P=P-P*H'/(hR+H*P*H')*H*P;
    K=P*H'/(hR);
    hx=hx+K*(y-H*hx);
end