function [hx,P,z] = hinf_filter(A,H,G,Q,R,x0,P0,y,L,gamma)
% function [hx, P, z] = hinf_filter(A,H,G,Q,R,x0,P0,y,L,gamma)
% System: x(k+1)=Ax(k)+Gw(k)
%         y(k)=Hx(k)+e(k)
%         w(k)~N(0,Q) -> process noise
%         e(k)~N(0,I) -> measurement noise
%
% One Prediction and Update step of Hinf Filter for the nominal 
% system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the nominal model
%         x0, P0 -> initial conditions      
%         L -> Output matrix
%         G -> Disturbance matrix
%        
% output: hx -> final estimate
%         P -> final Covariance
%         z -> output (z=L*hx)
%
% Implementation based on the paper
% Shaked, U., & Theodor, Y. (1992). Hinf Optimal Estimation: A tutorial. 
% 31st Conference on Decision and Control, 2278–2286. Tucson, Arizona: IEEE.
%
%
% Author: Marcos Rogerio Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 25/09/2018

n=size(A,1);
p=size(y,1);
%%
if min(eig(inv(P0)-L'*L/gamma^2))>0
    %% Existence condition satisfied
    tildeP=(inv(P0)-L'*L/gamma^2)\eye(n);
    hx=A*x0+A*tildeP*H'/(eye(p)+H*tildeP*H')*(y-H*x0);
    Re=blkdiag(eye(p),-gamma^2*eye(n))+[H;L]*P0*[H' L'];
    K=A*P0*[H' L'];
    P=A*P0*A'+G*Q*G-K/Re*K';
    z=L*hx;
else
    hx=inf;
    P=inf;
    z=inf;
end
end