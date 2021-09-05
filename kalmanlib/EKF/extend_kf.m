function [hx,P,v,Ppred,K,hy] = extend_kf(f,h,Ai,Hi,Q,R,x0,P0,y,u)
% function [hx, P, v] = extend_kalman_filter(f,h,A,H,Q,R,x0,P0,y,u)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal
%system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions
%
% output: hx -> final estimate
%         P -> final Covariance

%
% Author: Marcos Rogério Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 07/05/2019

%% prediction step
hx=f(x0,u);
%Jacobianos
A=Ai(x0,u);
H=Hi(hx);
Ppred=A*P0*A'+Q;
%% update step
P=inv(inv(Ppred)+H'*inv(R)*H);
K=P*H'*inv(R);
hy=h(hx);
v=K*(y-hy);
hx=hx+v;
end