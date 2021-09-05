function [hx,P,v] = kalman_filter(A,H,Q,R,x0,P0,y)
% function [hx, P, v] = kalman_filter(A,H,Q,R,x0,P0,y)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal 
% system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  A,H,Q,R -> matrices of the system
%         x0, P0 -> initial conditions      
%        
% output: hx -> final estimate (a posteriori)
%         P -> final Covariance
%         v -> estimation variation
%
%
% Last Update: 01/04/2018
% Author: Marcos Rogerio Fernandes
% E-mail: eng.marofe@hotmail.com
% Personal Site: https://marofe.github.io


 %% Prediction Step
    hx=A*x0;
    P=A*P0*A'+Q;  
 %% Update Step
    P=P-P*H'/(R+H*P*H')*H*P;
    K=P*H'/R;
    v=K*(y-H*hx);
    hx=hx+v;
end
