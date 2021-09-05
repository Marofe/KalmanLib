function [ms0,Ps0,v,P,K] = kalman_smooth(A,H,Q,R,m0,P0,ms,Ps)
% function [hx, P, v] = kalman_smooth(A,H,Q,R,m0,P0,ms,Ps)
%
% One Prediction and Update step of Classic Kalman Smooth forward-backward 
% for the nominal system (A,H,Q,R)
%
% input:  A,H,Q,R -> matrices of the system
%         m0, P0 -> mean and covariance from Kalman Filter (forward)
%         ms, Ps -> mean and covariance from Kalman Smoother (backward)
%        
% output: ms0 -> smoothed estimate of x
%         Ps0 -> smoothed Covariance
%         v -> smoothed estimation variation
%
%
% Last Update: 27/08/2019
% Author: Marcos Rogerio Fernandes
% E-mail: eng.marofe@hotmail.com
% Site: https://marofe.github.com
% Reference: 
% Särkkä, S. (2013). Bayesian Filtering and Smoothing. Ch8.


 %% Prediction Step
    m=A*m0;
    P=A*P0*A'+Q;  
 %% Update Step
    K=P0*A'/P;
    v=K*(ms-m);
    ms0=m+v;
    Ps0=P0+K*(Ps-P)*K';
end

