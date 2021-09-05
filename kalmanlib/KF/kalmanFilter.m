function [x,P] = kalmanFilter(x,P,y,sys)
% function [hx, P] = kalman_filter(x0,P0,y,sys)
%
% One Prediction and Update step of Classic Kalman Filter for the nominal
% system (A,H,Q,R) with x0,P0 initial condition and y measurement
%
% input:  x0, P0        -> initial conditions
%         sys={A,H,Q,R} -> matrices of the system as a struct
%
%
% output: hx -> posterior estimate of x(k) given y(1:k)
%         P  -> posterior Covariance

%% Prediction Step (time-update)
x=sys.A*x;
P=sys.A*P*sys.A'+sys.Q;
%% Update Step (measurement-update)
P=P-P*sys.H'/(sys.R+sys.H*P*sys.H')*sys.H*P;
K=P*sys.H'/sys.R;
x=x+K*(y-sys.H*x);
end
