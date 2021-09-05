function [hx,P,A] = EKF_prediction(f,hx0,P0,Ax,Q,varargin)
A=Ax(hx0); %Jacobian at hx0
hx=f(hx0);
P=A*P0*A'+Q;
P=0.5*(P+P'); %compensate numeric issues
end

