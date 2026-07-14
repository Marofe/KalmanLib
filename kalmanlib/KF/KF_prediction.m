function [hx,P] = KF_prediction(hx0,u,P0,A,B,Q,varargin)
hx=A*hx0+B*u;
P=A*P0*A'+Q;
P=0.5*(P+P'); %compensate numeric issues
end

