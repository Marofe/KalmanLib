function [hx,P] = EKS(f,Ax,Q,hx0,P0,x1,P1)
[hx,P,Ak] = EKF_prediction(f,hx0,P0,Ax,Q);
G=P0*Ak'/P;
hx=hx0+(x1-hx)*G';
P=P0+G*(P1-P)*G';
end

