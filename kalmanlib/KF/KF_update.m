function [hx,P,z,v] = KF_update(hx0,P0,H,R,y)
hy=H*hx0;
K=P0*H'/(R+H*P0*H');%Kalman Gain
z=y-hy; %innovation
v=K*z; % estimation variation
hx=hx0+v; %mean update
P=(eye(size(P0,1))-K*H)*P0; %covariance update
P=0.5*(P+P'); %compensate numeric issues
end

