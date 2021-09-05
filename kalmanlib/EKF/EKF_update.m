function [hx,P,z,H,v] = EKF_update(h,hx0,P0,Hx,R,y)
H=Hx(hx0); %Measurement Jacobian at hx0
hy=h(hx0);
K=P0*H'/(R+H*P0*H');%Kalman Gain
z=y-hy; %innovation
v=K*z'; % estimation variation
hx=hx0+v'; %mean update
if sum(hx<0)~=0
    warning("negative state!")
    hx(hx<0)=0;
end
P=(eye(size(P0,1))-K*H)*P0; %covariance update
P=0.5*(P+P'); %compensate numeric issues
end

