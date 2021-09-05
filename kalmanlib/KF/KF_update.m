function [hx,P] = KF_update(hx0,P0,H,R,y)
K=P0*H'/(R+H*P0*H');
hx=(hx0+K*(y-H*hx0))';
P=(eye(size(P0,1))-K*H)*P0;
P=0.5*(P+P'); %compensate numeric issues
end

