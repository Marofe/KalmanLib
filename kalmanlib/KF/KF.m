function [hx,P] = KF(A,H,Q,R,hx0,P0,y,varargin)
if nargin==7
    B=zeros(size(hx0'));
    u=0;
else
    B=varargin{1};
    u=varargin{2};
end
[hx,P] = KF_prediction(hx0',P0,A,Q,B,u);
[hx,P] = KF_update(hx',P,H,R,y');
end

