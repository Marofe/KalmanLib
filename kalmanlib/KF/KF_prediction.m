function [hx,P] = KF_prediction(hx0,P0,A,Q,varargin)
if nargin==4
    B=zeros(size(hx0));
    u=0;
else
    B=varargin{1};
    u=varargin{2};
end
hx=(A*hx0+B*u)';
P=A*P0*A'+Q;
P=0.5*(P+P'); %compensate numeric issues
end

