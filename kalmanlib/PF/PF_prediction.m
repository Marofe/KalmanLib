function [Xi,hx,P] = PF_prediction(f,Xi,u,Q,varargin)
[n,Np]=size(Xi);
for i=1:Np
    wi=mvnrnd(zeros(n,1),Q)';
    Xi(:,i)=f(Xi(:,i),u)+wi;
end
hx=mean(Xi,2);
errX=Xi-hx;
P=1/(Np-1)*(errX*errX');
P=0.5*(P+P'); %compensate numeric issues
end

