function [Xi,hx,P] = PF_prediction(Xi,f,hx0,P0,Q,varargin)
[Np,n]=size(Xi);
for i=1:Np
    wi=mvnrnd(zeros(n,1),Q);
    Xi(i,:)=f(Xi(i,:))+wi;
end
hx=mean(Xi);
errX=Xi-hx;
P=1/(Np-1)*(errX'*errX);
P=0.5*(P+P'); %compensate numeric issues
end

